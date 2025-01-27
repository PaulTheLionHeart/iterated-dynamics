/* (C) 1990, Mark C. Peterson, CompuServe [70441,3353]
   All rights reserved.

   Code may be used in any program provided the author is credited
    either during program execution or in the documentation.  Source
    code may be distributed only in combination with public domain or
    shareware source code.  Source code may be modified provided the
    copyright notice and this message is left unchanged and all
    modifications are clearly documented.

    I would appreciate a copy of any work which incorporates this code,
    however this is optional.

    Mark C. Peterson
    405-C Queen St. Suite #181
    Southington, CT 06489
    (203) 276-9721
*/
#include "port.h"
#include "prototyp.h" // for stricmp

#include "parser.h"

#include "calcfrac.h"
#include "cmdfiles.h"
#include "convert_center_mag.h"
#include "debug_flags.h"
#include "drivers.h"
#include "file_item.h"
#include "fpu087.h"
#include "fractalp.h"
#include "fractals.h"
#include "id_data.h"
#include "jiim.h"
#include "mpmath_c.h"
#include "newton.h"
#include "pixel_grid.h"
#include "sign.h"
#include "stop_msg.h"
#include "trig_fns.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iterator>
#include <stdexcept>
#include <string>

// ** Formula Declarations **
enum MATH_TYPE { D_MATH, M_MATH, L_MATH };

using Function = void();
using FunctionPtr = Function *;

enum
{
    MAX_OPS = 250,
    MAX_ARGS = 100
};

unsigned int g_max_function_ops  = MAX_OPS;
unsigned int g_max_function_args = MAX_ARGS;

static MATH_TYPE s_math_type = D_MATH;
static unsigned long s_num_ops{};
static unsigned long s_num_loads{};
static unsigned long s_num_stores{};
static unsigned long s_num_jumps{};

struct PEND_OP
{
    FunctionPtr f;
    int p;
};

#define MAX_JUMPS 200  // size of JUMP_CONTROL array

struct JUMP_PTRS_ST
{
    int      JumpOpPtr;
    int      JumpLodPtr;
    int      JumpStoPtr;
};

enum class jump_control_type
{
    NONE = 0,
    IF = 1,
    ELSE_IF = 2,
    ELSE = 3,
    END_IF = 4
};

struct JUMP_CONTROL_ST
{
    jump_control_type type;
    JUMP_PTRS_ST ptrs;
    int DestJumpIndex;
};

std::vector<JUMP_CONTROL_ST> g_jump_control;
int g_jump_index{};
static int s_init_jump_index{};

inline void push_jump(jump_control_type type)
{
    JUMP_CONTROL_ST value{};
    value.type = type;
    g_jump_control.push_back(value);
    ++g_jump_index;
}

static bool frm_prescan(std::FILE * open_file);

#define CASE_TERMINATOR case',':\
                        case '\n':\
                        case '(':\
                        case ')':\
                        case '!':\
                        case '=':\
                        case '<':\
                        case '>':\
                        case '|':\
                        case '&':\
                        case '}':\
                        case ':':\
                        case '+':\
                        case '-':\
                        case '*':\
                        case '/':\
                        case '^'

#define CASE_ALPHA      case 'a':\
                        case 'b':\
                        case 'c':\
                        case 'd':\
                        case 'e':\
                        case 'f':\
                        case 'g':\
                        case 'h':\
                        case 'i':\
                        case 'j':\
                        case 'k':\
                        case 'l':\
                        case 'm':\
                        case 'n':\
                        case 'o':\
                        case 'p':\
                        case 'q':\
                        case 'r':\
                        case 's':\
                        case 't':\
                        case 'u':\
                        case 'v':\
                        case 'w':\
                        case 'x':\
                        case 'y':\
                        case 'z'

#define CASE_NUM        case '0':\
                        case '1':\
                        case '2':\
                        case '3':\
                        case '4':\
                        case '5':\
                        case '6':\
                        case '7':\
                        case '8':\
                        case '9'

// token_type definitions
enum class token_type
{
    NOT_A_TOKEN = 0,
    PARENS = 1,
    PARAM_VARIABLE = 2,
    USER_NAMED_VARIABLE = 3,
    PREDEFINED_VARIABLE = 4,
    REAL_CONSTANT = 5,
    COMPLEX_CONSTANT = 6,
    FUNCTION = 7,
    PARAM_FUNCTION = 8,
    FLOW_CONTROL = 9,
    OPERATOR = 10,
    END_OF_FORMULA = 11
};
inline int operator+(token_type value)
{
    return static_cast<int>(value);
}

// token IDs
enum class token_id
{
    NONE = 0,                   // no relevant id for token
    END_OF_FILE = 1,            // begin error condition ids
    ILLEGAL_CHARACTER = 2,      //
    ILLEGAL_VARIABLE_NAME = 3,  //
    TOKEN_TOO_LONG = 4,         //
    FUNC_USED_AS_VAR = 5,       //
    JUMP_MISSING_BOOLEAN = 6,   //
    JUMP_WITH_ILLEGAL_CHAR = 7, //
    UNDEFINED_FUNCTION = 8,     //
    ILLEGAL_OPERATOR = 9,       //
    ILL_FORMED_CONSTANT = 10,   // end error condition ids
    OPEN_PARENS = 1,            //
    CLOSE_PARENS = -1,          //
    FUNC_SIN = 0,               // these ids must match up with the entries in s_func_list
    FUNC_SINH = 1,              //
    FUNC_COS = 2,               //
    FUNC_COSH = 3,              //
    FUNC_SQR = 4,               //
    FUNC_LOG = 5,               //
    FUNC_EXP = 6,               //
    FUNC_ABS = 7,               //
    FUNC_CONJ = 8,              //
    FUNC_REAL = 9,              //
    FUNC_IMAG = 10,             //
    FUNC_FN1 = 11,              //
    FUNC_FN2 = 12,              //
    FUNC_FN3 = 13,              //
    FUNC_FN4 = 14,              //
    FUNC_FLIP = 15,             //
    FUNC_TAN = 16,              //
    FUNC_TANH = 17,             //
    FUNC_COTAN = 18,            //
    FUNC_COTANH = 19,           //
    FUNC_COSXX = 20,            //
    FUNC_SRAND = 21,            //
    FUNC_ASIN = 22,             //
    FUNC_ASINH = 23,            //
    FUNC_ACOS = 24,             //
    FUNC_ACOSH = 25,            //
    FUNC_ATAN = 26,             //
    FUNC_ATANH = 27,            //
    FUNC_SQRT = 28,             //
    FUNC_CABS = 29,             //
    FUNC_FLOOR = 30,            //
    FUNC_CEIL = 31,             //
    FUNC_TRUNC = 32,            //
    FUNC_ROUND = 33,            // end of function name ids
    OP_COMMA = 0,               // these ids must match up with the entries in s_op_list
    OP_NOT_EQUAL = 1,           //
    OP_ASSIGN = 2,              //
    OP_EQUAL = 3,               //
    OP_LT = 4,                  //
    OP_LE = 5,                  //
    OP_GT = 6,                  //
    OP_GE = 7,                  //
    OP_MODULUS = 8,             //
    OP_OR = 9,                  //
    OP_AND = 10,                //
    OP_COLON = 11,              //
    OP_PLUS = 12,               //
    OP_MINUS = 13,              //
    OP_MULTIPLY = 14,           //
    OP_DIVIDE = 15,             //
    OP_POWER = 16,              // end of operator name ids
    VAR_PIXEL = 0,              // these ids must match up with the entries in s_variables
    VAR_P1 = 1,                 //
    VAR_P2 = 2,                 //
    VAR_Z = 3,                  //
    VAR_LAST_SQR = 4,           //
    VAR_PI = 5,                 //
    VAR_E = 6,                  //
    VAR_RAND = 7,               //
    VAR_P3 = 8,                 //
    VAR_WHITE_SQUARE = 9,       //
    VAR_SCREEN_PIX = 10,        //
    VAR_SCREEN_MAX = 11,        //
    VAR_MAX_ITER = 12,          //
    VAR_IS_MANDEL = 13,         //
    VAR_CENTER = 14,            //
    VAR_MAG_X_MAG = 15,         //
    VAR_ROTATION_SKEW = 16,     //
    VAR_P4 = 17,                //
    VAR_P5 = 18,                // end of predefined variable name ids
    JUMP_IF = 1,                // these ids must match up with the entries in s_jump_list, offset by one
    JUMP_ELSE_IF = 2,            //
    JUMP_ELSE = 3,              //
    JUMP_END_IF = 4,             // end of jump list name ids
};
inline int operator+(token_id value)
{
    return static_cast<int>(value);
}

struct token_st
{
    char str[80];
    token_type type;
    token_id id;
    DComplex constant;
};


// MAX_STORES must be even to make Unix alignment work

#define MAX_STORES ((g_max_function_ops/4)*2)  // at most only half the ops can be stores
#define MAX_LOADS  ((unsigned)(g_max_function_ops*.8))  // and 80% can be loads

static std::vector<PEND_OP> s_op;

static void parser_allocate();

Arg *Arg1{};
Arg *Arg2{};

// Some of these variables should be renamed for safety
std::array<Arg, 20> g_stack{};
std::vector<Arg *> Store;
std::vector<Arg *> Load;
int g_op_ptr{};
std::vector<FunctionPtr> f;
std::vector<ConstArg> v;
int g_store_index{};
int g_load_index{};
bool g_is_mandelbrot{true};
unsigned int g_operation_index{};
unsigned int g_variable_index{};
unsigned int g_last_op{};
static unsigned int s_n{};
static unsigned int s_next_op{};
static unsigned int s_init_n{};
static int s_paren{};
static bool s_expecting_arg{};
int InitLodPtr{};
int InitStoPtr{};
int InitOpPtr{};
int g_last_init_op{};
static int s_delta16{};
double g_fudge_limit{};
static double s_fudge{};
static int s_shift_back{};
static bool s_set_random{};
static bool s_randomized{};
static unsigned long s_rand_num{};
bool g_frm_uses_p1{};
bool g_frm_uses_p2{};
bool g_frm_uses_p3{};
bool g_frm_uses_p4{};
bool g_frm_uses_p5{};
bool g_uses_jump{};
bool g_frm_uses_ismand{};
static unsigned int s_chars_in_formula{};

inline bool check_denom(long denom)
{
    if (denom == 0 || g_overflow)
    {
        g_overflow = true;
        return true;
    }
    if (denom == 0)
    {
        return true;
    }
    return false;
}

inline bool check_denom(double denom)
{
    if (std::fabs(denom) <= DBL_MIN)
    {
        g_overflow = true;
        return true;
    }
    return false;
}

#define LastSqr v[4].a

/* ParseErrs() defines; all calls to ParseErrs(), or any variable which will
   be used as the argument in a call to ParseErrs(), should use one of these
   defines.
*/
enum
{
    PE_NO_ERRORS_FOUND = -1,
    PE_SHOULD_BE_ARGUMENT = 0,
    PE_SHOULD_BE_OPERATOR = 1,
    PE_NEED_A_MATCHING_OPEN_PARENS = 2,
    PE_NEED_MORE_CLOSE_PARENS = 3,
    PE_UNDEFINED_OPERATOR = 4,
    PE_UNDEFINED_FUNCTION = 5,
    PE_TABLE_OVERFLOW = 6,
    PE_NO_MATCH_RIGHT_PAREN = 7,
    PE_NO_LEFT_BRACKET_FIRST_LINE = 8,
    PE_UNEXPECTED_EOF = 9,
    PE_INVALID_SYM_USING_NOSYM = 10,
    PE_FORMULA_TOO_LARGE = 11,
    PE_INSUFFICIENT_MEM_FOR_TYPE_FORMULA = 12,
    PE_COULD_NOT_OPEN_FILE_WHERE_FORMULA_LOCATED = 13,
    PE_JUMP_NOT_FIRST = 14,
    PE_NO_CHAR_AFTER_THIS_JUMP = 15,
    PE_JUMP_NEEDS_BOOLEAN = 16,
    PE_ENDIF_REQUIRED_AFTER_ELSE = 17,
    PE_ENDIF_WITH_NO_IF = 18,
    PE_MISPLACED_ELSE_OR_ELSEIF = 19,
    PE_UNMATCHED_IF_IN_INIT_SECTION = 20,
    PE_IF_WITH_NO_ENDIF = 21,
    PE_ERROR_IN_PARSING_JUMP_STATEMENTS = 22,
    PE_TOO_MANY_JUMPS = 23,
    PE_FORMULA_NAME_TOO_LARGE = 24,
    PE_ILLEGAL_ASSIGNMENT = 25,
    PE_ILLEGAL_VAR_NAME = 26,
    PE_INVALID_CONST = 27,
    PE_ILLEGAL_CHAR = 28,
    PE_NESTING_TO_DEEP = 29,
    PE_UNMATCHED_MODULUS = 30,
    PE_FUNC_USED_AS_VAR = 31,
    PE_NO_NEG_AFTER_EXPONENT = 32,
    PE_TOKEN_TOO_LONG = 33,
    PE_SECOND_COLON = 34,
    PE_INVALID_CALL_TO_PARSEERRS = 35
};

static char const *ParseErrs(int which)
{
    static constexpr const char *const ErrStrings[] = {
        "Should be an Argument",
        "Should be an Operator",
        "')' needs a matching '('",
        "Need more ')'",
        "Undefined Operator",
        "Undefined Function",
        "Table overflow",
        "Didn't find matching ')' in symmetry declaration",
        "No '{' found on first line",
        "Unexpected EOF!",
        "Symmetry below is invalid, will use NOSYM",
        "Formula is too large",
        "Insufficient memory to run fractal type 'formula'",
        "Could not open file where formula located",
        "No characters may precede jump instruction",
        "No characters may follow this jump instruction",
        "Jump instruction missing required (boolean argument)",
        "Next jump after \"else\" must be \"endif\"",
        R"("endif" has no matching "if")",
        R"msg(Misplaced "else" or "elseif()")msg",
        R"msg("if ()" in initialization has no matching "endif")msg",
        R"msg("if ()" has no matching "endif")msg",
        "Error in parsing jump statements",
        "Formula has too many jump commands",
        "Formula name has too many characters",
        "Only variables are allowed to left of assignment",
        "Illegal variable name",
        "Invalid constant expression",
        "This character not supported by parser",
        "Nesting of parentheses exceeds maximum depth",
        R"(Unmatched modulus operator "|" in this expression)", //last one
        "Can't use function name as variable",
        "Negative exponent must be enclosed in parens",
        "Variable or constant exceeds 32 character limit",
        R"(Only one ":" permitted in a formula)",
        "Invalid ParseErrs code",
    };
    constexpr int lasterr = std::size(ErrStrings) - 1;
    if (which > lasterr)
    {
        which = lasterr;
    }
    return ErrStrings[which];
}

/* use the following when only float functions are implemented to
   get MP math and Integer math */

static void mStkFunct(FunctionPtr fct)   // call mStk via dStk
{
    Arg1->d = MPC2cmplx(Arg1->m);
    (*fct)();
    Arg1->m = cmplx2MPC(Arg1->d);
}

static void lStkFunct(FunctionPtr fct)   // call lStk via dStk
{
    double y;
    /*
       intermediate variable needed for safety because of
       different size of double and long in Arg union
    */
    y = (double)Arg1->l.y / s_fudge;
    Arg1->d.x = (double)Arg1->l.x / s_fudge;
    Arg1->d.y = y;
    (*fct)();
    if (std::fabs(Arg1->d.x) < g_fudge_limit && std::fabs(Arg1->d.y) < g_fudge_limit)
    {
        Arg1->l.x = (long)(Arg1->d.x * s_fudge);
        Arg1->l.y = (long)(Arg1->d.y * s_fudge);
    }
    else
    {
        g_overflow = true;
    }
}

static unsigned long NewRandNum()
{
    return s_rand_num = ((s_rand_num << 15) + rand15()) ^ s_rand_num;
}

static void lRandom()
{
    v[7].a.l.x = NewRandNum() >> (32 - g_bit_shift);
    v[7].a.l.y = NewRandNum() >> (32 - g_bit_shift);
}

static void dRandom()
{
    long x;
    long y;

    /* Use the same algorithm as for fixed math so that they will generate
           the same fractals when the srand() function is used. */
    x = NewRandNum() >> (32 - g_bit_shift);
    y = NewRandNum() >> (32 - g_bit_shift);
    v[7].a.d.x = ((double)x / (1L << g_bit_shift));
    v[7].a.d.y = ((double)y / (1L << g_bit_shift));

}

static void mRandom()
{
    long x;
    long y;

    /* Use the same algorithm as for fixed math so that they will generate
       the same fractals when the srand() function is used. */
    x = NewRandNum() >> (32 - g_bit_shift);
    y = NewRandNum() >> (32 - g_bit_shift);
    v[7].a.m.x = *fg2MP(x, g_bit_shift);
    v[7].a.m.y = *fg2MP(y, g_bit_shift);
}

static void SetRandFnct()
{
    unsigned Seed;

    if (!s_set_random)
    {
        s_rand_num = Arg1->l.x ^ Arg1->l.y;
    }

    Seed = (unsigned)s_rand_num ^ (unsigned)(s_rand_num >> 16);
    srand(Seed);
    s_set_random = true;

    // Clear out the seed
    NewRandNum();
    NewRandNum();
    NewRandNum();
}

static void RandomSeed()
{
    std::time_t ltime;

    // Use the current time to randomize the random number sequence.
    std::time(&ltime);
    srand((unsigned int)ltime);

    NewRandNum();
    NewRandNum();
    NewRandNum();
    s_randomized = true;
}

static void lStkSRand()
{
    SetRandFnct();
    lRandom();
    Arg1->l = v[7].a.l;
}

static void mStkSRand()
{
    Arg1->l.x = Arg1->m.x.Mant ^ (long)Arg1->m.x.Exp;
    Arg1->l.y = Arg1->m.y.Mant ^ (long)Arg1->m.y.Exp;
    SetRandFnct();
    mRandom();
    Arg1->m = v[7].a.m;
}

static void dStkSRand()
{
    Arg1->l.x = (long)(Arg1->d.x * (1L << g_bit_shift));
    Arg1->l.y = (long)(Arg1->d.y * (1L << g_bit_shift));
    SetRandFnct();
    dRandom();
    Arg1->d = v[7].a.d;
}

static FunctionPtr StkSRand{dStkSRand};

void dStkLodDup()
{
    Arg1 += 2;
    Arg2 += 2;
    *Arg1 = *Load[g_load_index];
    *Arg2 = *Arg1;
    g_load_index += 2;
}

void dStkLodSqr()
{
    Arg1++;
    Arg2++;
    Arg1->d.y = Load[g_load_index]->d.x * Load[g_load_index]->d.y * 2.0;
    Arg1->d.x = (Load[g_load_index]->d.x * Load[g_load_index]->d.x) - (Load[g_load_index]->d.y * Load[g_load_index]->d.y);
    g_load_index++;
}

void dStkLodSqr2()
{
    Arg1++;
    Arg2++;
    LastSqr.d.x = Load[g_load_index]->d.x * Load[g_load_index]->d.x;
    LastSqr.d.y = Load[g_load_index]->d.y * Load[g_load_index]->d.y;
    Arg1->d.y = Load[g_load_index]->d.x * Load[g_load_index]->d.y * 2.0;
    Arg1->d.x = LastSqr.d.x - LastSqr.d.y;
    LastSqr.d.x += LastSqr.d.y;
    LastSqr.d.y = 0;
    g_load_index++;
}

void dStkLodDbl()
{
    Arg1++;
    Arg2++;
    Arg1->d.x = Load[g_load_index]->d.x * 2.0;
    Arg1->d.y = Load[g_load_index]->d.y * 2.0;
    g_load_index++;
}

void dStkSqr0()
{
    LastSqr.d.y = Arg1->d.y * Arg1->d.y; // use LastSqr as temp storage
    Arg1->d.y = Arg1->d.x * Arg1->d.y * 2.0;
    Arg1->d.x = Arg1->d.x * Arg1->d.x - LastSqr.d.y;
}

void dStkSqr3()
{
    Arg1->d.x = Arg1->d.x * Arg1->d.x;
}

void dStkAbs()
{
    Arg1->d.x = std::fabs(Arg1->d.x);
    Arg1->d.y = std::fabs(Arg1->d.y);
}

void mStkAbs()
{
    if (Arg1->m.x.Exp < 0)
    {
        Arg1->m.x.Exp = -Arg1->m.x.Exp;
    }
    if (Arg1->m.y.Exp < 0)
    {
        Arg1->m.y.Exp = -Arg1->m.y.Exp;
    }
}

void lStkAbs()
{
    Arg1->l.x = labs(Arg1->l.x);
    Arg1->l.y = labs(Arg1->l.y);
}

static FunctionPtr StkAbs{dStkAbs};

void dStkSqr()
{
    LastSqr.d.x = Arg1->d.x * Arg1->d.x;
    LastSqr.d.y = Arg1->d.y * Arg1->d.y;
    Arg1->d.y = Arg1->d.x * Arg1->d.y * 2.0;
    Arg1->d.x = LastSqr.d.x - LastSqr.d.y;
    LastSqr.d.x += LastSqr.d.y;
    LastSqr.d.y = 0;
}

void mStkSqr()
{
    LastSqr.m.x = *MPmul(Arg1->m.x, Arg1->m.x);
    LastSqr.m.y = *MPmul(Arg1->m.y, Arg1->m.y);
    Arg1->m.y = *MPmul(Arg1->m.x, Arg1->m.y);
    Arg1->m.y.Exp++;
    Arg1->m.x = *MPsub(LastSqr.m.x, LastSqr.m.y);
    LastSqr.m.x = *MPadd(LastSqr.m.x, LastSqr.m.y);
    LastSqr.m.y.Exp = 0;
    LastSqr.m.y.Mant = 0;
}

void lStkSqr()
{
    LastSqr.l.x = multiply(Arg1->l.x, Arg1->l.x, g_bit_shift);
    LastSqr.l.y = multiply(Arg1->l.y, Arg1->l.y, g_bit_shift);
    Arg1->l.y = multiply(Arg1->l.x, Arg1->l.y, g_bit_shift) << 1;
    Arg1->l.x = LastSqr.l.x - LastSqr.l.y;
    LastSqr.l.x += LastSqr.l.y;
    LastSqr.l.y = 0L;
}

static FunctionPtr StkSqr{dStkSqr};

static void dStkAdd()
{
    Arg2->d.x += Arg1->d.x;
    Arg2->d.y += Arg1->d.y;
    Arg1--;
    Arg2--;
}

static void mStkAdd()
{
    Arg2->m = MPCadd(Arg2->m, Arg1->m);
    Arg1--;
    Arg2--;
}

static void lStkAdd()
{
    Arg2->l.x += Arg1->l.x;
    Arg2->l.y += Arg1->l.y;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkAdd{dStkAdd};

static void dStkSub()
{
    Arg2->d.x -= Arg1->d.x;
    Arg2->d.y -= Arg1->d.y;
    Arg1--;
    Arg2--;
}

static void mStkSub()
{
    Arg2->m = MPCsub(Arg2->m, Arg1->m);
    Arg1--;
    Arg2--;
}

static void lStkSub()
{
    Arg2->l.x -= Arg1->l.x;
    Arg2->l.y -= Arg1->l.y;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkSub{dStkSub};

void dStkConj()
{
    Arg1->d.y = -Arg1->d.y;
}

void mStkConj()
{
    Arg1->m.y.Exp ^= 0x8000;
}

void lStkConj()
{
    Arg1->l.y = -Arg1->l.y;
}

static FunctionPtr StkConj{dStkConj};

void dStkFloor()
{
    Arg1->d.x = floor(Arg1->d.x);
    Arg1->d.y = floor(Arg1->d.y);
}

void mStkFloor()
{
    mStkFunct(dStkFloor);   // call lStk via dStk
}

void lStkFloor()
{
    /*
     * Kill fractional part. This operation truncates negative numbers
     * toward negative infinity as desired.
     */
    Arg1->l.x = (Arg1->l.x) >> g_bit_shift;
    Arg1->l.y = (Arg1->l.y) >> g_bit_shift;
    Arg1->l.x = (Arg1->l.x) << g_bit_shift;
    Arg1->l.y = (Arg1->l.y) << g_bit_shift;
}

static FunctionPtr StkFloor{dStkFloor};

void dStkCeil()
{
    Arg1->d.x = ceil(Arg1->d.x);
    Arg1->d.y = ceil(Arg1->d.y);
}

void mStkCeil()
{
    mStkFunct(dStkCeil);   // call lStk via dStk
}

void lStkCeil()
{
    /* the shift operation does the "floor" operation, so we
       negate everything before the operation */
    Arg1->l.x = (-Arg1->l.x) >> g_bit_shift;
    Arg1->l.y = (-Arg1->l.y) >> g_bit_shift;
    Arg1->l.x = -((Arg1->l.x) << g_bit_shift);
    Arg1->l.y = -((Arg1->l.y) << g_bit_shift);
}

static FunctionPtr StkCeil{dStkCeil};

void dStkTrunc()
{
    Arg1->d.x = (int)(Arg1->d.x);
    Arg1->d.y = (int)(Arg1->d.y);
}

void mStkTrunc()
{
    mStkFunct(dStkTrunc);   // call lStk via dStk
}

void lStkTrunc()
{
    /* shifting and shifting back truncates positive numbers,
       so we make the numbers positive */
    int signx;
    int signy;
    signx = sign(Arg1->l.x);
    signy = sign(Arg1->l.y);
    Arg1->l.x = labs(Arg1->l.x);
    Arg1->l.y = labs(Arg1->l.y);
    Arg1->l.x = (Arg1->l.x) >> g_bit_shift;
    Arg1->l.y = (Arg1->l.y) >> g_bit_shift;
    Arg1->l.x = (Arg1->l.x) << g_bit_shift;
    Arg1->l.y = (Arg1->l.y) << g_bit_shift;
    Arg1->l.x = signx*Arg1->l.x;
    Arg1->l.y = signy*Arg1->l.y;
}

static FunctionPtr StkTrunc{dStkTrunc};

void dStkRound()
{
    Arg1->d.x = floor(Arg1->d.x+.5);
    Arg1->d.y = floor(Arg1->d.y+.5);
}

void mStkRound()
{
    mStkFunct(dStkRound);   // call lStk via dStk
}

void lStkRound()
{
    // Add .5 then truncate
    Arg1->l.x += (1L << g_bit_shift_less_1);
    Arg1->l.y += (1L << g_bit_shift_less_1);
    lStkFloor();
}

static FunctionPtr StkRound{dStkRound};

void dStkZero()
{
    Arg1->d.x = 0.0;
    Arg1->d.y = Arg1->d.x;
}

void mStkZero()
{
    Arg1->m.x.Exp = 0;
    Arg1->m.x.Mant = 0;
    Arg1->m.y.Exp = 0;
    Arg1->m.y.Mant = 0;
}

void lStkZero()
{
    Arg1->l.x = 0;
    Arg1->l.y = Arg1->l.x;
}

static FunctionPtr StkZero{dStkZero};

void dStkOne()
{
    Arg1->d.x = 1.0;
    Arg1->d.y = 0.0;
}

void mStkOne()
{
    Arg1->m = g_mpc_one;
}

void lStkOne()
{
    Arg1->l.x = (long) s_fudge;
    Arg1->l.y = 0L;
}

static FunctionPtr StkOne{dStkOne};

static void dStkReal()
{
    Arg1->d.y = 0.0;
}

static void mStkReal()
{
    Arg1->m.y.Exp = 0;
    Arg1->m.y.Mant = 0;
}

static void lStkReal()
{
    Arg1->l.y = 0l;
}

static FunctionPtr StkReal{dStkReal};

static void dStkImag()
{
    Arg1->d.x = Arg1->d.y;
    Arg1->d.y = 0.0;
}

static void mStkImag()
{
    Arg1->m.x = Arg1->m.y;
    Arg1->m.y.Exp = 0;
    Arg1->m.y.Mant = 0;
}

static void lStkImag()
{
    Arg1->l.x = Arg1->l.y;
    Arg1->l.y = 0l;
}

static FunctionPtr StkImag{dStkImag};

static void dStkNeg()
{
    Arg1->d.x = -Arg1->d.x;
    Arg1->d.y = -Arg1->d.y;
}

static void mStkNeg()
{
    Arg1->m.x.Exp ^= 0x8000;
    Arg1->m.y.Exp ^= 0x8000;
}

static void lStkNeg()
{
    Arg1->l.x = -Arg1->l.x;
    Arg1->l.y = -Arg1->l.y;
}

static FunctionPtr StkNeg{dStkNeg};

void dStkMul()
{
    FPUcplxmul(&Arg2->d, &Arg1->d, &Arg2->d);
    Arg1--;
    Arg2--;
}

static void mStkMul()
{
    Arg2->m = MPCmul(Arg2->m, Arg1->m);
    Arg1--;
    Arg2--;
}

void lStkMul()
{
    long x;
    long y;

    x = multiply(Arg2->l.x, Arg1->l.x, g_bit_shift) -
        multiply(Arg2->l.y, Arg1->l.y, g_bit_shift);
    y = multiply(Arg2->l.y, Arg1->l.x, g_bit_shift) +
        multiply(Arg2->l.x, Arg1->l.y, g_bit_shift);
    Arg2->l.x = x;
    Arg2->l.y = y;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkMul{dStkMul};

static void dStkDiv()
{
    FPUcplxdiv(&Arg2->d, &Arg1->d, &Arg2->d);
    Arg1--;
    Arg2--;
}

static void mStkDiv()
{
    Arg2->m = MPCdiv(Arg2->m, Arg1->m);
    Arg1--;
    Arg2--;
}

static void lStkDiv()
{
    long x;
    long y;
    long mod;
    long x2;
    long y2;

    mod = multiply(Arg1->l.x, Arg1->l.x, g_bit_shift) +
          multiply(Arg1->l.y, Arg1->l.y, g_bit_shift);
    x = divide(Arg1->l.x, mod, g_bit_shift);
    y = -divide(Arg1->l.y, mod, g_bit_shift);
    x2 = multiply(Arg2->l.x, x, g_bit_shift) - multiply(Arg2->l.y, y, g_bit_shift);
    y2 = multiply(Arg2->l.y, x, g_bit_shift) + multiply(Arg2->l.x, y, g_bit_shift);
    Arg2->l.x = x2;
    Arg2->l.y = y2;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkDiv{dStkDiv};

static void dStkMod()
{
    Arg1->d.x = (Arg1->d.x * Arg1->d.x) + (Arg1->d.y * Arg1->d.y);
    Arg1->d.y = 0.0;
}

static void mStkMod()
{
    Arg1->m.x = MPCmod(Arg1->m);
    Arg1->m.y.Exp = 0;
    Arg1->m.y.Mant = 0;
}

static void lStkMod()
{
    //   Arg1->l.x = multiply(Arg2->l.x, Arg1->l.x, bitshift) +
    //   multiply(Arg2->l.y, Arg1->l.y, bitshift);
    Arg1->l.x = multiply(Arg1->l.x, Arg1->l.x, g_bit_shift) +
                multiply(Arg1->l.y, Arg1->l.y, g_bit_shift);
    if (Arg1->l.x < 0)
    {
        g_overflow = true;
    }
    Arg1->l.y = 0L;
}

static FunctionPtr StkMod{dStkMod};

static void StkSto()
{
    assert(Store[g_store_index] != nullptr);
    *Store[g_store_index++] = *Arg1;
}

static void StkLod()
{
    Arg1++;
    Arg2++;
    *Arg1 = *Load[g_load_index++];
}

static void StkClr()
{
    g_stack[0] = *Arg1;
    Arg1 = &g_stack[0];
    Arg2 = &g_stack[0];
    Arg2--;
}

void dStkFlip()
{
    double t;

    t = Arg1->d.x;
    Arg1->d.x = Arg1->d.y;
    Arg1->d.y = t;
}

void mStkFlip()
{
    MP t;

    t = Arg1->m.x;
    Arg1->m.x = Arg1->m.y;
    Arg1->m.y = t;
}

void lStkFlip()
{
    long t;

    t = Arg1->l.x;
    Arg1->l.x = Arg1->l.y;
    Arg1->l.y = t;
}

static FunctionPtr StkFlip{dStkFlip};

void dStkSin()
{
    double sinx;
    double cosx;
    double sinhy;
    double coshy;

    FPUsincos(&Arg1->d.x, &sinx, &cosx);
    FPUsinhcosh(&Arg1->d.y, &sinhy, &coshy);
    Arg1->d.x = sinx*coshy;
    Arg1->d.y = cosx*sinhy;
}

void mStkSin()
{
    mStkFunct(dStkSin);   // call lStk via dStk
}

void lStkSin()
{
    long x;
    long y;
    long sinx;
    long cosx;
    long sinhy;
    long coshy;
    x = Arg1->l.x >> s_delta16;
    y = Arg1->l.y >> s_delta16;
    SinCos086(x, &sinx, &cosx);
    SinhCosh086(y, &sinhy, &coshy);
    Arg1->l.x = multiply(sinx, coshy, s_shift_back);
    Arg1->l.y = multiply(cosx, sinhy, s_shift_back);
}

static FunctionPtr StkSin{dStkSin};

/* The following functions are supported by both the parser and for fn
   variable replacement.
*/
void dStkTan()
{
    double sinx;
    double cosx;
    double sinhy;
    double coshy;
    double denom;
    Arg1->d.x *= 2;
    Arg1->d.y *= 2;
    FPUsincos(&Arg1->d.x, &sinx, &cosx);
    FPUsinhcosh(&Arg1->d.y, &sinhy, &coshy);
    denom = cosx + coshy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->d.x = sinx/denom;
    Arg1->d.y = sinhy/denom;
}

void mStkTan()
{
    mStkFunct(dStkTan);   // call lStk via dStk
}

void lStkTan()
{
    long x;
    long y;
    long sinx;
    long cosx;
    long sinhy;
    long coshy;
    long denom;
    x = Arg1->l.x >> s_delta16;
    x = x << 1;
    y = Arg1->l.y >> s_delta16;
    y = y << 1;
    SinCos086(x, &sinx, &cosx);
    SinhCosh086(y, &sinhy, &coshy);
    denom = cosx + coshy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->l.x = divide(sinx, denom, g_bit_shift);
    Arg1->l.y = divide(sinhy, denom, g_bit_shift);
}

static FunctionPtr StkTan{dStkTan};

void dStkTanh()
{
    double siny;
    double cosy;
    double sinhx;
    double coshx;
    double denom;
    Arg1->d.x *= 2;
    Arg1->d.y *= 2;
    FPUsincos(&Arg1->d.y, &siny, &cosy);
    FPUsinhcosh(&Arg1->d.x, &sinhx, &coshx);
    denom = coshx + cosy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->d.x = sinhx/denom;
    Arg1->d.y = siny/denom;
}

void mStkTanh()
{
    mStkFunct(dStkTanh);   // call lStk via dStk
}

void lStkTanh()
{
    long x;
    long y;
    long siny;
    long cosy;
    long sinhx;
    long coshx;
    long denom;
    x = Arg1->l.x >> s_delta16;
    x = x << 1;
    y = Arg1->l.y >> s_delta16;
    y = y << 1;
    SinCos086(y, &siny, &cosy);
    SinhCosh086(x, &sinhx, &coshx);
    denom = coshx + cosy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->l.x = divide(sinhx, denom, g_bit_shift);
    Arg1->l.y = divide(siny, denom, g_bit_shift);
}

static FunctionPtr StkTanh{dStkTanh};

void dStkCoTan()
{
    double sinx;
    double cosx;
    double sinhy;
    double coshy;
    double denom;
    Arg1->d.x *= 2;
    Arg1->d.y *= 2;
    FPUsincos(&Arg1->d.x, &sinx, &cosx);
    FPUsinhcosh(&Arg1->d.y, &sinhy, &coshy);
    denom = coshy - cosx;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->d.x = sinx/denom;
    Arg1->d.y = -sinhy/denom;
}

void mStkCoTan()
{
    mStkFunct(dStkCoTan);   // call lStk via dStk
}

void lStkCoTan()
{
    long x;
    long y;
    long sinx;
    long cosx;
    long sinhy;
    long coshy;
    long denom;
    x = Arg1->l.x >> s_delta16;
    x = x << 1;
    y = Arg1->l.y >> s_delta16;
    y = y << 1;
    SinCos086(x, &sinx, &cosx);
    SinhCosh086(y, &sinhy, &coshy);
    denom = coshy - cosx;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->l.x = divide(sinx, denom, g_bit_shift);
    Arg1->l.y = -divide(sinhy, denom, g_bit_shift);
}

static FunctionPtr StkCoTan{dStkCoTan};

void dStkCoTanh()
{
    double siny;
    double cosy;
    double sinhx;
    double coshx;
    double denom;
    Arg1->d.x *= 2;
    Arg1->d.y *= 2;
    FPUsincos(&Arg1->d.y, &siny, &cosy);
    FPUsinhcosh(&Arg1->d.x, &sinhx, &coshx);
    denom = coshx - cosy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->d.x = sinhx/denom;
    Arg1->d.y = -siny/denom;
}

void mStkCoTanh()
{
    mStkFunct(dStkCoTanh);   // call lStk via dStk
}

void lStkCoTanh()
{
    long x;
    long y;
    long siny;
    long cosy;
    long sinhx;
    long coshx;
    long denom;
    x = Arg1->l.x >> s_delta16;
    x = x << 1;
    y = Arg1->l.y >> s_delta16;
    y = y << 1;
    SinCos086(y, &siny, &cosy);
    SinhCosh086(x, &sinhx, &coshx);
    denom = coshx - cosy;
    if (check_denom(denom))
    {
        return;
    }
    Arg1->l.x = divide(sinhx, denom, g_bit_shift);
    Arg1->l.y = -divide(siny, denom, g_bit_shift);
}

static FunctionPtr StkCoTanh{dStkCoTanh};

/* The following functions are not directly used by the parser - support
   for the parser was not provided because the existing parser language
   represents these quite easily. They are used for fn variable support
   in miscres.c but are placed here because they follow the pattern of
   the other parser functions.
*/

void dStkRecip()
{
    double mod;
    mod =Arg1->d.x * Arg1->d.x + Arg1->d.y * Arg1->d.y;
    if (check_denom(mod))
    {
        return;
    }
    Arg1->d.x =  Arg1->d.x/mod;
    Arg1->d.y = -Arg1->d.y/mod;
}

void mStkRecip()
{
    MP mod;
    mod = *MPadd(*MPmul(Arg1->m.x, Arg1->m.x), *MPmul(Arg1->m.y, Arg1->m.y));
    if (mod.Mant == 0L)
    {
        g_overflow = true;
        return;
    }
    Arg1->m.x = *MPdiv(Arg1->m.x, mod);
    Arg1->m.y = *MPdiv(Arg1->m.y, mod);
    Arg1->m.y.Exp ^= 0x8000;
}

void lStkRecip()
{
    long mod;
    mod = multiply(Arg1->l.x, Arg1->l.x, g_bit_shift)
          + multiply(Arg1->l.y, Arg1->l.y, g_bit_shift);
    if (check_denom(mod))
    {
        return;
    }
    Arg1->l.x =  divide(Arg1->l.x, mod, g_bit_shift);
    Arg1->l.y = -divide(Arg1->l.y, mod, g_bit_shift);
}

void StkIdent()
{
    // do nothing - the function Z
}

void dStkSinh()
{
    double siny;
    double cosy;
    double sinhx;
    double coshx;

    FPUsincos(&Arg1->d.y, &siny, &cosy);
    FPUsinhcosh(&Arg1->d.x, &sinhx, &coshx);
    Arg1->d.x = sinhx*cosy;
    Arg1->d.y = coshx*siny;
}

void mStkSinh()
{
    mStkFunct(dStkSinh);   // call lStk via dStk
}

void lStkSinh()
{
    long x;
    long y;
    long sinhx;
    long coshx;
    long siny;
    long cosy;

    x = Arg1->l.x >> s_delta16;
    y = Arg1->l.y >> s_delta16;
    SinCos086(y, &siny, &cosy);
    SinhCosh086(x, &sinhx, &coshx);
    Arg1->l.x = multiply(cosy, sinhx, s_shift_back);
    Arg1->l.y = multiply(siny, coshx, s_shift_back);
}

static FunctionPtr StkSinh{dStkSinh};

void dStkCos()
{
    double sinx;
    double cosx;
    double sinhy;
    double coshy;

    FPUsincos(&Arg1->d.x, &sinx, &cosx);
    FPUsinhcosh(&Arg1->d.y, &sinhy, &coshy);
    Arg1->d.x = cosx*coshy;
    Arg1->d.y = -sinx*sinhy;
}

void mStkCos()
{
    mStkFunct(dStkCos);   // call lStk via dStk
}

void lStkCos()
{
    long x;
    long y;
    long sinx;
    long cosx;
    long sinhy;
    long coshy;

    x = Arg1->l.x >> s_delta16;
    y = Arg1->l.y >> s_delta16;
    SinCos086(x, &sinx, &cosx);
    SinhCosh086(y, &sinhy, &coshy);
    Arg1->l.x = multiply(cosx, coshy, s_shift_back);
    Arg1->l.y = -multiply(sinx, sinhy, s_shift_back);
}

static FunctionPtr StkCos{dStkCos};

// Bogus version of cos, to replicate bug which was in regular cos till v16:

void dStkCosXX()
{
    dStkCos();
    Arg1->d.y = -Arg1->d.y;
}

void mStkCosXX()
{
    mStkFunct(dStkCosXX);   // call lStk via dStk
}

void lStkCosXX()
{
    lStkCos();
    Arg1->l.y = -Arg1->l.y;
}

static FunctionPtr StkCosXX{dStkCosXX};

void dStkCosh()
{
    double siny;
    double cosy;
    double sinhx;
    double coshx;

    FPUsincos(&Arg1->d.y, &siny, &cosy);
    FPUsinhcosh(&Arg1->d.x, &sinhx, &coshx);
    Arg1->d.x = coshx*cosy;
    Arg1->d.y = sinhx*siny;
}

void mStkCosh()
{
    mStkFunct(dStkCosh);   // call lStk via dStk
}

void lStkCosh()
{
    long x;
    long y;
    long sinhx;
    long coshx;
    long siny;
    long cosy;

    x = Arg1->l.x >> s_delta16;
    y = Arg1->l.y >> s_delta16;
    SinCos086(y, &siny, &cosy);
    SinhCosh086(x, &sinhx, &coshx);
    Arg1->l.x = multiply(cosy, coshx, s_shift_back);
    Arg1->l.y = multiply(siny, sinhx, s_shift_back);
}

static FunctionPtr StkCosh{dStkCosh};

void dStkASin()
{
    Arcsinz(Arg1->d, &(Arg1->d));
}

void mStkASin()
{
    mStkFunct(dStkASin);
}

void lStkASin()
{
    lStkFunct(dStkASin);
}

static FunctionPtr StkASin{dStkASin};

void dStkASinh()
{
    Arcsinhz(Arg1->d, &(Arg1->d));
}

void mStkASinh()
{
    mStkFunct(dStkASinh);
}

void lStkASinh()
{
    lStkFunct(dStkASinh);
}

static FunctionPtr StkASinh{dStkASinh};

void dStkACos()
{
    Arccosz(Arg1->d, &(Arg1->d));
}

void mStkACos()
{
    mStkFunct(dStkACos);
}

void lStkACos()
{
    lStkFunct(dStkACos);
}

static FunctionPtr StkACos{dStkACos};

void dStkACosh()
{
    Arccoshz(Arg1->d, &(Arg1->d));
}

void mStkACosh()
{
    mStkFunct(dStkACosh);
}

void lStkACosh()
{
    lStkFunct(dStkACosh);
}

static FunctionPtr StkACosh{dStkACosh};

void dStkATan()
{
    Arctanz(Arg1->d, &(Arg1->d));
}

void mStkATan()
{
    mStkFunct(dStkATan);
}

void lStkATan()
{
    lStkFunct(dStkATan);
}

static FunctionPtr StkATan{dStkATan};

void dStkATanh()
{
    Arctanhz(Arg1->d, &(Arg1->d));
}

void mStkATanh()
{
    mStkFunct(dStkATanh);
}

void lStkATanh()
{
    lStkFunct(dStkATanh);
}

static FunctionPtr StkATanh{dStkATanh};

void dStkSqrt()
{
    Arg1->d = ComplexSqrtFloat(Arg1->d.x, Arg1->d.y);
}

void mStkSqrt()
{
    mStkFunct(dStkSqrt);
}

void lStkSqrt()
{
    // lStkFunct(dStkSqrt);
    Arg1->l = ComplexSqrtLong(Arg1->l.x, Arg1->l.y);
}

static FunctionPtr StkSqrt{dStkSqrt};

void dStkCAbs()
{
    Arg1->d.x = std::sqrt(sqr(Arg1->d.x)+sqr(Arg1->d.y));
    Arg1->d.y = 0.0;
}

void mStkCAbs()
{
    mStkFunct(dStkCAbs);
}

void lStkCAbs()
{
    lStkFunct(dStkCAbs);
}

static FunctionPtr StkCAbs{dStkCAbs};

static void dStkLT()
{
    Arg2->d.x = (double)(Arg2->d.x < Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkLT()
{
    Arg2->m.x = *fg2MP((long)(MPcmp(Arg2->m.x, Arg1->m.x) == -1), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkLT()
{
    Arg2->l.x = (long)(Arg2->l.x < Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkLT{dStkLT};

static void dStkGT()
{
    Arg2->d.x = (double)(Arg2->d.x > Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkGT()
{
    Arg2->m.x = *fg2MP((long)(MPcmp(Arg2->m.x, Arg1->m.x) == 1), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkGT()
{
    Arg2->l.x = (long)(Arg2->l.x > Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkGT{dStkGT};

static void dStkLTE()
{
    Arg2->d.x = (double)(Arg2->d.x <= Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkLTE()
{
    int comp;

    comp = MPcmp(Arg2->m.x, Arg1->m.x);
    Arg2->m.x = *fg2MP((long)(comp == -1 || comp == 0), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkLTE()
{
    Arg2->l.x = (long)(Arg2->l.x <= Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkLTE{dStkLTE};

static void dStkGTE()
{
    Arg2->d.x = (double)(Arg2->d.x >= Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkGTE()
{
    int comp;

    comp = MPcmp(Arg2->m.x, Arg1->m.x);
    Arg2->m.x = *fg2MP((long)(comp == 1 || comp == 0), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkGTE()
{
    Arg2->l.x = (long)(Arg2->l.x >= Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkGTE{dStkGTE};

static void dStkEQ()
{
    Arg2->d.x = (double)(Arg2->d.x == Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkEQ()
{
    int comp;

    comp = MPcmp(Arg2->m.x, Arg1->m.x);
    Arg2->m.x = *fg2MP((long)(comp == 0), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkEQ()
{
    Arg2->l.x = (long)(Arg2->l.x == Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkEQ{dStkEQ};

static void dStkNE()
{
    Arg2->d.x = (double)(Arg2->d.x != Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkNE()
{
    int comp;

    comp = MPcmp(Arg2->m.x, Arg1->m.x);
    Arg2->m.x = *fg2MP((long)(comp != 0), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkNE()
{
    Arg2->l.x = (long)(Arg2->l.x != Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkNE{dStkNE};

static void dStkOR()
{
    Arg2->d.x = (double)(Arg2->d.x || Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkOR()
{
    Arg2->m.x = *fg2MP((long)(Arg2->m.x.Mant || Arg1->m.x.Mant), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkOR()
{
    Arg2->l.x = (long)(Arg2->l.x || Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkOR{dStkOR};

static void dStkAND()
{
    Arg2->d.x = (double)(Arg2->d.x && Arg1->d.x);
    Arg2->d.y = 0.0;
    Arg1--;
    Arg2--;
}

static void mStkAND()
{
    Arg2->m.x = *fg2MP((long)(Arg2->m.x.Mant && Arg1->m.x.Mant), 0);
    Arg2->m.y.Exp = 0;
    Arg2->m.y.Mant = 0;
    Arg1--;
    Arg2--;
}

static void lStkAND()
{
    Arg2->l.x = (long)(Arg2->l.x && Arg1->l.x) << g_bit_shift;
    Arg2->l.y = 0l;
    Arg1--;
    Arg2--;
}

static FunctionPtr StkAND{dStkAND};

void dStkLog()
{
    FPUcplxlog(&Arg1->d, &Arg1->d);
}

void mStkLog()
{
    mStkFunct(dStkLog);   // call lStk via dStk
}

void lStkLog()
{
    lStkFunct(dStkLog);
}

static FunctionPtr StkLog{dStkLog};

void dStkExp()
{
    FPUcplxexp(&Arg1->d, &Arg1->d);
}

void mStkExp()
{
    mStkFunct(dStkExp);   // call lStk via dStk
}

void lStkExp()
{
    lStkFunct(dStkExp);
}

static FunctionPtr StkExp{dStkExp};

void dStkPwr()
{
    Arg2->d = ComplexPower(Arg2->d, Arg1->d);
    Arg1--;
    Arg2--;
}

void mStkPwr()
{
    DComplex x;
    DComplex y;

    x = MPC2cmplx(Arg2->m);
    y = MPC2cmplx(Arg1->m);
    x = ComplexPower(x, y);
    Arg2->m = cmplx2MPC(x);
    Arg1--;
    Arg2--;
}

void lStkPwr()
{
    DComplex x;
    DComplex y;

    x.x = (double)Arg2->l.x / s_fudge;
    x.y = (double)Arg2->l.y / s_fudge;
    y.x = (double)Arg1->l.x / s_fudge;
    y.y = (double)Arg1->l.y / s_fudge;
    x = ComplexPower(x, y);
    if (std::fabs(x.x) < g_fudge_limit && std::fabs(x.y) < g_fudge_limit)
    {
        Arg2->l.x = (long)(x.x * s_fudge);
        Arg2->l.y = (long)(x.y * s_fudge);
    }
    else
    {
        g_overflow = true;
    }
    Arg1--;
    Arg2--;
}

static FunctionPtr StkPwr{dStkPwr};

static void EndInit()
{
    g_last_init_op = g_op_ptr;
    s_init_jump_index = g_jump_index;
}

static FunctionPtr PtrEndInit{EndInit};

void StkJump()
{
    g_op_ptr =  g_jump_control[g_jump_index].ptrs.JumpOpPtr;
    g_load_index = g_jump_control[g_jump_index].ptrs.JumpLodPtr;
    g_store_index = g_jump_control[g_jump_index].ptrs.JumpStoPtr;
    g_jump_index = g_jump_control[g_jump_index].DestJumpIndex;
}

void dStkJumpOnFalse()
{
    if (Arg1->d.x == 0)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

void mStkJumpOnFalse()
{
    if (Arg1->m.x.Mant == 0)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

void lStkJumpOnFalse()
{
    if (Arg1->l.x == 0)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

static FunctionPtr StkJumpOnFalse{dStkJumpOnFalse};

void dStkJumpOnTrue()
{
    if (Arg1->d.x)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

void mStkJumpOnTrue()
{
    if (Arg1->m.x.Mant)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

void lStkJumpOnTrue()
{
    if (Arg1->l.x)
    {
        StkJump();
    }
    else
    {
        g_jump_index++;
    }
}

static FunctionPtr StkJumpOnTrue{dStkJumpOnTrue};

void StkJumpLabel()
{
    g_jump_index++;
}

static unsigned int SkipWhiteSpace(char const *Str)
{
    unsigned n;
    bool Done = false;
    for (n = 0; !Done; n++)
    {
        switch (Str[n])
        {
        case ' ':
        case '\t':
        case '\n':
        case '\r':
            break;
        default:
            Done = true;
        }
    }
    return n - 1;
}

// detect if constant is part of a (a,b) construct
static bool isconst_pair(char const *Str)
{
    int n;
    int j;
    bool answer = false;
    // skip past first number
    for (n = 0; std::isdigit(Str[n]) || Str[n] == '.'; n++)
    {
    }
    if (Str[n] == ',')
    {
        j = n + SkipWhiteSpace(&Str[n+1]) + 1;
        if (std::isdigit(Str[j])
            || (Str[j] == '-' && (std::isdigit(Str[j+1]) || Str[j+1] == '.'))
            || Str[j] == '.')
        {
            answer = true;
        }
    }
    return answer;
}

static ConstArg *is_const(char const *Str, int Len)
{
    DComplex z;
    // next line enforces variable vs constant naming convention
    for (unsigned n = 0U; n < g_variable_index; n++)
    {
        if (v[n].len == Len)
        {
            if (!strnicmp(v[n].s, Str, Len))
            {
                if (n == 1)          // The formula uses 'p1'.
                {
                    g_frm_uses_p1 = true;
                }
                if (n == 2)          // The formula uses 'p2'.
                {
                    g_frm_uses_p2 = true;
                }
                if (n == 7)          // The formula uses 'rand'.
                {
                    RandomSeed();
                }
                if (n == 8)          // The formula uses 'p3'.
                {
                    g_frm_uses_p3 = true;
                }
                if (n == 13)          // The formula uses 'ismand'.
                {
                    g_frm_uses_ismand = true;
                }
                if (n == 17)          // The formula uses 'p4'.
                {
                    g_frm_uses_p4 = true;
                }
                if (n == 18)          // The formula uses 'p5'.
                {
                    g_frm_uses_p5 = true;
                }
                if (n == 10 || n == 11 || n == 12)
                {
                    if (s_math_type == L_MATH)
                    {
                        driver_unget_key('f');
                    }
                }
                if (!isconst_pair(Str))
                {
                    return &v[n];
                }
            }
        }
    }
    v[g_variable_index].s = Str;
    v[g_variable_index].len = Len;
    v[g_variable_index].a.d.y = 0.0;
    v[g_variable_index].a.d.x = v[g_variable_index].a.d.y;

    // v[vsp].a should already be zeroed out
    switch (s_math_type)
    {
    case M_MATH:
        v[g_variable_index].a.m.x.Exp = 0;
        v[g_variable_index].a.m.x.Mant = 0;
        v[g_variable_index].a.m.y.Exp = 0;
        v[g_variable_index].a.m.y.Mant = 0;
        break;
    case L_MATH:
        v[g_variable_index].a.l.y = 0;
        v[g_variable_index].a.l.x = v[g_variable_index].a.l.y;
        break;
    case D_MATH:
        break;
    }

    if (std::isdigit(Str[0])
        || (Str[0] == '-' && (std::isdigit(Str[1]) || Str[1] == '.'))
        || Str[0] == '.')
    {
        assert(g_operation_index > 0);
        assert(g_operation_index == s_op.size());
        if (s_op.back().f == StkNeg)
        {
            s_op.pop_back();
            g_operation_index--;
            Str = Str - 1;
            s_init_n--;
            v[g_variable_index].len++;
        }
        unsigned n;
        for (n = 1; std::isdigit(Str[n]) || Str[n] == '.'; n++)
        {
        }
        if (Str[n] == ',')
        {
            unsigned j = n + SkipWhiteSpace(&Str[n+1]) + 1;
            if (std::isdigit(Str[j])
                || (Str[j] == '-' && (std::isdigit(Str[j+1]) || Str[j+1] == '.'))
                || Str[j] == '.')
            {
                z.y = std::atof(&Str[j]);
                for (; std::isdigit(Str[j]) || Str[j] == '.' || Str[j] == '-'; j++)
                {
                }
                v[g_variable_index].len = j;
            }
            else
            {
                z.y = 0.0;
            }
        }
        else
        {
            z.y = 0.0;
        }
        z.x = std::atof(Str);
        switch (s_math_type)
        {
        case D_MATH:
            v[g_variable_index].a.d = z;
            break;
        case M_MATH:
            v[g_variable_index].a.m = cmplx2MPC(z);
            break;
        case L_MATH:
            v[g_variable_index].a.l.x = (long)(z.x * s_fudge);
            v[g_variable_index].a.l.y = (long)(z.y * s_fudge);
            break;
        }
        v[g_variable_index].s = Str;
    }
    return &v[g_variable_index++];
}

namespace
{

struct FNCT_LIST
{
    char const *s;
    FunctionPtr *ptr;
};

} // namespace

static constexpr std::array<char const *, 4> s_jump_list
{
    "if",
    "elseif",
    "else",
    "endif"
};

static FunctionPtr StkTrig0{dStkSin};
static FunctionPtr StkTrig1{dStkSqr};
static FunctionPtr StkTrig2{dStkSinh};
static FunctionPtr StkTrig3{dStkCosh};

/* return values
    0 - Not a jump
    1 - if
    2 - elseif
    3 - else
    4 - endif
*/
static jump_control_type is_jump(char const *Str, int Len)
{
    for (int i = 0; i < static_cast<int>(s_jump_list.size()); i++)
    {
        if ((int) std::strlen(s_jump_list[i]) == Len && strnicmp(s_jump_list[i], Str, Len) == 0)
        {
            return static_cast<jump_control_type>(i + 1);
        }
    }
    return jump_control_type::NONE;
}

char g_max_function{};

static constexpr std::array<FNCT_LIST, 34> s_func_list
{
    FNCT_LIST{"sin",   &StkSin},
    FNCT_LIST{"sinh",  &StkSinh},
    FNCT_LIST{"cos",   &StkCos},
    FNCT_LIST{"cosh",  &StkCosh},
    FNCT_LIST{"sqr",   &StkSqr},
    FNCT_LIST{"log",   &StkLog},
    FNCT_LIST{"exp",   &StkExp},
    FNCT_LIST{"abs",   &StkAbs},
    FNCT_LIST{"conj",  &StkConj},
    FNCT_LIST{"real",  &StkReal},
    FNCT_LIST{"imag",  &StkImag},
    FNCT_LIST{"fn1",   &StkTrig0},
    FNCT_LIST{"fn2",   &StkTrig1},
    FNCT_LIST{"fn3",   &StkTrig2},
    FNCT_LIST{"fn4",   &StkTrig3},
    FNCT_LIST{"flip",  &StkFlip},
    FNCT_LIST{"tan",   &StkTan},
    FNCT_LIST{"tanh",  &StkTanh},
    FNCT_LIST{"cotan", &StkCoTan},
    FNCT_LIST{"cotanh", &StkCoTanh},
    FNCT_LIST{"cosxx", &StkCosXX},
    FNCT_LIST{"srand", &StkSRand},
    FNCT_LIST{"asin",  &StkASin},
    FNCT_LIST{"asinh", &StkASinh},
    FNCT_LIST{"acos",  &StkACos},
    FNCT_LIST{"acosh", &StkACosh},
    FNCT_LIST{"atan",  &StkATan},
    FNCT_LIST{"atanh", &StkATanh},
    FNCT_LIST{"sqrt",  &StkSqrt},
    FNCT_LIST{"cabs",  &StkCAbs},
    FNCT_LIST{"floor", &StkFloor},
    FNCT_LIST{"ceil",  &StkCeil},
    FNCT_LIST{"trunc", &StkTrunc},
    FNCT_LIST{"round", &StkRound},
};

static std::array<char const *, 17> s_op_list
{
    ",",    //  0
    "!=",   //  1
    "=",    //  2
    "==",   //  3
    "<",    //  4
    "<=",   //  5
    ">",    //  6
    ">=",   //  7
    "|",    //  8
    "||",   //  9
    "&&",   // 10
    ":",    // 11
    "+",    // 12
    "-",    // 13
    "*",    // 14
    "/",    // 15
    "^"     // 16
};

static void NotAFnct()
{
}

static void FnctNotFound()
{
}

// determine if s names a function and if so which one
static int whichfn(char const *s, int len)
{
    int out;
    if (len != 3)
    {
        out = 0;
    }
    else if (strnicmp(s, "fn", 2))
    {
        out = 0;
    }
    else
    {
        out = std::atoi(s+2);
    }
    if (out < 1 || out > 4)
    {
        out = 0;
    }
    return out;
}

static FunctionPtr is_func(char const *Str, int Len)
{
    unsigned n = SkipWhiteSpace(&Str[Len]);
    if (Str[Len+n] == '(')
    {
        for (n = 0; n < static_cast<unsigned>(s_func_list.size()); n++)
        {
            if (!strnicmp(s_func_list[n].s, Str, Len))
            {
                // count function variables
                int functnum = whichfn(Str, Len);
                if (functnum != 0)
                {
                    if (functnum > g_max_function)
                    {
                        g_max_function = (char)functnum;
                    }
                }
                return *s_func_list[n].ptr;
            }
        }
        return FnctNotFound;
    }
    return NotAFnct;
}

static void RecSortPrec()
{
    int ThisOp = s_next_op++;
    while (s_op[ThisOp].p > s_op[s_next_op].p && s_next_op < g_operation_index)
    {
        RecSortPrec();
    }
    if (g_op_ptr > static_cast<int>(f.size()))
    {
        throw std::runtime_error(
            "OpPtr (" + std::to_string(g_op_ptr) + ") exceeds size of f[] (" + std::to_string(f.size()) + ")");
    }
    f.push_back(s_op[ThisOp].f);
    ++g_op_ptr;
}

static constexpr std::array<char const *, 19> s_variables
{
    "pixel",        // v[0]
    "p1",           // v[1]
    "p2",           // v[2]
    "z",            // v[3]
    "LastSqr",      // v[4]
    "pi",           // v[5]
    "e",            // v[6]
    "rand",         // v[7]
    "p3",           // v[8]
    "whitesq",      // v[9]
    "scrnpix",      // v[10]
    "scrnmax",      // v[11]
    "maxit",        // v[12]
    "ismand",       // v[13]
    "center",       // v[14]
    "magxmag",      // v[15]
    "rotskew",      // v[16]
    "p4",           // v[17]
    "p5"            // v[18]
};

struct SymmetryName
{
    char const *s;
    symmetry_type n;
};

static constexpr std::array<SymmetryName, 14> s_symmetry_names
{
    SymmetryName{ "NOSYM",         symmetry_type::NONE },
    SymmetryName{ "XAXIS_NOPARM",  symmetry_type::X_AXIS_NO_PARAM },
    SymmetryName{ "XAXIS",         symmetry_type::X_AXIS },
    SymmetryName{ "YAXIS_NOPARM",  symmetry_type::Y_AXIS_NO_PARAM },
    SymmetryName{ "YAXIS",         symmetry_type::Y_AXIS },
    SymmetryName{ "XYAXIS_NOPARM", symmetry_type::XY_AXIS_NO_PARAM },
    SymmetryName{ "XYAXIS",        symmetry_type::XY_AXIS },
    SymmetryName{ "ORIGIN_NOPARM", symmetry_type::ORIGIN_NO_PARAM },
    SymmetryName{ "ORIGIN",        symmetry_type::ORIGIN },
    SymmetryName{ "PI_SYM_NOPARM", symmetry_type::PI_SYM_NO_PARAM },
    SymmetryName{ "PI_SYM",        symmetry_type::PI_SYM },
    SymmetryName{ "XAXIS_NOIMAG",  symmetry_type::X_AXIS_NO_IMAG },
    SymmetryName{ "XAXIS_NOREAL",  symmetry_type::X_AXIS_NO_REAL },
    SymmetryName{ "NOPLOT",        symmetry_type::NO_PLOT },
};

inline void push_pending_op(FunctionPtr f, int p)
{
    s_op.push_back(PEND_OP{f, p});
    ++g_operation_index;
    assert(g_operation_index == s_op.size());
}

static bool parse_formula_text(char const *text)
{
    ConstArg *c;
    int ModFlag = 999;
    int Len;
    int Equals = 0;
    std::array<int, 20> Mods{};
    int mdstk = 0;
    double const_pi;
    double const_e;
    double Xctr;
    double Yctr;
    double Xmagfactor;
    double Rotation;
    double Skew;
    LDBL Magnification;
    s_set_random = false;
    s_randomized = false;
    g_uses_jump = false;
    g_jump_index = 0;
    g_jump_control.clear();

    switch (s_math_type)
    {
    case D_MATH:
        StkAdd = dStkAdd;
        StkSub = dStkSub;
        StkNeg = dStkNeg;
        StkMul = dStkMul;
        StkSin = dStkSin;
        StkSinh = dStkSinh;
        StkLT = dStkLT;
        StkLTE = dStkLTE;
        StkMod = dStkMod;
        StkSqr = dStkSqr;
        StkCos = dStkCos;
        StkCosh = dStkCosh;
        StkLog = dStkLog;
        StkExp = dStkExp;
        StkPwr = dStkPwr;
        StkDiv = dStkDiv;
        StkAbs = dStkAbs;
        StkReal = dStkReal;
        StkImag = dStkImag;
        StkConj = dStkConj;
        StkTrig0 = dtrig0;
        StkTrig1 = dtrig1;
        StkTrig2 = dtrig2;
        StkTrig3 = dtrig3;
        StkFlip = dStkFlip;
        StkTan = dStkTan;
        StkTanh = dStkTanh;
        StkCoTan = dStkCoTan;
        StkCoTanh = dStkCoTanh;
        StkCosXX = dStkCosXX;
        StkGT  = dStkGT;
        StkGTE = dStkGTE;
        StkEQ  = dStkEQ;
        StkNE  = dStkNE;
        StkAND = dStkAND;
        StkOR  = dStkOR ;
        StkSRand = dStkSRand;
        StkASin = dStkASin;
        StkASinh = dStkASinh;
        StkACos = dStkACos;
        StkACosh = dStkACosh;
        StkATan = dStkATan;
        StkATanh = dStkATanh;
        StkCAbs = dStkCAbs;
        StkSqrt = dStkSqrt;
        StkZero = dStkZero;
        StkFloor = dStkFloor;
        StkCeil = dStkCeil;
        StkTrunc = dStkTrunc;
        StkRound = dStkRound;
        StkJumpOnTrue  = dStkJumpOnTrue;
        StkJumpOnFalse = dStkJumpOnFalse;
        StkOne = dStkOne;
        break;
    case M_MATH:
        StkAdd = mStkAdd;
        StkSub = mStkSub;
        StkNeg = mStkNeg;
        StkMul = mStkMul;
        StkSin = mStkSin;
        StkSinh = mStkSinh;
        StkLT = mStkLT;
        StkLTE = mStkLTE;
        StkMod = mStkMod;
        StkSqr = mStkSqr;
        StkCos = mStkCos;
        StkCosh = mStkCosh;
        StkLog = mStkLog;
        StkExp = mStkExp;
        StkPwr = mStkPwr;
        StkDiv = mStkDiv;
        StkAbs = mStkAbs;
        StkReal = mStkReal;
        StkImag = mStkImag;
        StkConj = mStkConj;
        StkTrig0 = mtrig0;
        StkTrig1 = mtrig1;
        StkTrig2 = mtrig2;
        StkTrig3 = mtrig3;
        StkFlip = mStkFlip;
        StkTan  = mStkTan;
        StkTanh  = mStkTanh;
        StkCoTan  = mStkCoTan;
        StkCoTanh  = mStkCoTanh;
        StkCosXX = mStkCosXX;
        StkGT  = mStkGT;
        StkGTE = mStkGTE;
        StkEQ  = mStkEQ;
        StkNE  = mStkNE;
        StkAND = mStkAND;
        StkOR  = mStkOR ;
        StkSRand = mStkSRand;
        StkASin = mStkASin;
        StkACos = mStkACos;
        StkACosh = mStkACosh;
        StkATan = mStkATan;
        StkATanh = mStkATanh;
        StkCAbs = mStkCAbs;
        StkSqrt = mStkSqrt;
        StkZero = mStkZero;
        StkFloor = mStkFloor;
        StkCeil = mStkCeil;
        StkTrunc = mStkTrunc;
        StkRound = mStkRound;
        StkJumpOnTrue  = mStkJumpOnTrue;
        StkJumpOnFalse = mStkJumpOnFalse;
        StkOne = mStkOne;
        break;
    case L_MATH:
        s_delta16 = g_bit_shift - 16;
        s_shift_back = 32 - g_bit_shift;
        StkAdd = lStkAdd;
        StkSub = lStkSub;
        StkNeg = lStkNeg;
        StkMul = lStkMul;
        StkSin = lStkSin;
        StkSinh = lStkSinh;
        StkLT = lStkLT;
        StkLTE = lStkLTE;
        StkMod = lStkMod;
        StkSqr = lStkSqr;
        StkCos = lStkCos;
        StkCosh = lStkCosh;
        StkLog = lStkLog;
        StkExp = lStkExp;
        StkPwr = lStkPwr;
        StkDiv = lStkDiv;
        StkAbs = lStkAbs;
        StkReal = lStkReal;
        StkImag = lStkImag;
        StkConj = lStkConj;
        StkTrig0 = ltrig0;
        StkTrig1 = ltrig1;
        StkTrig2 = ltrig2;
        StkTrig3 = ltrig3;
        StkFlip = lStkFlip;
        StkTan  = lStkTan;
        StkTanh  = lStkTanh;
        StkCoTan  = lStkCoTan;
        StkCoTanh  = lStkCoTanh;
        StkCosXX = lStkCosXX;
        StkGT  = lStkGT;
        StkGTE = lStkGTE;
        StkEQ  = lStkEQ;
        StkNE  = lStkNE;
        StkAND = lStkAND;
        StkOR  = lStkOR ;
        StkSRand = lStkSRand;
        StkASin = lStkASin;
        StkACos = lStkACos;
        StkACosh = lStkACosh;
        StkATan = lStkATan;
        StkATanh = lStkATanh;
        StkCAbs = lStkCAbs;
        StkSqrt = lStkSqrt;
        StkZero = lStkZero;
        StkFloor = lStkFloor;
        StkCeil = lStkCeil;
        StkTrunc = lStkTrunc;
        StkRound = lStkRound;
        StkJumpOnTrue  = lStkJumpOnTrue;
        StkJumpOnFalse = lStkJumpOnFalse;
        StkOne = lStkOne;
        break;
    }
    g_max_function = 0;
    for (g_variable_index = 0; g_variable_index < static_cast<unsigned>(s_variables.size()); g_variable_index++)
    {
        v[g_variable_index].s = s_variables[g_variable_index];
        v[g_variable_index].len = (int) std::strlen(s_variables[g_variable_index]);
    }
    cvtcentermag(&Xctr, &Yctr, &Magnification, &Xmagfactor, &Rotation, &Skew);
    const_pi = std::atan(1.0) * 4;
    const_e  = std::exp(1.0);
    v[7].a.d.y = 0.0;
    v[7].a.d.x = v[7].a.d.y;
    v[11].a.d.x = (double)g_logical_screen_x_dots;
    v[11].a.d.y = (double)g_logical_screen_y_dots;
    v[12].a.d.x = (double)g_max_iterations;
    v[12].a.d.y = 0;
    v[13].a.d.x = g_is_mandelbrot ? 1.0 : 0.0;
    v[13].a.d.y = 0;
    v[14].a.d.x = Xctr;
    v[14].a.d.y = Yctr;
    v[15].a.d.x = (double)Magnification;
    v[15].a.d.y = Xmagfactor;
    v[16].a.d.x = Rotation;
    v[16].a.d.y = Skew;

    switch (s_math_type)
    {
    case D_MATH:
        v[1].a.d.x = g_params[0];
        v[1].a.d.y = g_params[1];
        v[2].a.d.x = g_params[2];
        v[2].a.d.y = g_params[3];
        v[5].a.d.x = const_pi;
        v[5].a.d.y = 0.0;
        v[6].a.d.x = const_e;
        v[6].a.d.y = 0.0;
        v[8].a.d.x = g_params[4];
        v[8].a.d.y = g_params[5];
        v[17].a.d.x = g_params[6];
        v[17].a.d.y = g_params[7];
        v[18].a.d.x = g_params[8];
        v[18].a.d.y = g_params[9];
        break;
    case M_MATH:
        v[1].a.m.x = *d2MP(g_params[0]);
        v[1].a.m.y = *d2MP(g_params[1]);
        v[2].a.m.x = *d2MP(g_params[2]);
        v[2].a.m.y = *d2MP(g_params[3]);
        v[5].a.m.x = *d2MP(const_pi);
        v[5].a.m.y = *d2MP(0.0);
        v[6].a.m.x = *d2MP(const_e);
        v[6].a.m.y = *d2MP(0.0);
        v[8].a.m.x = *d2MP(g_params[4]);
        v[8].a.m.y = *d2MP(g_params[5]);
        v[11].a.m  = cmplx2MPC(v[11].a.d);
        v[12].a.m  = cmplx2MPC(v[12].a.d);
        v[13].a.m  = cmplx2MPC(v[13].a.d);
        v[14].a.m  = cmplx2MPC(v[14].a.d);
        v[15].a.m  = cmplx2MPC(v[15].a.d);
        v[16].a.m  = cmplx2MPC(v[16].a.d);
        v[17].a.m.x = *d2MP(g_params[6]);
        v[17].a.m.y = *d2MP(g_params[7]);
        v[18].a.m.x = *d2MP(g_params[8]);
        v[18].a.m.y = *d2MP(g_params[9]);
        break;
    case L_MATH:
        v[1].a.l.x = (long)(g_params[0] * s_fudge);
        v[1].a.l.y = (long)(g_params[1] * s_fudge);
        v[2].a.l.x = (long)(g_params[2] * s_fudge);
        v[2].a.l.y = (long)(g_params[3] * s_fudge);
        v[5].a.l.x = (long)(const_pi * s_fudge);
        v[5].a.l.y = 0L;
        v[6].a.l.x = (long)(const_e * s_fudge);
        v[6].a.l.y = 0L;
        v[8].a.l.x = (long)(g_params[4] * s_fudge);
        v[8].a.l.y = (long)(g_params[5] * s_fudge);
        v[11].a.l.x = g_logical_screen_x_dots;
        v[11].a.l.x <<= g_bit_shift;
        v[11].a.l.y = g_logical_screen_y_dots;
        v[11].a.l.y <<= g_bit_shift;
        v[12].a.l.x = g_max_iterations;
        v[12].a.l.x <<= g_bit_shift;
        v[12].a.l.y = 0L;
        v[13].a.l.x = g_is_mandelbrot ? 1 : 0;
        v[13].a.l.x <<= g_bit_shift;
        v[13].a.l.y = 0L;
        v[14].a.l.x = (long)(v[14].a.d.x * s_fudge);
        v[14].a.l.y = (long)(v[14].a.d.y * s_fudge);
        v[15].a.l.x = (long)(v[15].a.d.x * s_fudge);
        v[15].a.l.y = (long)(v[15].a.d.y * s_fudge);
        v[16].a.l.x = (long)(v[16].a.d.x * s_fudge);
        v[16].a.l.y = (long)(v[16].a.d.y * s_fudge);
        v[17].a.l.x = (long)(g_params[6] * s_fudge);
        v[17].a.l.y = (long)(g_params[7] * s_fudge);
        v[18].a.l.x = (long)(g_params[8] * s_fudge);
        v[18].a.l.y = (long)(g_params[9] * s_fudge);
        break;
    }

    g_operation_index = 0;
    s_op.clear();
    g_store_index = 0;
    g_load_index = 0;
    g_op_ptr = 0;
    s_paren = 0;
    g_last_init_op = 0;
    s_expecting_arg = true;
    for (s_n = 0; text[s_n]; s_n++)
    {
        if (!text[s_n])
        {
            break;
        }
        s_init_n = s_n;
        switch (text[s_n])
        {
        case ' ':
        case '\t':
        case '\r':
        case '\n':
            break;
        case '(':
            s_paren++;
            break;
        case ')':
            s_paren--;
            break;
        case '|':
            if (text[s_n+1] == '|')
            {
                s_expecting_arg = true;
                s_n++;
                push_pending_op(StkOR, 7 - (s_paren + Equals) * 15);
            }
            else if (ModFlag == s_paren-1)
            {
                s_paren--;
                ModFlag = Mods[--mdstk];
            }
            else
            {
                assert(mdstk < Mods.size());
                Mods[mdstk++] = ModFlag;
                push_pending_op(StkMod, 2 - (s_paren + Equals) * 15);
                ModFlag = s_paren++;
            }
            break;
        case ',':
        case ';':
            if (!s_expecting_arg)
            {
                s_expecting_arg = true;
                push_pending_op(nullptr, 15);
                push_pending_op(StkClr, -30000);
                s_paren = 0;
                Equals = s_paren;
            }
            break;
        case ':':
            s_expecting_arg = true;
            push_pending_op(nullptr, 15);
            push_pending_op(EndInit, -30000);
            s_paren = 0;
            Equals = s_paren;
            g_last_init_op = 10000;
            break;
        case '+':
            s_expecting_arg = true;
            push_pending_op(StkAdd, 4 - (s_paren + Equals)*15);
            break;
        case '-':
            if (s_expecting_arg)
            {
                push_pending_op(StkNeg, 2 - (s_paren + Equals)*15);
            }
            else
            {
                push_pending_op(StkSub, 4 - (s_paren + Equals)*15);
                s_expecting_arg = true;
            }
            break;
        case '&':
            s_expecting_arg = true;
            s_n++;
            push_pending_op(StkAND, 7 - (s_paren + Equals)*15);
            break;
        case '!':
            s_expecting_arg = true;
            s_n++;
            push_pending_op(StkNE, 6 - (s_paren + Equals)*15);
            break;
        case '<':
            s_expecting_arg = true;
            {
                FunctionPtr fn;
                if (text[s_n + 1] == '=')
                {
                    s_n++;
                    fn = StkLTE;
                }
                else
                {
                    fn = StkLT;
                }
                push_pending_op(fn, 6 - (s_paren + Equals) * 15);
            }
            break;
        case '>':
            s_expecting_arg = true;
            {
                FunctionPtr fn;
                if (text[s_n + 1] == '=')
                {
                    s_n++;
                    fn = StkGTE;
                }
                else
                {
                    fn = StkGT;
                }
                push_pending_op(fn, 6 - (s_paren + Equals) * 15);
            }
            break;
        case '*':
            s_expecting_arg = true;
            push_pending_op(StkMul, 3 - (s_paren + Equals)*15);
            break;
        case '/':
            s_expecting_arg = true;
            push_pending_op(StkDiv, 3 - (s_paren + Equals)*15);
            break;
        case '^':
            s_expecting_arg = true;
            push_pending_op(StkPwr, 2 - (s_paren + Equals)*15);
            break;
        case '=':
            s_expecting_arg = true;
            if (text[s_n+1] == '=')
            {
                s_n++;
                push_pending_op(StkEQ, 6 - (s_paren + Equals)*15);
            }
            else
            {
                s_op[g_operation_index-1].f = StkSto;
                s_op[g_operation_index-1].p = 5 - (s_paren + Equals)*15;
                Store[g_store_index++] = Load[--g_load_index];
                Equals++;
            }
            break;
        default:
            while (std::isalnum(text[s_n+1]) || text[s_n+1] == '.' || text[s_n+1] == '_')
            {
                s_n++;
            }
            Len = (s_n+1)-s_init_n;
            s_expecting_arg = false;
            if (const jump_control_type type = is_jump(&text[s_init_n], Len); type != jump_control_type::NONE)
            {
                g_uses_jump = true;
                switch (type)
                {
                case jump_control_type::IF:
                    s_expecting_arg = true;
                    push_jump(jump_control_type::IF);
                    push_pending_op(StkJumpOnFalse, 1);
                    break;
                case jump_control_type::ELSE_IF:
                    s_expecting_arg = true;
                    push_jump(jump_control_type::ELSE_IF);
                    push_jump(jump_control_type::ELSE_IF);
                    push_pending_op(StkJump, 1);
                    push_pending_op(nullptr, 15);
                    push_pending_op(StkClr, -30000);
                    push_pending_op(StkJumpOnFalse, 1);
                    break;
                case jump_control_type::ELSE:
                    push_jump(jump_control_type::ELSE);
                    push_pending_op(StkJump, 1);
                    break;
                case jump_control_type::END_IF:
                    push_jump(jump_control_type::END_IF);
                    push_pending_op(StkJumpLabel, 1);
                    break;
                default:
                    break;
                }
            }
            else
            {
                if (const FunctionPtr fn = is_func(&text[s_init_n], Len); fn != NotAFnct)
                {
                    push_pending_op(fn,  1 - (s_paren + Equals)*15);
                    s_expecting_arg = true;
                }
                else
                {
                    c = is_const(&text[s_init_n], Len);
                    Load[g_load_index++] = &(c->a);
                    push_pending_op(StkLod, 1 - (s_paren + Equals)*15);
                    s_n = s_init_n + c->len - 1;
                }
            }
            break;
        }
    }
    push_pending_op(nullptr, 16);
    s_next_op = 0;
    g_last_op = g_operation_index;
    while (s_next_op < g_operation_index)
    {
        if (s_op[s_next_op].f)
        {
            RecSortPrec();
        }
        else
        {
            s_next_op++;
            g_last_op--;
        }
    }
    return false;
}

int Formula()
{
    if (g_formula_name.empty() || g_overflow)
    {
        return 1;
    }

    g_load_index = InitLodPtr;
    g_store_index = InitStoPtr;
    g_op_ptr = InitOpPtr;
    g_jump_index = s_init_jump_index;
    // Set the random number
    if (s_set_random || s_randomized)
    {
        switch (s_math_type)
        {
        case D_MATH:
            dRandom();
            break;
        case L_MATH:
            lRandom();
            break;
        case M_MATH:
            mRandom();
        }
    }

    Arg1 = &g_stack[0];
    Arg2 = &g_stack[0];
    --Arg2;
    while (g_op_ptr < (int)g_last_op)
    {
        f[g_op_ptr]();
        g_op_ptr++;
    }

    switch (s_math_type)
    {
    case D_MATH:
        g_new_z = v[3].a.d;
        g_old_z = g_new_z;
        return Arg1->d.x == 0.0;
    case M_MATH:
        g_new_z = MPC2cmplx(v[3].a.m);
        g_old_z = g_new_z;
        return Arg1->m.x.Exp == 0 && Arg1->m.x.Mant == 0;
    case L_MATH:
        g_l_new_z = v[3].a.l;
        g_l_old_z = g_l_new_z;
        if (g_overflow)
        {
            return 1;
        }
        return Arg1->l.x == 0L;
    }
    return 1;
}

int form_per_pixel()
{
    if (g_formula_name.empty())
    {
        return 1;
    }
    g_overflow = false;
    g_jump_index = 0;
    g_op_ptr = 0;
    g_store_index = 0;
    g_load_index = 0;
    Arg1 = &g_stack[0];
    Arg2 = &g_stack[0];
    Arg2--;


    v[10].a.d.x = (double)g_col;
    v[10].a.d.y = (double)g_row;

    switch (s_math_type)
    {
    case D_MATH:
        if ((g_row+g_col)&1)
        {
            v[9].a.d.x = 1.0;
        }
        else
        {
            v[9].a.d.x = 0.0;
        }
        v[9].a.d.y = 0.0;
        break;
    
    case M_MATH:
        if ((g_row+g_col)&1)
        {
            v[9].a.m = g_mpc_one;
        }
        else
        {
            v[9].a.m.x.Exp = 0;
            v[9].a.m.x.Mant = 0;
            v[9].a.m.y.Exp = 0;
            v[9].a.m.y.Mant = 0;
        }
        v[10].a.m = cmplx2MPC(v[10].a.d);
        break;

    case L_MATH:
        v[9].a.l.x = (long)(((g_row+g_col)&1) * s_fudge);
        v[9].a.l.y = 0L;
        v[10].a.l.x = g_col;
        v[10].a.l.x <<= g_bit_shift;
        v[10].a.l.y = g_row;
        v[10].a.l.y <<= g_bit_shift;
        break;
    }

    {
        if (g_invert != 0)
        {
            invertz2(&g_old_z);
            switch (s_math_type)
            {
            case D_MATH:
                v[0].a.d.x = g_old_z.x;
                v[0].a.d.y = g_old_z.y;
                break;
            case M_MATH:
                v[0].a.m.x = *d2MP(g_old_z.x);
                v[0].a.m.y = *d2MP(g_old_z.y);
                break;
            case L_MATH:
                // watch out for overflow
                if (sqr(g_old_z.x)+sqr(g_old_z.y) >= 127)
                {
                    g_old_z.x = 8;  // value to bail out in one iteration
                    g_old_z.y = 8;
                }
                // convert to fudged longs
                v[0].a.l.x = (long)(g_old_z.x*s_fudge);
                v[0].a.l.y = (long)(g_old_z.y*s_fudge);
                break;
            }
        }
        else
        {
            switch (s_math_type)
            {
            case D_MATH:
                v[0].a.d.x = g_dx_pixel();
                v[0].a.d.y = g_dy_pixel();
                break;
            case M_MATH:
                v[0].a.m.x = *d2MP(g_dx_pixel());
                v[0].a.m.y = *d2MP(g_dy_pixel());
                break;
            case L_MATH:
                v[0].a.l.x = g_l_x_pixel();
                v[0].a.l.y = g_l_y_pixel();
                break;
            }
        }
    }

    if (g_last_init_op)
    {
        g_last_init_op = g_last_op;
    }
    while (g_op_ptr < g_last_init_op)
    {
        f[g_op_ptr]();
        g_op_ptr++;
    }
    InitLodPtr = g_load_index;
    InitStoPtr = g_store_index;
    InitOpPtr = g_op_ptr;
    // Set old variable for orbits
    switch (s_math_type)
    {
    case D_MATH:
        g_old_z = v[3].a.d;
        break;
    case M_MATH:
        g_old_z = MPC2cmplx(v[3].a.m);
        break;
    case L_MATH:
        g_l_old_z = v[3].a.l;
        break;
    }

    return g_overflow ? 0 : 1;
}

static int fill_if_group(int endif_index, JUMP_PTRS_ST* jump_data)
{
    int i   = endif_index;
    int ljp = endif_index; // ljp means "last jump processed"
    while (i > 0)
    {
        i--;
        switch (g_jump_control[i].type)
        {
        case jump_control_type::IF:    //if (); this concludes processing of this group
            g_jump_control[i].ptrs = jump_data[ljp];
            g_jump_control[i].DestJumpIndex = ljp + 1;
            return i;
        case jump_control_type::ELSE_IF:    //elseif* ( 2 jumps, the else and the if
            // first, the "if" part
            g_jump_control[i].ptrs = jump_data[ljp];
            g_jump_control[i].DestJumpIndex = ljp + 1;

            // then, the else part
            i--; //fall through to "else" is intentional
        case jump_control_type::ELSE:
            g_jump_control[i].ptrs = jump_data[endif_index];
            g_jump_control[i].DestJumpIndex = endif_index + 1;
            ljp = i;
            break;
        case jump_control_type::END_IF:    //endif
            i = fill_if_group(i, jump_data);
            break;
        default:
            break;
        }
    }
    return -1; // should never get here
}

static bool fill_jump_struct()
{
    // Completes all entries in jump structure. Returns 1 on error)
    // On entry, jump_index is the number of jump functions in the formula
    int i = 0;
    int loadcount = 0;
    int storecount = 0;
    bool checkforelse = false;
    FunctionPtr JumpFunc = nullptr;
    bool find_new_func = true;

    std::vector<JUMP_PTRS_ST> jump_data;

    for (g_op_ptr = 0; g_op_ptr < (int) g_last_op; g_op_ptr++)
    {
        if (find_new_func)
        {
            if (i < static_cast<int>(g_jump_control.size()))
            {
                switch (g_jump_control[i].type)
                {
                case jump_control_type::IF:
                    JumpFunc = StkJumpOnFalse;
                    break;
                case jump_control_type::ELSE_IF:
                    checkforelse = !checkforelse;
                    if (checkforelse)
                    {
                        JumpFunc = StkJump;
                    }
                    else
                    {
                        JumpFunc = StkJumpOnFalse;
                    }
                    break;
                case jump_control_type::ELSE:
                    JumpFunc = StkJump;
                    break;
                case jump_control_type::END_IF:
                    JumpFunc = StkJumpLabel;
                    break;
                default:
                    break;
                }
            }
            find_new_func = false;
        }
        if (*(f[g_op_ptr]) == StkLod)
        {
            loadcount++;
        }
        else if (*(f[g_op_ptr]) == StkSto)
        {
            storecount++;
        }
        else if (*(f[g_op_ptr]) == JumpFunc)
        {
            JUMP_PTRS_ST value{};
            value.JumpOpPtr = g_op_ptr;
            value.JumpLodPtr = loadcount;
            value.JumpStoPtr = storecount;
            jump_data.push_back(value);
            i++;
            find_new_func = true;
        }
    }

    // Following for safety only; all should always be false
    if (i != g_jump_index
        || g_jump_control[i - 1].type != jump_control_type::END_IF
        || g_jump_control[0].type != jump_control_type::IF)
    {
        return true;
    }

    while (i > 0)
    {
        i--;
        i = fill_if_group(i, jump_data.data());
    }
    return i < 0;
}

static std::string s_formula;

static int frmgetchar(std::FILE * openfile)
{
    int c;
    bool done = false;
    bool linewrap = false;
    while (!done)
    {
        c = getc(openfile);
        switch (c)
        {
        case '\r':
        case ' ' :
        case '\t' :
            break;
        case '\\':
            linewrap = true;
            break;
        case ';' :
            while ((c = getc(openfile)) != '\n' && c != EOF)
            {
            }
            if (c == EOF)
            {
                done = true;
            }
        case '\n' :
            if (!linewrap)
            {
                done = true;
            }
            linewrap = false;
            break;
        default:
            done = true;
            break;
        }
    }
    return std::tolower(c);
}

// This function also gets flow control info

static void getfuncinfo(token_st * tok)
{
    for (int i = 0; i < static_cast<int>(s_func_list.size()); i++)
    {
        if (!std::strcmp(s_func_list[i].s, tok->str))
        {
            tok->id = static_cast<token_id>(i);
            if (tok->id >= token_id::FUNC_FN1 && tok->id <= token_id::FUNC_FN4)
            {
                tok->type = token_type::PARAM_FUNCTION;
            }
            else
            {
                tok->type = token_type::FUNCTION;
            }
            return;
        }
    }

    for (int i = 0; i < static_cast<int>(s_jump_list.size()); i++) // pick up flow control
    {
        if (std::strcmp(s_jump_list[i], tok->str) == 0)
        {
            tok->type = token_type::FLOW_CONTROL;
            tok->id   = static_cast<token_id>(i + 1);
            return;
        }
    }
    tok->type = token_type::NOT_A_TOKEN;
    tok->id   = token_id::UNDEFINED_FUNCTION;
}

static void getvarinfo(token_st * tok)
{
    for (int i = 0; i < static_cast<int>(s_variables.size()); i++)
    {
        if (!std::strcmp(s_variables[i], tok->str))
        {
            tok->id = static_cast<token_id>(i);
            switch (tok->id)
            {
            case token_id::VAR_P1:
            case token_id::VAR_P2:
            case token_id::VAR_P3:
            case token_id::VAR_IS_MANDEL:
            case token_id::VAR_P4:
            case token_id::VAR_P5:
                tok->type = token_type::PARAM_VARIABLE;
                break;
            default:
                tok->type = token_type::PREDEFINED_VARIABLE;
                break;
            }
            return;
        }
    }
    tok->type = token_type::USER_NAMED_VARIABLE;
    tok->id   = token_id::NONE;
}

// fills in token structure where numeric constant is indicated
/* Note - this function will be called twice to fill in the components
        of a complex constant. See is_complex_constant() below. */

// returns 1 on success, 0 on NOT_A_TOKEN
static bool frmgetconstant(std::FILE * openfile, token_st * tok)
{
    int c;
    int i = 1;
    bool getting_base = true;
    long filepos = ftell(openfile);
    bool got_decimal_already = false;
    bool done = false;
    tok->constant.x = 0.0;          //initialize values to 0
    tok->constant.y = 0.0;
    if (tok->str[0] == '.')
    {
        got_decimal_already = true;
    }
    while (!done)
    {
        c = frmgetchar(openfile);
        switch (c)
        {
        case EOF:
            tok->str[i] = (char) 0;
            tok->type = token_type::NOT_A_TOKEN;
            tok->id   = token_id::END_OF_FILE;
            return false;
        CASE_NUM:
            tok->str[i++] = (char) c;
            filepos = ftell(openfile);
            break;
        case '.':
            if (got_decimal_already || !getting_base)
            {
                tok->str[i++] = (char) c;
                tok->str[i++] = (char) 0;
                tok->type = token_type::NOT_A_TOKEN;
                tok->id = token_id::ILL_FORMED_CONSTANT;
                return false;
            }
            tok->str[i++] = (char) c;
            got_decimal_already = true;
            filepos = ftell(openfile);
            break;
        default :
            if (c == 'e' && getting_base && (std::isdigit(tok->str[i-1]) || (tok->str[i-1] == '.' && i > 1)))
            {
                tok->str[i++] = (char) c;
                getting_base = false;
                got_decimal_already = false;
                filepos = ftell(openfile);
                c = frmgetchar(openfile);
                if (c == '-' || c == '+')
                {
                    tok->str[i++] = (char) c;
                    filepos = ftell(openfile);
                }
                else
                {
                    fseek(openfile, filepos, SEEK_SET);
                }
            }
            else if (std::isalpha(c) || c == '_')
            {
                tok->str[i++] = (char) c;
                tok->str[i++] = (char) 0;
                tok->type = token_type::NOT_A_TOKEN;
                tok->id = token_id::ILL_FORMED_CONSTANT;
                return false;
            }
            else if (tok->str[i-1] == 'e' || (tok->str[i-1] == '.' && i == 1))
            {
                tok->str[i++] = (char) c;
                tok->str[i++] = (char) 0;
                tok->type = token_type::NOT_A_TOKEN;
                tok->id = token_id::ILL_FORMED_CONSTANT;
                return false;
            }
            else
            {
                fseek(openfile, filepos, SEEK_SET);
                tok->str[i++] = (char) 0;
                done = true;
            }
            break;
        }
        if (i == 33 && tok->str[32])
        {
            tok->str[33] = (char) 0;
            tok->type = token_type::NOT_A_TOKEN;
            tok->id = token_id::TOKEN_TOO_LONG;
            return false;
        }
    }    // end of while loop. Now fill in the value
    tok->constant.x = std::atof(tok->str);
    tok->type = token_type::REAL_CONSTANT;
    tok->id   = token_id::NONE;
    return true;
}

static void is_complex_constant(std::FILE * openfile, token_st * tok)
{
    assert(tok->str[0] == '(');
    token_st temp_tok;
    long filepos;
    int c;
    int sign_value = 1;
    bool done = false;
    bool getting_real = true;
    std::FILE * debug_token = nullptr;
    tok->str[1] = (char) 0;  // so we can concatenate later

    filepos = ftell(openfile);
    if (g_debug_flag == debug_flags::write_formula_debug_information)
    {
        debug_token = std::fopen("frmconst.txt", "at");
    }

    while (!done)
    {
        c = frmgetchar(openfile);
        switch (c)
        {
CASE_NUM :
        case '.':
            if (debug_token != nullptr)
            {
                std::fprintf(debug_token,  "Set temp_tok.token_str[0] to %c\n", c);
            }
            temp_tok.str[0] = (char) c;
            break;
        case '-' :
            if (debug_token != nullptr)
            {
                std::fprintf(debug_token,  "First char is a minus\n");
            }
            sign_value = -1;
            c = frmgetchar(openfile);
            if (c == '.' || std::isdigit(c))
            {
                if (debug_token != nullptr)
                {
                    std::fprintf(debug_token,  "Set temp_tok.token_str[0] to %c\n", c);
                }
                temp_tok.str[0] = (char) c;
            }
            else
            {
                if (debug_token != nullptr)
                {
                    std::fprintf(debug_token,  "First char not a . or NUM\n");
                }
                done = true;
            }
            break;
        default:
            if (debug_token != nullptr)
            {
                std::fprintf(debug_token,  "First char not a . or NUM\n");
            }
            done = true;
            break;
        }
        if (debug_token != nullptr)
        {
            std::fprintf(debug_token,  "Calling frmgetconstant unless done is true; done is %s\n", done ? "true" : "false");
        }
        if (!done && frmgetconstant(openfile, &temp_tok))
        {
            c = frmgetchar(openfile);
            if (debug_token != nullptr)
            {
                std::fprintf(debug_token, "frmgetconstant returned 1; next token is %c\n", c);
            }
            if (getting_real && c == ',')   // we have the real part now
            {
                if (sign_value == -1)
                {
                    std::strcat(tok->str, "-");
                }
                std::strcat(tok->str, temp_tok.str);
                std::strcat(tok->str, ",");
                tok->constant.x = temp_tok.constant.x * sign_value;
                getting_real = false;
                sign_value = 1;
            }
            else if (!getting_real && c == ')') // we have the complex part
            {
                if (sign_value == -1)
                {
                    std::strcat(tok->str, "-");
                }
                std::strcat(tok->str, temp_tok.str);
                std::strcat(tok->str, ")");
                tok->constant.y = temp_tok.constant.x * sign_value;
                tok->type = tok->constant.y ? token_type::COMPLEX_CONSTANT : token_type::REAL_CONSTANT;
                tok->id   = token_id::NONE;
                if (debug_token != nullptr)
                {
                    std::fprintf(debug_token, "Exiting with type set to %d\n",
                        +(tok->constant.y ? token_type::COMPLEX_CONSTANT : token_type::REAL_CONSTANT));
                    std::fclose(debug_token);
                }
                return;
            }
            else
            {
                done = true;
            }
        }
        else
        {
            done = true;
        }
    }
    fseek(openfile, filepos, SEEK_SET);
    tok->str[1] = (char) 0;
    tok->constant.x = 0.0;
    tok->constant.y = tok->constant.x;
    tok->type = token_type::PARENS;
    tok->id = token_id::OPEN_PARENS;
    if (debug_token != nullptr)
    {
        std::fprintf(debug_token,  "Exiting with ID set to OPEN_PARENS\n");
        std::fclose(debug_token);
    }
}

static bool frmgetalpha(std::FILE * openfile, token_st * tok)
{
    int c;
    int i = 1;
    bool var_name_too_long = false;
    long filepos;
    long last_filepos = ftell(openfile);
    while ((c = frmgetchar(openfile)) != EOF)
    {
        filepos = ftell(openfile);
        switch (c)
        {
CASE_ALPHA:
CASE_NUM:
        case '_':
            if (i < 79)
            {
                tok->str[i++] = (char) c;
            }
            else
            {
                tok->str[i] = (char) 0;
            }
            if (i == 33)
            {
                var_name_too_long = true;
            }
            last_filepos = filepos;
            break;
        default:
            if (c == '.')       // illegal character in variable or func name
            {
                tok->type = token_type::NOT_A_TOKEN;
                tok->id   = token_id::ILLEGAL_VARIABLE_NAME;
                tok->str[i++] = '.';
                tok->str[i] = (char) 0;
                return false;
            }
            else if (var_name_too_long)
            {
                tok->type = token_type::NOT_A_TOKEN;
                tok->id   = token_id::TOKEN_TOO_LONG;
                tok->str[i] = (char) 0;
                fseek(openfile, last_filepos, SEEK_SET);
                return false;
            }
            tok->str[i] = (char) 0;
            std::fseek(openfile, last_filepos, SEEK_SET);
            getfuncinfo(tok);
            if (c == '(') //getfuncinfo() correctly filled structure
            {
                if (tok->type == token_type::NOT_A_TOKEN)
                {
                    return false;
                }
                if (tok->type == token_type::FLOW_CONTROL &&
                    (tok->id == token_id::ILLEGAL_VARIABLE_NAME || tok->id == token_id::TOKEN_TOO_LONG))
                {
                    tok->type = token_type::NOT_A_TOKEN;
                    tok->id   = token_id::JUMP_WITH_ILLEGAL_CHAR;
                    return false;
                }
                return true;
            }
            //can't use function names as variables
            if (tok->type == token_type::FUNCTION || tok->type == token_type::PARAM_FUNCTION)
            {
                tok->type = token_type::NOT_A_TOKEN;
                tok->id   = token_id::FUNC_USED_AS_VAR;
                return false;
            }
            if (tok->type == token_type::FLOW_CONTROL &&
                (tok->id == token_id::END_OF_FILE || tok->id == token_id::ILLEGAL_CHARACTER))
            {
                tok->type = token_type::NOT_A_TOKEN;
                tok->id   = token_id::JUMP_MISSING_BOOLEAN;
                return false;
            }
            if (tok->type == token_type::FLOW_CONTROL &&
                (tok->id == token_id::ILLEGAL_VARIABLE_NAME || tok->id == token_id::TOKEN_TOO_LONG))
            {
                if (c == ',' || c == '\n' || c == ':')
                {
                    return true;
                }
                tok->type = token_type::NOT_A_TOKEN;
                tok->id   = token_id::JUMP_WITH_ILLEGAL_CHAR;
                return false;
            }
            getvarinfo(tok);
            return true;
        }
    }
    tok->str[0] = (char) 0;
    tok->type = token_type::NOT_A_TOKEN;
    tok->id   = token_id::END_OF_FILE;
    return false;
}

static void frm_get_eos(std::FILE * openfile, token_st * this_token)
{
    long last_filepos = ftell(openfile);
    int c;

    for (c = frmgetchar(openfile); (c == '\n' || c == ',' || c == ':'); c = frmgetchar(openfile))
    {
        if (c == ':')
        {
            this_token->str[0] = ':';
        }
        last_filepos = ftell(openfile);
    }
    if (c == '}')
    {
        this_token->str[0] = '}';
        this_token->type = token_type::END_OF_FORMULA;
        this_token->id   = token_id::NONE;
    }
    else
    {
        std::fseek(openfile, last_filepos, SEEK_SET);
        if (this_token->str[0] == '\n')
        {
            this_token->str[0] = ',';
        }
    }
}

/*frmgettoken fills token structure; returns 1 on success and 0 on
  NOT_A_TOKEN and END_OF_FORMULA
*/
static bool frmgettoken(std::FILE * openfile, token_st * this_token)
{
    int i = 1;
    long filepos;

    int c = frmgetchar(openfile);
    switch (c)
    {
CASE_NUM:
    case '.':
        this_token->str[0] = (char) c;
        return frmgetconstant(openfile, this_token);
CASE_ALPHA:
    case '_':
        this_token->str[0] = (char) c;
        return frmgetalpha(openfile, this_token);
CASE_TERMINATOR:
        this_token->type = token_type::OPERATOR; // this may be changed below
        this_token->str[0] = (char) c;
        filepos = ftell(openfile);
        if (c == '<' || c == '>' || c == '=')
        {
            c = frmgetchar(openfile);
            if (c == '=')
            {
                this_token->str[i++] = (char) c;
            }
            else
            {
                std::fseek(openfile, filepos, SEEK_SET);
            }
        }
        else if (c == '!')
        {
            c = frmgetchar(openfile);
            if (c == '=')
            {
                this_token->str[i++] = (char) c;
            }
            else
            {
                std::fseek(openfile, filepos, SEEK_SET);
                this_token->str[1] = (char) 0;
                this_token->type = token_type::NOT_A_TOKEN;
                this_token->id = token_id::ILLEGAL_OPERATOR;
                return false;
            }
        }
        else if (c == '|')
        {
            c = frmgetchar(openfile);
            if (c == '|')
            {
                this_token->str[i++] = (char) c;
            }
            else
            {
                std::fseek(openfile, filepos, SEEK_SET);
            }
        }
        else if (c == '&')
        {
            c = frmgetchar(openfile);
            if (c == '&')
            {
                this_token->str[i++] = (char) c;
            }
            else
            {
                std::fseek(openfile, filepos, SEEK_SET);
                this_token->str[1] = (char) 0;
                this_token->type = token_type::NOT_A_TOKEN;
                this_token->id = token_id::ILLEGAL_OPERATOR;
                return false;
            }
        }
        else if (this_token->str[0] == '}')
        {
            this_token->type = token_type::END_OF_FORMULA;
            this_token->id   = token_id::NONE;
        }
        else if (this_token->str[0] == '\n'
            || this_token->str[0] == ','
            || this_token->str[0] == ':')
        {
            frm_get_eos(openfile, this_token);
        }
        else if (this_token->str[0] == ')')
        {
            this_token->type = token_type::PARENS;
            this_token->id = token_id::CLOSE_PARENS;
        }
        else if (this_token->str[0] == '(')
        {
            /* the following function will set token_type to PARENS and
               token_id to OPEN_PARENS if this is not the start of a
               complex constant */
            is_complex_constant(openfile, this_token);
            return true;
        }
        this_token->str[i] = (char) 0;
        if (this_token->type == token_type::OPERATOR)
        {
            for (int j = 0; j < static_cast<int>(s_op_list.size()); j++)
            {
                if (std::strcmp(s_op_list[j], this_token->str) == 0)
                {
                    this_token->id = static_cast<token_id>(j);
                    break;
                }
            }
        }
        return this_token->str[0] != '}';
    case EOF:
        this_token->str[0] = (char) 0;
        this_token->type = token_type::NOT_A_TOKEN;
        this_token->id = token_id::END_OF_FILE;
        return false;
    default:
        this_token->str[0] = (char) c;
        this_token->str[1] = (char) 0;
        this_token->type = token_type::NOT_A_TOKEN;
        this_token->id = token_id::ILLEGAL_CHARACTER;
        return false;
    }
}

int frm_get_param_stuff(char const *Name)
{
    std::FILE *debug_token = nullptr;
    int c;
    token_st current_token;
    std::FILE * entry_file = nullptr;
    g_frm_uses_p1 = false;
    g_frm_uses_p2 = false;
    g_frm_uses_p3 = false;
    g_frm_uses_ismand = false;
    g_max_function = 0;
    g_frm_uses_p4 = false;
    g_frm_uses_p5 = false;

    if (g_formula_name.empty())
    {
        return 0;  //  and don't reset the pointers
    }
    if (find_file_item(g_formula_filename, Name, &entry_file, gfe_type::FORMULA))
    {
        stopmsg(ParseErrs(PE_COULD_NOT_OPEN_FILE_WHERE_FORMULA_LOCATED));
        return 0;
    }
    while ((c = frmgetchar(entry_file)) != '{' && c != EOF)
    {
    }
    if (c != '{')
    {
        stopmsg(ParseErrs(PE_UNEXPECTED_EOF));
        std::fclose(entry_file);
        return 0;
    }

    if (g_debug_flag == debug_flags::write_formula_debug_information)
    {
        debug_token = std::fopen("frmtokens.txt", "at");
        if (debug_token != nullptr)
        {
            std::fprintf(debug_token, "%s\n", Name);
        }
    }
    while (frmgettoken(entry_file, &current_token))
    {
        if (debug_token != nullptr)
        {
            std::fprintf(debug_token, "%s\n", current_token.str);
            std::fprintf(debug_token, "token_type is %d\n", +current_token.type);
            std::fprintf(debug_token, "token_id is %d\n", +current_token.id);
            if (current_token.type == token_type::REAL_CONSTANT || current_token.type == token_type::COMPLEX_CONSTANT)
            {
                std::fprintf(debug_token, "Real value is %f\n", current_token.constant.x);
                std::fprintf(debug_token, "Imag value is %f\n", current_token.constant.y);
            }
            std::fprintf(debug_token, "\n");
        }
        switch (current_token.type)
        {
        case token_type::PARAM_VARIABLE:
            if (current_token.id == token_id::VAR_P1)
            {
                g_frm_uses_p1 = true;
            }
            else if (current_token.id == token_id::VAR_P2)
            {
                g_frm_uses_p2 = true;
            }
            else if (current_token.id == token_id::VAR_P3)
            {
                g_frm_uses_p3 = true;
            }
            else if (current_token.id == token_id::VAR_IS_MANDEL)
            {
                g_frm_uses_ismand = true;
            }
            else if (current_token.id == token_id::VAR_P4)
            {
                g_frm_uses_p4 = true;
            }
            else if (current_token.id == token_id::VAR_P5)
            {
                g_frm_uses_p5 = true;
            }
            break;
        case token_type::PARAM_FUNCTION:
            if ((+current_token.id - 10) > g_max_function)
            {
                g_max_function = (char)(+current_token.id - 10);
            }
            break;
        default:
            break;
        }
    }
    std::fclose(entry_file);
    if (debug_token)
    {
        std::fclose(debug_token);
    }
    if (current_token.type != token_type::END_OF_FORMULA)
    {
        g_frm_uses_p1 = false;
        g_frm_uses_p2 = false;
        g_frm_uses_p3 = false;
        g_frm_uses_ismand = false;
        g_max_function = 0;
        g_frm_uses_p4 = false;
        g_frm_uses_p5 = false;
        return 0;
    }
    return 1;
}

/* frm_check_name_and_sym():
     error checking to the open brace on the first line; return true
     on success, and false if errors are found which should cause the
     formula not to be executed
*/
static bool frm_check_name_and_sym(std::FILE * open_file, bool report_bad_sym)
{
    long filepos = ftell(open_file);
    int c;
    int i;
    bool at_end_of_name = false;

    // first, test name
    i = 0;
    bool done = false;
    while (!done)
    {
        c = getc(open_file);
        switch (c)
        {
        case EOF:
            stopmsg(ParseErrs(PE_UNEXPECTED_EOF));
            return false;
        case '\r':
        case '\n':
            stopmsg(ParseErrs(PE_NO_LEFT_BRACKET_FIRST_LINE));
            return false;
        case ' ':
        case '\t':
            at_end_of_name = true;
            break;
        case '(':
        case '{':
            done = true;
            break;
        default :
            if (!at_end_of_name)
            {
                i++;
            }
            break;
        }
    }

    if (i > ITEM_NAME_LEN)
    {
        int j;
        int k = (int) std::strlen(ParseErrs(PE_FORMULA_NAME_TOO_LARGE));
        char msgbuf[100];
        std::strcpy(msgbuf, ParseErrs(PE_FORMULA_NAME_TOO_LARGE));
        std::strcat(msgbuf, ":\n   ");
        std::fseek(open_file, filepos, SEEK_SET);
        for (j = 0; j < i && j < 25; j++)
        {
            msgbuf[j+k+2] = (char) getc(open_file);
        }
        msgbuf[j+k+2] = (char) 0;
        stopmsg(stopmsg_flags::FIXED_FONT, msgbuf);
        return false;
    }
    // get symmetry
    g_symmetry = symmetry_type::NONE;
    if (c == '(')
    {
        char sym_buf[20];
        done = false;
        i = 0;
        while (!done)
        {
            c = getc(open_file);
            switch (c)
            {
            case EOF:
                stopmsg(ParseErrs(PE_UNEXPECTED_EOF));
                return false;
            case '\r':
            case '\n':
                stopmsg(stopmsg_flags::FIXED_FONT, ParseErrs(PE_NO_LEFT_BRACKET_FIRST_LINE));
                return false;
            case '{':
                stopmsg(stopmsg_flags::FIXED_FONT, ParseErrs(PE_NO_MATCH_RIGHT_PAREN));
                return false;
            case ' ':
            case '\t':
                break;
            case ')':
                done = true;
                break;
            default :
                if (i < 19)
                {
                    sym_buf[i++] = (char) std::toupper(c);
                }
                break;
            }
        }
        sym_buf[i] = (char) 0;
        constexpr int num_names{s_symmetry_names.size()};
        for (i = 0; i < num_names; i++)
        {
            if (stricmp(s_symmetry_names[i].s, sym_buf) == 0)
            {
                g_symmetry = s_symmetry_names[i].n;
                break;
            }
        }
        if (i == num_names && report_bad_sym)
        {
            std::string msgbuf{ParseErrs(PE_INVALID_SYM_USING_NOSYM)};
            msgbuf += ":\n   ";
            msgbuf += sym_buf;
            stopmsg(stopmsg_flags::FIXED_FONT, msgbuf);
        }
    }
    if (c != '{')
    {
        done = false;
        while (!done)
        {
            c = getc(open_file);
            switch (c)
            {
            case EOF:
                stopmsg(stopmsg_flags::FIXED_FONT, ParseErrs(PE_UNEXPECTED_EOF));
                return false;
            case '\r':
            case '\n':
                stopmsg(stopmsg_flags::FIXED_FONT, ParseErrs(PE_NO_LEFT_BRACKET_FIRST_LINE));
                return false;
            case '{':
                done = true;
                break;
            default :
                break;
            }
        }
    }
    return true;
}

/* This function sets the
    symmetry and converts a formula into a string  with no spaces,
    and one comma after each expression except where the ':' is placed
    and except the final expression in the formula. The open file passed
    as an argument is open in "rb" mode and is positioned at the first
    letter of the name of the formula to be prepared. This function
    is called from run_formula() below.
*/
static std::string PrepareFormula(std::FILE *file, bool report_bad_sym)
{
    long filepos = ftell(file);

    //Test for a repeat
    if (!frm_check_name_and_sym(file, report_bad_sym))
    {
        std::fseek(file, filepos, SEEK_SET);
        return {};
    }
    if (!frm_prescan(file))
    {
        std::fseek(file, filepos, SEEK_SET);
        return {};
    }

    if (s_chars_in_formula > 8190)
    {
        std::fseek(file, filepos, SEEK_SET);
        return {};
    }

    std::FILE *debug_fp = nullptr;
    if (g_debug_flag == debug_flags::write_formula_debug_information)
    {
        debug_fp = std::fopen("debugfrm.txt", "at");
        if (debug_fp != nullptr)
        {
            std::fprintf(debug_fp, "%s\n", g_formula_name.c_str());
            if (g_symmetry != symmetry_type::NONE)
            {
                auto it = std::find_if(std::begin(s_symmetry_names), std::end(s_symmetry_names),
                    [](const SymmetryName &item) { return item.n == g_symmetry; });
                if (it != std::end(s_symmetry_names))
                {
                    std::fprintf(debug_fp, "%s\n", it->s);
                }
            }
        }
    }

    std::string FormulaStr;

    bool Done = false;

    //skip opening end-of-lines
    token_st temp_tok;
    while (!Done)
    {
        frmgettoken(file, &temp_tok);
        if (temp_tok.type == token_type::NOT_A_TOKEN)
        {
            stopmsg(stopmsg_flags::FIXED_FONT, "Unexpected token error in PrepareFormula\n");
            std::fseek(file, filepos, SEEK_SET);
            if (debug_fp != nullptr)
            {
                std::fclose(debug_fp);
            }
            return {};
        }
        if (temp_tok.type == token_type::END_OF_FORMULA)
        {
            stopmsg(stopmsg_flags::FIXED_FONT, "Formula has no executable instructions\n");
            std::fseek(file, filepos, SEEK_SET);
            if (debug_fp != nullptr)
            {
                std::fclose(debug_fp);
            }
            return {};
        }
        if (temp_tok.str[0] != ',')
        {
            FormulaStr += temp_tok.str;
            Done = true;
        }
    }

    Done = false;
    while (!Done)
    {
        frmgettoken(file, &temp_tok);
        switch (temp_tok.type)
        {
        case token_type::NOT_A_TOKEN:
            stopmsg(stopmsg_flags::FIXED_FONT, "Unexpected token error in PrepareFormula\n");
            std::fseek(file, filepos, SEEK_SET);
            if (debug_fp != nullptr)
            {
                std::fclose(debug_fp);
            }
            return {};
        case token_type::END_OF_FORMULA:
            Done = true;
            std::fseek(file, filepos, SEEK_SET);
            break;
        default:
            FormulaStr += temp_tok.str;
            break;
        }
    }

    if (debug_fp != nullptr && !FormulaStr.empty())
    {
        std::fprintf(debug_fp, "   %s\n", FormulaStr.c_str());
    }
    if (debug_fp != nullptr)
    {
        std::fclose(debug_fp);
    }

    return FormulaStr;
}

int BadFormula()
{
    //  this is called when a formula is bad, instead of calling
    //     the normal functions which will produce undefined results
    return 1;
}

//  returns true if an error occurred
bool run_formula(const std::string &name, bool report_bad_sym)
{
    std::FILE *entry_file = nullptr;

    //  first set the pointers so they point to a fn which always returns 1
    g_cur_fractal_specific->per_pixel = BadFormula;
    g_cur_fractal_specific->orbitcalc = BadFormula;

    if (g_formula_name.empty())
    {
        return true;  //  and don't reset the pointers
    }

    if (find_file_item(g_formula_filename, name.c_str(), &entry_file, gfe_type::FORMULA))
    {
        stopmsg(ParseErrs(PE_COULD_NOT_OPEN_FILE_WHERE_FORMULA_LOCATED));
        return true;
    }

    s_formula = PrepareFormula(entry_file, report_bad_sym);
    std::fclose(entry_file);

    if (!s_formula.empty())  //  No errors while making string
    {
        parser_allocate();  //  ParseStr() will test if this alloc worked
        if (parse_formula_text(s_formula.c_str()))
        {
            return true;   //  parse failed, don't change fn pointers
        }
        if (g_uses_jump && fill_jump_struct())
        {
            stopmsg(ParseErrs(PE_ERROR_IN_PARSING_JUMP_STATEMENTS));
            return true;
        }

        // all parses succeeded so set the pointers back to good functions
        g_cur_fractal_specific->per_pixel = form_per_pixel;
        g_cur_fractal_specific->orbitcalc = Formula;
        return false;
    }
    return true; // error in making string
}

bool fpFormulaSetup()
{
    s_math_type = D_MATH;
    return !run_formula(g_formula_name, false); // run_formula() returns true for failure
}

#undef FORMULA_INTEGER_MATH
bool intFormulaSetup()
{
#ifndef FORMULA_INTEGER_MATH
    static bool been_here = false;
    if (!been_here)
    {
        stopmsg("This integer fractal type is unimplemented;\n"
                "Use float=yes to get a real image.");
        been_here = true;
    }
    return false;
#else
    s_math_type = L_MATH;
    s_fudge = (double)(1L << g_bit_shift);
    g_fudge_limit = (double)0x7fffffffL / s_fudge;
    s_shift_back = 32 - g_bit_shift;
    return !run_formula(g_formula_name, false);
#endif
}

void init_misc()
{
    static Arg argfirst;
    static Arg argsecond;
    if (v.empty())
    {
        v.resize(5);
    }
    Arg1 = &argfirst;
    Arg2 = &argsecond; // needed by all the ?Stk* functions
    s_fudge = (double)(1L << g_bit_shift);
    g_fudge_limit = (double)0x7fffffffL / s_fudge;
    s_shift_back = 32 - g_bit_shift;
    s_delta16 = g_bit_shift - 16;
    g_bit_shift_less_1 = g_bit_shift-1;
    g_frm_uses_p1 = false;
    g_frm_uses_p2 = false;
    g_frm_uses_p3 = false;
    g_uses_jump = false;
    g_frm_uses_ismand = false;
    g_frm_uses_p4 = false;
    g_frm_uses_p5 = false;
}

static void parser_allocate()
{
    free_workarea();
    g_max_function_ops = 2300;
    g_max_function_args = (unsigned)(g_max_function_ops/2.5);

    f.reserve(g_max_function_ops);
    Store.resize(MAX_STORES);
    Load.resize(MAX_LOADS);
    v.resize(g_max_function_args);

    if (!parse_formula_text(s_formula.c_str()))
    {
        // per Chuck Ebbert, fudge these up a little
        g_max_function_ops = g_operation_index + 4;
        g_max_function_args = g_variable_index + 4;
    }
    g_frm_uses_p1 = false;
    g_frm_uses_p2 = false;
    g_frm_uses_p3 = false;
    g_frm_uses_p4 = false;
    g_frm_uses_p5 = false;
}

void free_workarea()
{
    Store.clear();
    Load.clear();
    v.clear();
    f.clear();
}

struct error_data_st
{
    long start_pos;
    long error_pos;
    int error_number;
};

static std::array<error_data_st, 3> s_errors{};

static void frm_error(std::FILE * open_file, long begin_frm)
{
    token_st tok;
    int chars_to_error = 0;
    int chars_in_error = 0;
    int token_count;
    int statement_len;
    int line_number;
    char msgbuf[900];
    long filepos;
    std::strcpy(msgbuf, "\n");

    for (int j = 0; j < 3 && s_errors[j].start_pos; j++)
    {
        bool const initialization_error = s_errors[j].error_number == PE_SECOND_COLON;
        std::fseek(open_file, begin_frm, SEEK_SET);
        line_number = 1;
        while (ftell(open_file) != s_errors[j].error_pos)
        {
            int i = fgetc(open_file);
            if (i == '\n')
            {
                line_number++;
            }
            else if (i == EOF || i == '}')
            {
                stopmsg("Unexpected EOF or end-of-formula in error function.\n");
                std::fseek(open_file, s_errors[j].error_pos, SEEK_SET);
                frmgettoken(open_file, &tok); //reset file to end of error token
                return;
            }
        }
        std::sprintf(&msgbuf[(int) std::strlen(msgbuf)], "Error(%d) at line %d:  %s\n  ", s_errors[j].error_number, line_number, ParseErrs(s_errors[j].error_number));
        int i = (int) std::strlen(msgbuf);
        std::fseek(open_file, s_errors[j].start_pos, SEEK_SET);
        token_count = 0;
        statement_len = token_count;
        bool done = false;
        while (!done)
        {
            filepos = ftell(open_file);
            if (filepos == s_errors[j].error_pos)
            {
                chars_to_error = statement_len;
                frmgettoken(open_file, &tok);
                chars_in_error = (int) std::strlen(tok.str);
                statement_len += chars_in_error;
                token_count++;
            }
            else
            {
                frmgettoken(open_file, &tok);
                statement_len += (int) std::strlen(tok.str);
                token_count++;
            }
            if (tok.type == token_type::END_OF_FORMULA ||
                (tok.type == token_type::OPERATOR && (tok.id == token_id::OP_COMMA || tok.id == token_id::OP_COLON)) ||
                (tok.type == token_type::NOT_A_TOKEN && tok.id == token_id::END_OF_FILE))
            {
                done = true;
                if (token_count > 1 && !initialization_error)
                {
                    token_count--;
                }
            }
        }
        std::fseek(open_file, s_errors[j].start_pos, SEEK_SET);
        if (chars_in_error < 74)
        {
            while (chars_to_error + chars_in_error > 74)
            {
                frmgettoken(open_file, &tok);
                chars_to_error -= (int) std::strlen(tok.str);
                token_count--;
            }
        }
        else
        {
            std::fseek(open_file, s_errors[j].error_pos, SEEK_SET);
            chars_to_error = 0;
            token_count = 1;
        }
        while ((int) std::strlen(&msgbuf[i]) <=74 && token_count--)
        {
            frmgettoken(open_file, &tok);
            std::strcat(msgbuf, tok.str);
        }
        std::fseek(open_file, s_errors[j].error_pos, SEEK_SET);
        frmgettoken(open_file, &tok);
        if ((int) std::strlen(&msgbuf[i]) > 74)
        {
            msgbuf[i + 74] = (char) 0;
        }
        std::strcat(msgbuf, "\n");
        i = (int) std::strlen(msgbuf);
        while (chars_to_error-- > -2)
        {
            std::strcat(msgbuf, " ");
        }
        if (s_errors[j].error_number == PE_TOKEN_TOO_LONG)
        {
            chars_in_error = 33;
        }
        while (chars_in_error-- && (int) std::strlen(&msgbuf[i]) <=74)
        {
            std::strcat(msgbuf, "^");
        }
        std::strcat(msgbuf, "\n");
    }
    stopmsg(stopmsg_flags::FIXED_FONT, msgbuf);
}

/*frm_prescan() takes an open file with the file pointer positioned at
  the beginning of the relevant formula, and parses the formula, token
  by token, for syntax errors. The function also accumulates data for
  memory allocation to be done later.

  The function returns 1 if success, and 0 if errors are found.
*/
static bool frm_prescan(std::FILE * open_file)
{
    long filepos;
    long statement_pos;
    long orig_pos;
    bool done = false;
    token_st this_token;
    int errors_found = 0;
    bool ExpectingArg = true;
    bool NewStatement = true;
    bool assignment_ok = true;
    bool already_got_colon = false;
    unsigned long else_has_been_used = 0;
    unsigned long waiting_for_mod = 0;
    int waiting_for_endif = 0;
    int max_parens = sizeof(long) * 8;

    s_num_jumps = 0UL;
    s_num_stores = 0UL;
    s_num_loads = 0UL;
    s_num_ops = 0UL;
    s_chars_in_formula = 0U;
    g_uses_jump = false;
    s_paren = 0;

    statement_pos = ftell(open_file);
    orig_pos = statement_pos;
    for (error_data_st &error : s_errors)
    {
        error.start_pos    = 0L;
        error.error_pos    = 0L;
        error.error_number = 0;
    }

    while (!done)
    {
        filepos = ftell(open_file);
        frmgettoken(open_file, &this_token);
        s_chars_in_formula += (int) std::strlen(this_token.str);
        switch (this_token.type)
        {
        case token_type::NOT_A_TOKEN:
            assignment_ok = false;
            switch (this_token.id)
            {
            case token_id::END_OF_FILE:
                stopmsg(ParseErrs(PE_UNEXPECTED_EOF));
                std::fseek(open_file, orig_pos, SEEK_SET);
                return false;
            case token_id::ILLEGAL_CHARACTER:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_ILLEGAL_CHAR;
                }
                break;
            case token_id::ILLEGAL_VARIABLE_NAME:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_ILLEGAL_VAR_NAME;
                }
                break;
            case token_id::TOKEN_TOO_LONG:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_TOKEN_TOO_LONG;
                }
                break;
            case token_id::FUNC_USED_AS_VAR:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_FUNC_USED_AS_VAR;
                }
                break;
            case token_id::JUMP_MISSING_BOOLEAN:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_JUMP_NEEDS_BOOLEAN;
                }
                break;
            case token_id::JUMP_WITH_ILLEGAL_CHAR:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_NO_CHAR_AFTER_THIS_JUMP;
                }
                break;
            case token_id::UNDEFINED_FUNCTION:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_UNDEFINED_FUNCTION;
                }
                break;
            case token_id::ILLEGAL_OPERATOR:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_UNDEFINED_OPERATOR;
                }
                break;
            case token_id::ILL_FORMED_CONSTANT:
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_INVALID_CONST;
                }
                break;
            default:
                stopmsg("Unexpected arrival at default case in prescan()");
                std::fseek(open_file, orig_pos, SEEK_SET);
                return false;
            }
            break;
        case token_type::PARENS:
            assignment_ok = false;
            NewStatement = false;
            switch (this_token.id)
            {
            case token_id::OPEN_PARENS:
                if (++s_paren > max_parens)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_NESTING_TO_DEEP;
                    }
                }
                else if (!ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                    }
                }
                waiting_for_mod = waiting_for_mod << 1;
                break;
            case token_id::CLOSE_PARENS:
                if (s_paren)
                {
                    s_paren--;
                }
                else
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_NEED_A_MATCHING_OPEN_PARENS;
                    }
                    s_paren = 0;
                }
                if (waiting_for_mod & 1L)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_UNMATCHED_MODULUS;
                    }
                }
                else
                {
                    waiting_for_mod = waiting_for_mod >> 1;
                }
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                break;
            default:
                break;
            }
            break;
        case token_type::PARAM_VARIABLE: //i.e. p1, p2, p3, p4 or p5
            s_num_ops++;
            s_num_loads++;
            NewStatement = false;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            ExpectingArg = false;
            break;
        case token_type::USER_NAMED_VARIABLE: // i.e. c, iter, etc.
            s_num_ops++;
            s_num_loads++;
            NewStatement = false;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            ExpectingArg = false;
            break;
        case token_type::PREDEFINED_VARIABLE: // i.e. z, pixel, whitesq, etc.
            s_num_ops++;
            s_num_loads++;
            NewStatement = false;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            ExpectingArg = false;
            break;
        case token_type::REAL_CONSTANT: // i.e. 4, (4,0), etc.
            assignment_ok = false;
            s_num_ops++;
            s_num_loads++;
            NewStatement = false;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            ExpectingArg = false;
            break;
        case token_type::COMPLEX_CONSTANT: // i.e. (1,2) etc.
            assignment_ok = false;
            s_num_ops++;
            s_num_loads++;
            NewStatement = false;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            ExpectingArg = false;
            break;
        case token_type::FUNCTION:
            assignment_ok = false;
            NewStatement = false;
            s_num_ops++;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            break;
        case token_type::PARAM_FUNCTION:
            assignment_ok = false;
            s_num_ops++;
            if (!ExpectingArg)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                }
            }
            NewStatement = false;
            break;
        case token_type::FLOW_CONTROL:
            assignment_ok = false;
            s_num_ops++;
            s_num_jumps++;
            if (!NewStatement)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_JUMP_NOT_FIRST;
                }
            }
            else
            {
                g_uses_jump = true;
                switch (this_token.id)
                {
                case token_id::JUMP_IF:  // if
                    else_has_been_used = else_has_been_used << 1;
                    waiting_for_endif++;
                    break;
                case token_id::JUMP_ELSE_IF: //ELSEIF
                    s_num_ops += 3; //else + two clear statements
                    s_num_jumps++;  // this involves two jumps
                    if (else_has_been_used % 2)
                    {
                        if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                        {
                            s_errors[errors_found].start_pos      = statement_pos;
                            s_errors[errors_found].error_pos      = filepos;
                            s_errors[errors_found++].error_number = PE_ENDIF_REQUIRED_AFTER_ELSE;
                        }
                    }
                    else if (!waiting_for_endif)
                    {
                        if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                        {
                            s_errors[errors_found].start_pos      = statement_pos;
                            s_errors[errors_found].error_pos      = filepos;
                            s_errors[errors_found++].error_number = PE_MISPLACED_ELSE_OR_ELSEIF;
                        }
                    }
                    break;
                case token_id::JUMP_ELSE: //ELSE
                    if (else_has_been_used % 2)
                    {
                        if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                        {
                            s_errors[errors_found].start_pos      = statement_pos;
                            s_errors[errors_found].error_pos      = filepos;
                            s_errors[errors_found++].error_number = PE_ENDIF_REQUIRED_AFTER_ELSE;
                        }
                    }
                    else if (!waiting_for_endif)
                    {
                        if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                        {
                            s_errors[errors_found].start_pos      = statement_pos;
                            s_errors[errors_found].error_pos      = filepos;
                            s_errors[errors_found++].error_number = PE_MISPLACED_ELSE_OR_ELSEIF;
                        }
                    }
                    else_has_been_used = else_has_been_used | 1;
                    break;
                case token_id::JUMP_END_IF: //ENDIF
                    else_has_been_used = else_has_been_used >> 1;
                    waiting_for_endif--;
                    if (waiting_for_endif < 0)
                    {
                        if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                        {
                            s_errors[errors_found].start_pos      = statement_pos;
                            s_errors[errors_found].error_pos      = filepos;
                            s_errors[errors_found++].error_number = PE_ENDIF_WITH_NO_IF;
                        }
                        waiting_for_endif = 0;
                    }
                    break;
                default:
                    break;
                }
            }
            break;
        case token_type::OPERATOR:
            s_num_ops++; //This will be corrected below in certain cases
            switch (this_token.id)
            {
            case token_id::OP_COMMA:
            case token_id::OP_COLON:    // end of statement and :
                s_num_ops++; // ParseStr inserts a dummy op
                if (s_paren)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_NEED_MORE_CLOSE_PARENS;
                    }
                    s_paren = 0;
                }
                if (waiting_for_mod)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_UNMATCHED_MODULUS;
                    }
                    waiting_for_mod = 0;
                }
                if (!ExpectingArg)
                {
                    if (this_token.id == token_id::OP_COLON)
                    {
                        s_num_ops += 2;
                    }
                    else
                    {
                        s_num_ops++;
                    }
                }
                else if (!NewStatement)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                if (this_token.id == token_id::OP_COLON && waiting_for_endif)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_UNMATCHED_IF_IN_INIT_SECTION;
                    }
                    waiting_for_endif = 0;
                }
                if (this_token.id == token_id::OP_COLON && already_got_colon)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SECOND_COLON;
                    }
                }
                if (this_token.id == token_id::OP_COLON)
                {
                    already_got_colon = true;
                }
                NewStatement = true;
                ExpectingArg = true;
                assignment_ok = true;
                statement_pos = ftell(open_file);
                break;
            case token_id::OP_NOT_EQUAL:     // !=
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_ASSIGN:     // =
                s_num_ops--; //this just converts a load to a store
                s_num_loads--;
                s_num_stores++;
                if (!assignment_ok)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_ILLEGAL_ASSIGNMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_EQUAL:     // ==
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_LT:     // <
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_LE:     // <=
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_GT:     // >
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_GE:     // >=
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_MODULUS:     // | (half of the modulus operator)
                assignment_ok = false;
                if (!(waiting_for_mod & 1L))
                {
                    s_num_ops--;
                }
                if (!(waiting_for_mod & 1L) && !ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_OPERATOR;
                    }
                }
                else if ((waiting_for_mod & 1L) && ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                waiting_for_mod = waiting_for_mod ^ 1L; //switch right bit
                break;
            case token_id::OP_OR:     // ||
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_AND:    // &&
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_PLUS:    // + case 11 (":") is up with case 0
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_MINUS:    // -
                assignment_ok = false;
                ExpectingArg = true;
                break;
            case token_id::OP_MULTIPLY:    // *
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_DIVIDE:    // /
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                ExpectingArg = true;
                break;
            case token_id::OP_POWER:    // ^
                assignment_ok = false;
                if (ExpectingArg)
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                    }
                }
                filepos = ftell(open_file);
                frmgettoken(open_file, &this_token);
                if (this_token.str[0] == '-')
                {
                    if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                    {
                        s_errors[errors_found].start_pos      = statement_pos;
                        s_errors[errors_found].error_pos      = filepos;
                        s_errors[errors_found++].error_number = PE_NO_NEG_AFTER_EXPONENT;
                    }
                }
                else
                {
                    std::fseek(open_file, filepos, SEEK_SET);
                }
                ExpectingArg = true;
                break;
            default:
                break;
            }
            break;
        case token_type::END_OF_FORMULA:
            s_num_ops += 3; // Just need one, but a couple of extra just for the heck of it
            if (s_paren)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_NEED_MORE_CLOSE_PARENS;
                }
                s_paren = 0;
            }
            if (waiting_for_mod)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_UNMATCHED_MODULUS;
                }
                waiting_for_mod = 0;
            }
            if (waiting_for_endif)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_IF_WITH_NO_ENDIF;
                }
                waiting_for_endif = 0;
            }
            if (ExpectingArg && !NewStatement)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_SHOULD_BE_ARGUMENT;
                }
                statement_pos = ftell(open_file);
            }

            if (s_num_jumps >= MAX_JUMPS)
            {
                if (!errors_found || s_errors[errors_found-1].start_pos != statement_pos)
                {
                    s_errors[errors_found].start_pos      = statement_pos;
                    s_errors[errors_found].error_pos      = filepos;
                    s_errors[errors_found++].error_number = PE_TOO_MANY_JUMPS;
                }
            }
            done = true;
            break;
        default:
            break;
        }
        if (errors_found == 3)
        {
            done = true;
        }
    }
    if (s_errors[0].start_pos)
    {
        frm_error(open_file, orig_pos);
        std::fseek(open_file, orig_pos, SEEK_SET);
        return false;
    }
    std::fseek(open_file, orig_pos, SEEK_SET);

    return true;
}
