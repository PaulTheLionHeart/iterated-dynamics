#pragma once

#include "big.h"
#include "calcfrac.h"
#include "fractype.h"
#include "helpdefs.h"

struct AlternateMath
{
    fractal_type type;                  // index in fractalname of the fractal
    bf_math_type math;                  // kind of math used
    int (*orbitcalc)();                 // function that calculates one orbit
    int (*per_pixel)();                 // once-per-pixel init
    bool (*per_image)();                // once-per-image setup
};

struct MOREPARAMS
{
    fractal_type type;                      // index in fractalname of the fractal
    char const *param[MAX_PARAMS-4];     // name of the parameters
    double   paramvalue[MAX_PARAMS-4];   // default parameter values
};

struct fractalspecificstuff
{
    char const  *name;                  // name of the fractal
                                        // (leading "*" supresses name display)
    char const  *param[4];              // name of the parameters
    double paramvalue[4];               // default parameter values
    help_labels helptext;               // helpdefs.h HT_xxxx, -1 for none
    help_labels helpformula;            // helpdefs.h HF_xxxx, -1 for none
    unsigned flags;                     // constraints, bits defined below
    float xmin;                         // default XMIN corner
    float xmax;                         // default XMAX corner
    float ymin;                         // default YMIN corner
    float ymax;                         // default YMAX corner
    int   isinteger;                    // 1 if integerfractal, 0 otherwise
    fractal_type tojulia;               // mandel-to-julia switch
    fractal_type tomandel;              // julia-to-mandel switch
    fractal_type tofloat;               // integer-to-floating switch
    symmetry_type symmetry;             // applicable symmetry logic
                                        //  0 = no symmetry
                                        // -1 = y-axis symmetry (If No Params)
                                        //  1 = y-axis symmetry
                                        // -2 = x-axis symmetry (No Parms)
                                        //  2 = x-axis symmetry
                                        // -3 = y-axis AND x-axis (No Parms)
                                        //  3 = y-axis AND x-axis symmetry
                                        // -4 = polar symmetry (No Parms)
                                        //  4 = polar symmetry
                                        //  5 = PI (sin/cos) symmetry
                                        //  6 = NEWTON (power) symmetry
                                        //
    int (*orbitcalc)();                 // function that calculates one orbit
    int (*per_pixel)();                 // once-per-pixel init
    bool (*per_image)();                // once-per-image setup
    int (*calctype)();                  // name of main fractal function
    int orbit_bailout;                  // usual bailout value for orbit calc
};

// bitmask values for fractalspecific flags
enum
{
    NOZOOM = 1,      // zoombox not allowed at all
    NOGUESS = 2,     // solid guessing not allowed
    NOTRACE = 4,     // boundary tracing not allowed
    NOROTATE = 8,    // zoombox rotate/stretch not allowed
    NORESUME = 16,   // can't interrupt and resume
    INFCALC = 32,    // this type calculates forever
    TRIG1 = 64,      // number of trig functions in formula
    TRIG2 = 128,     //
    TRIG3 = 192,     //
    TRIG4 = 256,     //
    PARMS3D = 1024,  // uses 3d parameters
    OKJB = 2048,     // works with Julibrot
    MORE = 4096,     // more than 4 parms
    BAILTEST = 8192, // can use different bailout tests
    BF_MATH = 16384, // supports arbitrary precision
    LD_MATH = 32768  // supports long double
};

extern AlternateMath         g_alternate_math[];    // alternate math function pointers
extern fractalspecificstuff  g_fractal_specific[];
extern MOREPARAMS            g_more_fractal_params[];
extern int                   g_num_fractal_types;

inline bool per_image()
{
    return g_fractal_specific[static_cast<int>(g_fractal_type)].per_image();
}
inline int per_pixel()
{
    return g_fractal_specific[static_cast<int>(g_fractal_type)].per_pixel();
}
inline int orbit_calc()
{
    return g_fractal_specific[static_cast<int>(g_fractal_type)].orbitcalc();
}
