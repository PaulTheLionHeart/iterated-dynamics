/****************************************************
    PERTURBATION.CPP a module to explore Perturbation
    Thanks to Shirom Makkad fractaltodesktop@gmail.com
    Written in Microsoft Visual C++ by Paul de Leeuw.
****************************************************/

#include "fractalp.h"
#include "PertEngine.h"
#include "drivers.h"
#include "id_data.h"

// some early tests
//const char* xZoomPointString = "-1.25736802846652839265383159384773654166836713857126000896912753375688559878664765114255696457015368246531973104439755978333044015506759938503739206829441575363669402497147343368904702066174408250247081855416385744218741909521990441308969603994513271641284298225323509381146075334937409491188";
//const char* yZoomPointString = "0.37873083102862491151257052392242106932532193327534173605649141946411779667848532042309666819671311329800095959763206221251386819369531602358854394169140220049675504811341950346171196600590463661845947574424944950533273158558278821586333530950155398785389980291835095955110139525825547062070";
//const char* xZoomPointString = "-1.9428";
//const char* yZoomPointString = "0.0";
const char* xZoomPointString = "-0.75";
const char* yZoomPointString = "0.01";
//const char* xZoomPointString = "-2.0";
//const char* yZoomPointString = "0.0";
//const char* xZoomPointString = "-1.90750420816369968453818946601049";
//const char* yZoomPointString = "-8.007949072567580760039e-12";
//const char* xZoomPointString = "-1.74997338996683218759081143886";
//const char* yZoomPointString = "9.33818690200204924e-14";
//const char* xZoomPointString = "0.4114706336593580722758780101609024";
//const char* yZoomPointString = "0.6074067173514734090007173476657398";
//const char* xZoomPointString = "-1.6129965373041952026071917970859685";
//const char* yZoomPointString = "1.70637685183811519151570457965e-6";

//extern  int     driver_key_pressed();
extern  void    (*g_plot)(int, int, int); // function pointer
extern	int	    potential(double, long);
extern  int     g_screen_x_dots, g_screen_y_dots; // # of dots on the physical screen
extern calc_status_value g_calc_status;

extern	long	g_max_iterations;
extern  double  g_x_max, g_x_min, g_y_max, g_y_min;
extern  bool    g_potential_flag;

extern  double  g_magnitude_limit;          // bailout level
int decimals = /*bflength * 4*/ 24;         // we can sort this out later
extern	double	g_params[];
extern	int	    g_outside_color;			// outside filters
extern	int	    g_inside_color;			    // inside  filters
extern	int	    g_biomorph;			        // biomorph colour
//extern	RGBTRIPLE FilterRGB;			// for Tierazon filters
//extern	double	dStrands;
//extern	BOOL	UseCurrentPalette;		// do we use the ManpWIN palette? If false, generate internal filter palette

//////////////////////////////////////////////////////////////////////////////////////////////////////
extern	HWND	GlobalHwnd;			// This is the main windows handle
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*
static	CTZfilter	TZfilter;		// Tierazon filters
*/
CPertEngine PertEngine;
char		PertStatus[200];

int         DoPerturbation(int);
void        BigCvtcentermag(bf_t *Xctr, bf_t *Yctr, double *Magnification, double *Xmagfactor, double *Rotation, double *Skew);
void        bf2BigNum(BigDouble *BigNum, bf_t bfNum);
extern void cvtcentermag(double *Xctr, double *Yctr, LDBL *Magnification, double *Xmagfactor, double *Rotation, double *Skew);


/*************************************************************************
    Format string derived from a Bignum
    mpf_get_str() generates strings without decimal points and gives the exponent
    so we need to format it as a normal number
*************************************************************************/

void ConvertBignum2String(char *s, mpfr_t num)                                // can be used for debugging
    {
    char FormatString[24];

    sprintf(FormatString, "%%.%dRf", decimals + PRECISION_FACTOR);
    mpfr_sprintf(s, FormatString, num);
    }

/**************************************************************************
	Initialise Perturbation engine
**************************************************************************/

bool	InitPerturbation(int subtype)
    {
    double  mandel_width;    // width of display
    double  xCentre, yCentre, Xmagfactor, Rotation, Skew;
    char s[1200];
    bf_math_type TempMathType = g_bf_math;

    g_bf_math = bf_math_type::BIGFLT;       // pert engine needs this (for the moment)

/*
    int bitcount = decimals * 5;
    if (bitcount < 30)
        bitcount = 30;
    if (bitcount > SIZEOF_BF_VARS - 10)
        bitcount = SIZEOF_BF_VARS - 10;

//    mpfr_set_default_prec(bitcount);
        mpfr_set_default_prec(1600);
        xBigCentre.ChangePrecision(1600);
        yBigCentre.ChangePrecision(1600);
    */
    bf_t xBigCentre, yBigCentre, BigTmp;
    int saved;
    saved = save_stack();
    xBigCentre = alloc_stack(g_bf_length + 2);
    yBigCentre = alloc_stack(g_bf_length + 2);
    BigTmp = alloc_stack(g_bf_length + 2);

    double  Magnification;


    if (TempMathType != bf_math_type::NONE) // we assume bignum is flagged and bf variables are initialised
        {
        BigCvtcentermag(&xBigCentre, &yBigCentre, &Magnification, &Xmagfactor, &Rotation, &Skew);              // have to make a mpfr version of this
        }
    else
        {
        LDBL LDMagnification;
        cvtcentermag(&xCentre, &yCentre, &LDMagnification, &Xmagfactor, &Rotation, &Skew);
        floattobf(xBigCentre, xCentre);
        floattobf(yBigCentre, yCentre);
        }

    if (TempMathType == bf_math_type::NONE) 
        mandel_width = g_y_max - g_y_min;
    else
        {
        sub_bf(BigTmp, g_bf_y_max, g_bf_y_min);
        mandel_width = bftofloat(BigTmp);
        }

    /*
    if (method >= TIERAZONFILTERS)
	    TZfilter.InitFilter(method, threshold, dStrands, &FilterRGB, UseCurrentPalette);		// initialise the constants used by Tierazon fractals
    if (BigNumFlag)
	    mandel_width = mpfr_get_d(BigWidth.x, MPFR_RNDN);
    */
    neg_bf(yBigCentre, yBigCentre);
    PertEngine.initialiseCalculateFrame(g_screen_x_dots, g_screen_y_dots, g_max_iterations, xBigCentre, yBigCentre, mandel_width / 2, g_potential_flag/*, &TZfilter*/);
    restore_stack(saved);
    g_bf_math = TempMathType;               // better restore it before we return
    DoPerturbation(subtype);
    g_calc_status = calc_status_value::COMPLETED;
    return false;
    }

/**************************************************************************
	The Perturbation engine
**************************************************************************/

int DoPerturbation(int subtype)
    {
    int (*UserData)() = driver_key_pressed;
    Complex a = {0, 0};             // future TheRedshiftRider
    bool    IsPositive = false;
    int     degree;  // power
//    BYTE    subtype; // subtype

//    subtype = (int) g_params[0];
    degree = (int) g_params[2];
/*        
    if (subtype == 53)              // future TheRedshiftRider
	    {
	    a.x = g_params[2];
	    a.y = g_params[3];
	    IsPositive = (g_params[4] == 1.0);
	    }
*/
    if (PertEngine.calculateOneFrame(g_magnitude_limit, PertStatus, degree, g_inside_color, g_outside_color, g_biomorph, subtype, a, IsPositive, UserData, g_plot, potential/*, &TZfilter, &TrueCol*/) < 0)
        //    if (frameCalculator.calculateOneFrame(rqlim, PertStatus, degree, method, biomorph, subtype, a, IsPositive, UserData, plot, potential) < 0)
	    return -1;
    return 0;
    }

/*************************************************************************
    Format Bignum derived from a string
    mpf_set_str() Bignum from strings in the format:
    MeN, where M = mantissa and N exponent
*************************************************************************/

void	ConvertString2Bignum(mpfr_t num, char *s)
    {
    mpfr_set_str(num, s, 10, MPFR_RNDN);
    }

/**************************************************************************
	Convert bf to Bignum
**************************************************************************/

void bf2BigNum(BigDouble *BigNum, bf_t bfNum)
    {
    // let's just do a quick and dirty for now going via text.
    int     bfLengthNeeded = strlen_needed_bf();
    char    *bigstr = NULL;
    int     dec = g_decimals;

    bigstr = new char[bfLengthNeeded];

    bftostr(bigstr, dec, bfNum);
    ConvertString2Bignum(BigNum->x, bigstr);

    if (bigstr)  {delete[] bigstr; bigstr = NULL;}
    }
    
/**************************************************************************
	Convert Bignum to bf
**************************************************************************/

void BigNum2bf(bf_t *bfNum, BigDouble BigNum)
    {
    // let's just do a quick and dirty for now going via text.
    char    *bigstr = NULL;
    int     dec = g_decimals;

    bigstr = new char[SIZEOF_BF_VARS];
    ConvertBignum2String(bigstr, BigNum.x);
    strtobf(*bfNum, bigstr);
    if (bigstr)  {delete[] bigstr; bigstr = NULL;}
    }
    
/**************************************************************************
	Convert corners to centre/mag using BigNum
**************************************************************************/

void BigCvtcentermag(bf_t *Xctr, bf_t *Yctr, double *Magnification, double *Xmagfactor, double *Rotation, double *Skew)
    {
    // needs to be LDBL or won't work past 307 (-DBL_MIN_10_EXP) or so digits
    double      Height;
    bf_t BigWidth;
    bf_t BigHeight;
    bf_t BigTmpx;

    int saved;
    saved = save_stack();
    BigWidth = alloc_stack(g_bf_length + 2);
    BigHeight = alloc_stack(g_bf_length + 2);
    BigTmpx = alloc_stack(g_bf_length + 2);

    // simple normal case first
    // if (g_x_3rd == g_x_min && g_y_3rd == g_y_min)
    if (!cmp_bf(g_bf_x_3rd, g_bf_x_min) && !cmp_bf(g_bf_y_3rd, g_bf_y_min))
        {
        // no rotation or skewing, but stretching is allowed
        // Width  = g_x_max - g_x_min;
        sub_bf(BigWidth, g_bf_x_max, g_bf_x_min);
        double Width = bftofloat(BigWidth);
        // Height = g_y_max - g_y_min;
        sub_bf(BigHeight, g_bf_y_max, g_bf_y_min);
        Height = bftofloat(BigHeight);
        // *Xctr = (g_x_min + g_x_max)/2;
        add_bf(BigTmpx, g_bf_x_max, g_bf_x_min);
        half_bf(*Xctr, BigTmpx);
        // *Yctr = (g_y_min + g_y_max)/2;
        add_bf(BigTmpx, g_bf_y_max, g_bf_y_min);
        half_bf(*Yctr, BigTmpx);
        *Magnification = 2 / Height;
        *Xmagfactor = (double) (Height / (DEFAULT_ASPECT * Width));
        *Rotation = 0.0;
        *Skew = 0.0;
        }
/*              // something for later...
        else
        {
            bftmpx = alloc_stack(bflength + 2);
            bf_t bftmpy = alloc_stack(bflength + 2);

            // set up triangle ABC, having sides abc
            // side a = bottom, b = left, c = diagonal not containing (x3rd,y3rd)
            // IMPORTANT: convert from bf AFTER subtracting

            // tmpx = g_x_max - g_x_min;
            sub_bf(bftmpx, g_bf_x_max, g_bf_x_min);
            LDBL tmpx1 = bftofloat(bftmpx);
            // tmpy = g_y_max - g_y_min;
            sub_bf(bftmpy, g_bf_y_max, g_bf_y_min);
            LDBL tmpy1 = bftofloat(bftmpy);
            LDBL c2 = tmpx1 * tmpx1 + tmpy1 * tmpy1;

            // tmpx = g_x_max - g_x_3rd;
            sub_bf(bftmpx, g_bf_x_max, g_bf_x_3rd);
            tmpx1 = bftofloat(bftmpx);

            // tmpy = g_y_min - g_y_3rd;
            sub_bf(bftmpy, g_bf_y_min, g_bf_y_3rd);
            tmpy1 = bftofloat(bftmpy);
            LDBL a2 = tmpx1 * tmpx1 + tmpy1 * tmpy1;
            LDBL a = sqrtl(a2);

            // divide tmpx and tmpy by |tmpx| so that double version of atan2() can be used
            // atan2() only depends on the ratio, this puts it in double's range
            int signx = sign(tmpx1);
            LDBL tmpy = 0.0;
            if (signx)
            {
                tmpy = tmpy1 / tmpx1 * signx; // tmpy = tmpy / |tmpx|
            }
            *Rotation =
                (double) (-rad_to_deg(std::atan2((double) tmpy, signx))); // negative for image rotation

            // tmpx = g_x_min - g_x_3rd;
            sub_bf(bftmpx, g_bf_x_min, g_bf_x_3rd);
            LDBL tmpx2 = bftofloat(bftmpx);
            // tmpy = g_y_max - g_y_3rd;
            sub_bf(bftmpy, g_bf_y_max, g_bf_y_3rd);
            LDBL tmpy2 = bftofloat(bftmpy);
            LDBL b2 = tmpx2 * tmpx2 + tmpy2 * tmpy2;
            LDBL b = sqrtl(b2);

            double tmpa = std::acos((double) ((a2 + b2 - c2) / (2 * a * b))); // save tmpa for later use
            *Skew = 90 - rad_to_deg(tmpa);

            // these are the only two variables that must use big precision
            // *Xctr = (g_x_min + g_x_max)/2;
            add_bf(Xctr, g_bf_x_min, g_bf_x_max);
            half_a_bf(Xctr);
            // *Yctr = (g_y_min + g_y_max)/2;
            add_bf(Yctr, g_bf_y_min, g_bf_y_max);
            half_a_bf(Yctr);

            Height = b * std::sin(tmpa);
            *Magnification = 2 / Height; // 1/(h/2)
            *Xmagfactor = (double) (Height / (DEFAULT_ASPECT * a));

            // if vector_a cross vector_b is negative
            // then adjust for left-hand coordinate system
            if (tmpx1 * tmpy2 - tmpx2 * tmpy1 < 0 &&
                g_debug_flag != debug_flags::allow_negative_cross_product)
            {
                *Skew = -*Skew;
                *Xmagfactor = -*Xmagfactor;
                *Magnification = -*Magnification;
            }
        }
        if (*Magnification < 0)
        {
            *Magnification = -*Magnification;
            *Rotation += 180;
        }
*/
        restore_stack(saved);
    }


