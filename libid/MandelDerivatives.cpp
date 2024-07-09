/*
    MandelDerivatives.CPP a module for the per pixel calculations of fractals. 
    
    Written in Microsoft Visual 'C++' by Paul de Leeuw.

    This program is written in "standard" C. Hardware dependant code
    (console drivers & serial I/O) is in separate machine libraries.
*/

#include	<math.h>
#include	"Complex.h"
#include	"fractype.h"
#include	"resource.h"
#include	"cmplx.h"
#include    "bailout_formula.h"
#include    "pixel_grid.h"

extern  double g_params[];
extern  bailouts    g_bail_out_test;
extern  double g_magnitude_limit;
extern  DComplex    g_old_z, g_new_z;

static  Complex Sqr, sqrsqr, z, q;
static  double real_imag, RealImagSqr;
static  int degree, subtype;

bool BailoutTest(Complex *z, Complex SqrZ);

/**************************************************************************
	Initialise functions for each pixel
**************************************************************************/

int	   InitManDerFunctions(int subtype, Complex *z, Complex *q)
    {
    switch (subtype)
	    {
	    case 0:				// Perpendicular Mandelbrot
	    case 1:				// Burning Ship
	    case 3:				// Perpendicular Burning Ship
	    case 5:				// Perpendicular Buffalo
	    case 8:				// Mandelbar Celtic
	    case 9:				// Perpendicular Celtic
	    case 10:			// Cubic Flying Squirrel (Buffalo Imaginary)
	    case 11:			// Heart Mandelbrot
	    case 12:			// Celtic Heart
	    case 13:			// Partial Cubic Burning Ship Real
	    case 14:			// Partial Cubic Burning Ship Imaginary
	    case 15:			// Partial Cubic Buffalo Real (Celtic)
	    case 16:			// Cubic Quasi Burning Ship (Buffalo Burning Ship Hybrid)
	    case 17:			// Cubic Quasi Perpendicular
	    case 18:			// Cubic Quasi Heart
	    case 19:			// Mandelbrot 4th Order
	    case 20:			// Mandelbar 4th Order
	    case 21:			// Burning Ship 4th Order
	    case 22:			// Burning Ship 4th Partial Imag
	    case 23:			// Burning Ship 4th Partial Real
	    case 24:			// Burning Ship 4th Partial Real Mbar
	    case 25:			// Celtic Burning Ship 4th
	    case 26:			// Celtic Burning Ship 4th Partial Imag
	    case 27:			// Celtic Burning Ship 4th Partial Real
	    case 28:			// Celtic Burning Ship 4th Partial Real Mbar
	    case 29:			// Buffalo 4th Order
	    case 30:			// Buffalo 4th Partial Imag
	    case 31:			// Celtic (Buffalo 4th Partial Real)
	    case 32:			// Celtic 4th Mbar
	    case 33:			// False Quasi Perpendicular 4th
	    case 34:			// False Quasi Heart 4th
	    case 35:			// Celtic False Quasi Perpendicular 4th
	    case 36:			// Celtic False Quasi Heart 4th
	    case 37:			// Imag Quasi Perpendicular / Heart 4th
	    case 38:			// Real Quasi Perpendicular 4th
	    case 39:			// Real Quasi Heart 4th
	    case 40:			// Celtic Imag Quasi Perpendicular / Heart 4th
	    case 41:			// Celtic Real Quasi Perpendicular 4th
	    case 42:			// Celtic Real Quasi Heart 4th
	    case 43:			// Mandelbrot 5th
	    case 44:			// Mandelbar 5th (Vertical)
	    case 45:			// Mandelbar 5th (horizontal)
	    case 46:			// Burning Ship 5th
	    case 47:			// Buffalo 5th
	    case 48:			// Burning Ship 5th Partial
	    case 49:			// Burning Ship 5th Partial Mbar
	    case 50:			// Celtic 5th (Buffalo 5th Partial)
	    case 51:			// Celtic 5th Mbar
	    case 52:			// Quazi Burning Ship 5th (BS/Buffalo Hybrid)
	    case 53:			// Quazi Perpendicular 5th
	    case 54:			// Quazi Heart 5th
	    case 57:			// Kung Fu Panda
	    case 58:			// HPDZ Buffalo
	    case 59:			// SzegediButterfly 1
	    case 60:			// SzegediButterfly 2
//	        if (!juliaflag)                             // let's worry about Julai sets later
		        {
		        z->x = q->x + g_params[2];
		        z->y = q->y + g_params[3];
		        }
	        Sqr = 0;
	        real_imag = 0.0;

	        break;
	    case 2:				// Burning Ship of Higher Degree
	    case 4:				// Buffalo
	    case 6:				// Mandelbar (Tricorn)
	    case 7:				// Celtic
	    case 55:			// SimonBrot
	    case 56:			// SimonBrot2
//	        if (!juliaflag)
		        {
		        z->x = q->x + g_params[2];
		        z->y = q->y + g_params[3];
		        }
	        break;
	    }
    return 0;
    }

/**************************************************************************
	Run functions for each iteration
**************************************************************************/

int	   RunManDerFunctions(int subtype, Complex *z, Complex *q/*, BYTE *SpecialFlag, long *iteration*/)
    {
    switch (subtype)
	    {
	    case 0:						// Perpendicular Mandelbrot

/**************************************************************************
    Perpendicular Mandelbrot
    The equation is very similar to burning ship except you absolute just 
    the real factor when calculating the imaginary side.
    zi = abs(zr) * zi * -2.0 + JuliaI
    zr = zrsqr - zisqr + JuliaR

    Whereas burning ship is ...
    zi = abs(zr * zi) * 2.0 - JuliaI
    zr = zrsqr - zisqr - JuliaR
    Real -0.737,424,082,465,620,824,452,347,915,7­36,817,521,648,421,984,117,126,135,371,4
    Imaginary -0.355,631,230,436,806,548,631,431,830,9­06,449,574,310,522,006,013,120,497,532,0
    Zooms 182 magnification 6.13e54
***************************************************************************/

	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        real_imag = fabs(z->x) * z->y;
	        z->x = Sqr.x - Sqr.y + q->x;
	        z->y = -real_imag - real_imag + q->y;
	        return BailoutTest(z, Sqr);

	    case 1:						// Burning Ship
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        real_imag = fabs(z->x * z->y);
	        z->x = Sqr.x - Sqr.y + q->x;
	        z->y = real_imag + real_imag - q->y;
	        return BailoutTest(z, Sqr);

	    case 2:						// Burning Ship of Higher Degree
	        Sqr.x = z->x * z->x;    // only required to make bailout function work
	        Sqr.y = z->y * z->y;
	        z->x = fabs(z->x);
	        z->y = -fabs(z->y);
	        *z = z->CPolynomial(degree);
	        *z = *z + *q;
	        return BailoutTest(z, Sqr);
//	        return FractintBailoutTest(z);

	    case 3:						// Perpendicular Burning Ship

/**************************************************************************
    Perpendicular BurningShip
    The equation is very similar to burning ship except you absolute just 
    the imaginary factor instead of the real factor when calculating the imaginary side.
    zi = zr * abs(zi) * -2.0 + JuliaI
    zr = zrsqr - zisqr + JuliaR

    Whereas burning ship is ...
    zi = abs(zr * zi) * 2.0 - JuliaI
    zr = zrsqr - zisqr - JuliaR
    Real -0.737,424,082,465,620,824,452,347,915,7­36,817,521,648,421,984,117,126,135,371,4
    Imaginary -0.355,631,230,436,806,548,631,431,830,9­06,449,574,310,522,006,013,120,497,532,0
    Zooms 182 magnification 6.13e54
***************************************************************************/
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        real_imag = z->x * fabs(z->y);
	        z->x = Sqr.x - Sqr.y + q->x;
	        z->y = real_imag + real_imag + q->y;
	        return BailoutTest(z, Sqr);

	    case 4:						// Buffalo (works according to Kalles Fraktaller)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        if (degree == 2)
		        {
		        z->y = fabs(z->x) * fabs(z->y) * -2.0 - q->y;
		        z->x = fabs(Sqr.x - Sqr.y) + q->x;
		        }
    	    else if (degree == 3)
		        {
		        z->y = fabs(((Sqr.x * 3.0) - Sqr.y) * z->y) - q->y;
		        z->x = fabs((Sqr.x - (Sqr.y * 3.0)) * z->x) + q->x;
		        }
	        else	// degree > 3
		        {
		        *z = z->CPolynomial(degree);
		        z->x = fabs(z->x) + q->x;
		        z->y = fabs(z->y) - q->y;
		        }
	        return BailoutTest(z, Sqr);

	    case 5:						// Perpendicular Buffalo - (according to Kalles Fraktaller)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = (z->x) * fabs(z->y) * -2.0 - q->y;
	        z->x = fabs(Sqr.x - Sqr.y) + q->x;
	        return BailoutTest(z, Sqr);

	    case 6:						// Mandelbar (Tricorn)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        if (degree == 2)
		        {
		        real_imag = z->x * z->y;
		        z->x = Sqr.x - Sqr.y + q->x;
		        z->y = -real_imag - real_imag + q->y;
		        }
	        else
		        {
		        *z = z->CPolynomial(degree);
		        z->x = (g_params[3] == 1.0 ? -z->x : z->x) + q->x;
		        z->y = (g_params[3] == 1.0 ? z->y : -z->y) + q->y;
		        }
	        return BailoutTest(z, Sqr);

	    case 7:						// Celtic
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        if (degree == 2)
		        {
		        z->y = (z->x) * (z->y) * -2.0 - q->y;
		        z->x = -fabs(Sqr.x - Sqr.y) - q->x;
		        }
	        else	// degree > 2
		        {
		        *z = z->CPolynomial(degree);
		        z->x = fabs(z->x);
		        *z = *z + *q;
		        }
	        return BailoutTest(z, Sqr);

	    case 8:						// Mandelbar Celtic
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = z->x * z->y * -2.0 + q->y;
	        z->x = fabs(Sqr.x - Sqr.y) + q->x;
	        return BailoutTest(z, Sqr);

	    case 9:						// Perpendicular Celtic
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = fabs(z->x) * z->y * -2.0 + q->y;
	        z->x = fabs(Sqr.x - Sqr.y) + q->x;
	        return BailoutTest(z, Sqr);

	    case 10:					// Cubic Flying Squirrel (Buffalo Imaginary)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = fabs(((Sqr.x * 3.0) - Sqr.y) * z->y) - q->y;
	        z->x = ((Sqr.x - (Sqr.y * 3.0)) * z->x) + q->x;
	        return BailoutTest(z, Sqr);

	    case 11:					// Heart Mandelbrot
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = fabs(z->x) * z->y * 2.0 - q->y;
	        z->x = Sqr.x - Sqr.y + q->x;
	        return BailoutTest(z, Sqr);

	    case 12:					// Celtic Heart
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = fabs(z->x) * z->y * 2.0 - q->y;
	        z->x = fabs(Sqr.x - Sqr.y) + q->x;
	        return BailoutTest(z, Sqr);

	    case 13:					// Partial Cubic Burning Ship Real
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = ((Sqr.x * 3.0) - Sqr.y) * z->y + q->y;
	        z->x = (Sqr.x - (Sqr.y * 3.0)) * fabs(z->x) + q->x;
	        return BailoutTest(z, Sqr);

	    case 14:					// Partial Cubic Burning Ship Imaginary
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = ((Sqr.x * 3.0) - Sqr.y) * fabs(z->y) + q->y;
	        z->x = (Sqr.x - (Sqr.y * 3.0)) * z->x + q->x;
	        return BailoutTest(z, Sqr);

	    case 15:					// Partial Cubic Buffalo Real (Celtic)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = ((Sqr.x * 3.0) - Sqr.y) * z->y + q->y;
	        z->x = fabs((Sqr.x - (Sqr.y * 3.0)) * z->x) + q->x;
	        return BailoutTest(z, Sqr);

	    case 16:					// Cubic Quasi Burning Ship (Buffalo Burning Ship Hybrid)
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = -fabs(((Sqr.x * 3.0) - Sqr.y) * z->y) + q->y;
	        z->x = (Sqr.x - (Sqr.y * 3.0)) * fabs(z->x) + q->x;
	        return BailoutTest(z, Sqr);

	    case 17:					// Cubic Quasi Perpendicular
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = -fabs((Sqr.x * 3.0) - Sqr.y) * z->y + q->y;
	        z->x = (Sqr.x - (Sqr.y * 3.0)) * fabs(z->x) + q->x;
	        return BailoutTest(z, Sqr);

	    case 18:					// Cubic Quasi Heart
	        Sqr.x = z->x * z->x;
	        Sqr.y = z->y * z->y;
	        z->y = fabs((Sqr.x * 3.0) - Sqr.y) * z->y + q->y;
	        z->x = (Sqr.x - (Sqr.y * 3.0)) * fabs(z->x) + q->x;
	        return BailoutTest(z, Sqr);
//	    return ((Sqr.x + Sqr.y) >= rqlim);

    /****************************************************************
	4th Order Fractals:
    ****************************************************************/

    /****************************************************************
	Non ABS Variations (2)
    ****************************************************************/

	case 19:					// Mandelbrot 4th Order
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 20:					// Mandelbar 4th Order
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * z->x * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

    /****************************************************************
	***Straight ABS Variations (16)
    ****************************************************************/

	case 21:					// Burning Ship 4th Order
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x * z->y) * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 22:					// Burning Ship 4th Partial Imag
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * fabs(z->y) * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 23:					// Burning Ship 4th Partial Real
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x) * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 24:					// Burning Ship 4th Partial Real Mbar
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * fabs(z->x) * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 25:					// Celtic Burning Ship 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x * z->y) * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 26:					// Celtic Burning Ship 4th Partial Imag
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * fabs(z->y) * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 27:					// Celtic Burning Ship 4th Partial Real
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x) * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 28:					// Celtic Burning Ship 4th Partial Real Mbar
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * fabs(z->x) * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 29:					// Buffalo 4th Order
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x * z->y * (Sqr.x - Sqr.y)) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 30:					// Buffalo 4th Partial Imag
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x * z->y * (Sqr.x - Sqr.y)) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 31:					// Celtic (Buffalo 4th Partial Real)
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 32:					// Celtic 4th Mbar
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * z->x * z->y * (Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

    /**************************************************************************
	Quasi ABS Variations (10)
    ***************************************************************************/

	case 33:					// False Quasi Perpendicular 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * z->x * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 34:					// False Quasi Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 35:					// Celtic False Quasi Perpendicular 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * z->x * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 36:					// Celtic False Quasi Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 37:					// Imag Quasi Perpendicular / Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * fabs(z->y) * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 38:					// Real Quasi Perpendicular 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * fabs(z->x) * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);
//	    return ((Sqr.x + Sqr.y) >= rqlim);

	case 39:					// Real Quasi Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x) * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y + q->x;
	    return BailoutTest(z, Sqr);

	case 40:					// Celtic Imag Quasi Perpendicular / Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * z->x * fabs(z->y) * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 41:					// Celtic Real Quasi Perpendicular 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = -4.0 * fabs(z->x) * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 42:					// Celtic Real Quasi Heart 4th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    z->y = 4.0 * fabs(z->x) * z->y * fabs(Sqr.x - Sqr.y) + q->y;
	    z->x = fabs(Sqr.x * Sqr.x + Sqr.y * Sqr.y - 6.0 * Sqr.x * Sqr.y) + q->x;
	    return BailoutTest(z, Sqr);

    /****************************************************************
	5th Order Fractals:
    ****************************************************************/

	case 43:					// Mandelbrot 5th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 44:					// Mandelbar 5th (Vertical)
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = -z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 45:					// Mandelbar 5th (horizontal)
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = -z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 46:					// Burning Ship 5th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = fabs(z->y) * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 47:					// Buffalo 5th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = fabs(z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y)) + q->y;
	    z->x = fabs(z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y)) + q->x;
	    return BailoutTest(z, Sqr);

	case 48:					// Burning Ship 5th Partial
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 49:					// Burning Ship 5th Partial Mbar
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = -z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 50:					// Celtic 5th (Buffalo 5th Partial)
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y)) + q->x;
	    return BailoutTest(z, Sqr);

	case 51:					// Celtic 5th Mbar
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = -z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y)) + q->x;
	    return BailoutTest(z, Sqr);

	case 52:					// Quazi Burning Ship 5th (BS/Buffalo Hybrid)
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = -fabs(z->y * (5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y)) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 53:					// Quazi Perpendicular 5th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = -z->y * fabs(5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 54:					// Quazi Heart 5th
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;
	    sqrsqr.x = Sqr.x * Sqr.x;
	    sqrsqr.y = Sqr.y * Sqr.y;
	    RealImagSqr = Sqr.x * Sqr.y;
	    z->y = z->y * fabs(5.0 * sqrsqr.x - 10.0 * RealImagSqr + sqrsqr.y) + q->y;
	    z->x = fabs(z->x) * (sqrsqr.x - 10.0 * RealImagSqr + 5.0 * sqrsqr.y) + q->x;
	    return BailoutTest(z, Sqr);

	case 55:					// SimonBrot
/**************************************************************************
	Run SimonBrot type fractals
	z^n * |z|^2 + c	    (normal)

	It would also be nice to see Simonbrot 3rd, 5th, and 7th, but the formulas for these fractals are
	z1.5 * |z|2 + pixel, z2.5 * |z|2 + pixel, and z3.5 * |z|2 + pixel. Due to the fact that the powers of
	some of the terms are fractions, this could be difficult to implement. In a roundabout way, there is
	a power zero Simonbrot on Kalles Fraktaler, because the power zero Simonbrot is actually the Burning Ship.
**************************************************************************/
	    {
	    Complex	zabs, tempz, sqrtz;

	    zabs.x = fabs(z->x);
	    zabs.y = fabs(z->y);
	    tempz.y = z->y * zabs.x + z->x * zabs.y;
	    tempz.x = z->x * zabs.x - z->y * zabs.y;
	    sqrtz = (degree % 2 == 1) ? z->CSqrt() : 1.0;		// use square root power if degree is odd
	    *z = tempz.CPolynomial(degree / 2) * sqrtz + *q;
	    return BailoutTest(z, Sqr);
	    }

	case 56:					// SimonBrot2
/**************************************************************************
	Run SimonBrot2 type fractals
	z^n * |z^2| + c	    (SimonBrot2)
**************************************************************************/
	    {
	    Complex	zabs, tempz, sqrtz;

	    tempz = *z * *z;
	    zabs.x = fabs(tempz.x);
	    zabs.y = -fabs(tempz.y);
	    tempz = zabs;
	    sqrtz = (degree % 2 == 1) ? z->CSqrt() : 1.0;		// use square root power if degree is odd
	    *z = z->CPolynomial(degree / 2) * sqrtz * tempz + *q;
	    return BailoutTest(z, Sqr);
	    }

	case 57:					// Kung Fu Panda
/**************************************************************************
	Kung Fu Panda type fractals
	z = abs(z*z)
	z = z * z + p
**************************************************************************/
	    {
	    Complex	t1, t2;

	    t1 = *z * *z;
	    t2.x = fabs(t1.x);
	    t2.y = fabs(t1.y);
	    *z = t2 * t2 - *q;
	    return BailoutTest(z, Sqr);
	    }

	case 58:					// HPDZ Buffalo
/**************************************************************************
	HPDZ Buffalo type fractals
	z := (((x^2 - y^2) - |x|) + i (|2xy| - |y|)) + c
	or
	w := |x| + i |y|
	z := w^2 - w + c
**************************************************************************/
	    {
	    Complex w;

	    w.x = fabs(z->x);
	    w.y = fabs(z->y);
	    *z = w * w - w;
	    z->x += q->x;
	    z->y -= q->y;
	    return BailoutTest(z, Sqr);
	    }

	case 59:					// SzegediButterfly 1
/**************************************************************************
	Szegedi Butterfly
	double temp = complex[0].getIm();
	double temp2 = complex[0].getRe();
	complex[0] = new Complex(temp * temp - Math.sqrt(complex[0].getAbsRe()), temp2 * temp2 - Math.sqrt(complex[0].getAbsIm())).plus_mutable(complex[1]);
	complex[0] = new Complex(temp2 * temp2 - Math.sqrt(complex[0].getAbsIm()), temp * temp - Math.sqrt(complex[0].getAbsRe())).plus_mutable(complex[1]);
**************************************************************************/
	    {
	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;

	    z->x = Sqr.y - sqrt(fabs(z->x));
	    z->y = Sqr.x - sqrt(fabs(z->y));
	    *z += *q;

	    return BailoutTest(z, Sqr);
	    }

	case 60:					// SzegediButterfly 2
	    {
	    double  temp;

	    Sqr.x = z->x * z->x;
	    Sqr.y = z->y * z->y;

	    temp = Sqr.x - sqrt(fabs(z->y));
	    z->y = Sqr.y - sqrt(fabs(z->x));
	    z->x = temp;
	    *z += *q;

	    return BailoutTest(z, Sqr);
	    }


	}
    return 0;
    }

/**************************************************************************
    Bailout Test
**************************************************************************/

bool BailoutTest(Complex *z, Complex SqrZ)
    {
    //    Complex TempSqr;
    double magnitude;
    double manhmag;
    double manrmag;

    switch (g_bail_out_test)
        {
        case bailouts::Mod:
            magnitude = SqrZ.x + SqrZ.y;
            return (magnitude >= g_magnitude_limit);

        case bailouts::Real:
            return (SqrZ.x >= g_magnitude_limit);

        case bailouts::Imag:
            return (SqrZ.y >= g_magnitude_limit);

        case bailouts::Or:
            return (SqrZ.x >= g_magnitude_limit || SqrZ.y >= g_magnitude_limit);

        case bailouts::And:
            return (SqrZ.x >= g_magnitude_limit && SqrZ.y >= g_magnitude_limit);

        case bailouts::Manh:
            manhmag = fabs(z->x) + fabs(z->y);
            return ((manhmag * manhmag) >= g_magnitude_limit);

        case bailouts::Manr:
            manrmag = z->x + z->y; // don't need abs() since we square it next
            return ((manrmag * manrmag) >= g_magnitude_limit);

        default:
            magnitude = SqrZ.x + SqrZ.y;
            return (magnitude >= g_magnitude_limit);
        }
    }

/**************************************************************************
	Initialise functions for each pixel
**************************************************************************/

int init_mand_derivatives()
    {
    z.x = g_old_z.x;
    z.y = g_old_z.y;
    q.x = g_dx_pixel();
    q.y = g_dy_pixel();
    subtype = (int) g_params[0];
    degree = (int) g_params[1];
    g_bail_out_test = (bailouts) g_params[4];
    if (degree < 1)
        degree = 1;
    InitManDerFunctions(subtype, &z, &q);
    return 0;
    }

/**************************************************************************
	Run functions for each pixel
**************************************************************************/

int run_mand_derivatives()
    {
    int ReturnMode;
    ReturnMode = RunManDerFunctions(subtype, &z, &q);
    g_new_z.x = z.x;
    g_new_z.y = z.y;
    return ReturnMode;
    }

