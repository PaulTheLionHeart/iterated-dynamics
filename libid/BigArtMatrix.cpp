/*
    BIGFUNCTIONS.CPP a module for the per pixel calculations of Bignum fractals.
    
    Written in Microsoft Visual 'C++' by Paul de Leeuw.

    This program is written in "standard" C. Hardware dependant code
    (console drivers & serial I/O) is in separate machine libraries.
*/

#include	<math.h>
#include	"fractype.h"
#include	"Complex.h"
#include	"big.h"
#include	"BigDouble.h"
#include	"BigComplex.h"

extern double g_params[];
extern BFComplex bfparm, bfnew, bfold;
extern int g_row;
extern int g_col;
extern double g_magnitude_limit;
extern long g_max_iterations; // try this many iterations
extern long g_color_iter;

static BigComplex aBig, bBig, c1Big, c2Big, caa3Big, cbBig, zBig, z1Big, z2Big, z3Big, z4Big, zdBig, ztBig,
    qBig;
static BigDouble BigBailout, tBig, realimagBig, RealImagSqrBig;
static int type, subtype, special;
static double d;

// extern	void	ShowBignum(BigDouble x, char *Location);
extern void bf2BigNum(BigDouble *BigNum, bf_t bfNum);
extern void BigNum2bf(bf_t *bfNum, BigDouble BigNum);

extern bool juliaflag; // for the time being - declared in TierazonFunctions.cpp

extern void ShowBignum(BigDouble x, char *Location);


/**************************************************************************
	Initialise functions for each pixel
**************************************************************************/



#ifdef SPLAT

int	BigInitFunctions(WORD type, BigComplex *zBig, BigComplex *qBig)
    {
    switch (type)
	    {
	    case 0:					// Art Matrix Cubic
	        {
	        BigComplex	tempBig;

	        switch ((int) g_params[1])
		    {
		    case 0:
		        subtype = 'B';
		        break;
		    case 1:
		        subtype = 'C';
		        break;
		    case 2:
		        subtype = 'F';
		        break;
		    case 3:
		        subtype = 'K';
		        break;
		    default:
		        subtype = 'B';
		        break;
		    }
	        special = (int)g_params[2];
	        if (special < 0)
		        special = 0;

//	        period_level = FALSE;			// no periodicity checking
	        if (subtype == 'B')				// CBIN 
		    {
		    t3Big = *qBig * 3.0;			// T3 = 3*T
		    t2Big = qBig->CSqr();			// T2 = T*T
		    aBig = (t2Big + 1.0) / t3Big;		// A  = (T2 + 1)/T3
							    // B  = 2*A*A*A + (T2 - 2)/T3    
		    tempBig = aBig.CCube();			// A*A*A
		    tempBig = tempBig.CDouble();		// 2*A*A*A
		    bBig = (t2Big - 2.0) / t3Big + tempBig;	// B  = 2*A*A*A + (T2 - 2)/T3
		    }
	        else if (subtype == 'C' || subtype == 'F')	// CCIN or CFIN
		    {
		    aBig = *qBig;				// A = T
							    // find B = T + 2*T*T*T
		    tempBig = qBig->CCube();		// B = T*T*T
		    if (subtype == 'C')
		        bBig = tempBig.CDouble() + *qBig;	// B = B * 2 + T
		    else
		        {
		        bBig = (tempBig - *qBig);		// B = B - T
		        bBig = bBig.CDouble();		// B = B * 2 - 2 * T
		        a2Big = aBig.CDouble();
		        }
		    }
	        else if (subtype == 'K')			// CKIN 
		    {
		    aBig = 0;
		    vBig = 0;
		    bBig = *qBig;				// B = T
		    }
	        aa3Big = aBig.CSqr()*3.0;			// AA3 = A*A*3
	        if (!juliaflag)
		    *zBig = -aBig;				// Z = -A
	        }
	        break;

	    case 1:				// Art Matrix Newton
	        l2Big = qBig->CSqr();
	        aBig = -l2Big + 0.25;
	        bBig = -(l2Big + 0.75);
	        lm5Big = *qBig - 0.5;
	        lp5Big = *qBig + 0.5;
	        break;

	    case 2:					// Art Matriuc Matein fractal
	        {
//	        period_level = FALSE;			// no periodicity checking
	        double	one = 1.0;
	        if ((absolute = qBig->CSumSqr()) > one)
		    return(-1);				// not inside set
	        zBig->x = 1.0;
	        zBig->y = 0.0;
	        /* DO 300 I = 1,100 */
	        /* 300  Z = L*(Z + 1/Z) */
	        for (int i = 0; i < 100; ++i)
		    *zBig = (*zBig + zBig->CInvert()) * *qBig;

	        distance = 1.0;				// D = 1 
	        ozBig = zBig->CInvert();			// OZ = 1/Z
	        break;
	        }

	    }
    return 0;
    }


/**************************************************************************
	Run functions for each iteration
**************************************************************************/

int	BigRunFunctions(WORD type, BigComplex *zBig, BigComplex *qBig, BYTE *SpecialFlag, long *iteration)
    {
    switch (type)
	    {
	    case 0:					// Art Matrix Cubic
	        {
	        BigComplex	tempBig;

	        if (subtype == 'K')				// CKIN
		    *zBig = zBig->CCube() + bBig;		// Z = Z*Z*Z + B
	        else					// Z = Z*Z*Z - AA3*Z + B
		    {
		    tempBig = zBig->CCube() + bBig;		// Z = Z*Z*Z + B
		    *zBig = tempBig - aa3Big * *zBig;	// Z = Z*Z*Z - AA3*Z + B
		    }
	        if (zBig->CSumSqr() >= 100.0)
		    return (TRUE);
	        else
		    {
		    if (subtype == 'F')
		        {
		        if (qBig->CSumSqr() <= 0.111111)
			    {
			    *iteration = special;
			    *SpecialFlag = TRUE;		// for decomp and biomorph
			    return (TRUE);
			    }
		        vBig = *zBig + a2Big;
		        }
		    else if (subtype == 'K')
		        vBig = *zBig - vBig;
		    else
		        vBig = *zBig - aBig;
		    if (vBig.CSumSqr() <= 0.000001)
		        {
		        *iteration = special;
		        *SpecialFlag = TRUE;		// for decomp and biomorph
		        return (TRUE);
		        }
		    }
	        return(FALSE);
	        }

	    case 1:				// Art Matrix Newton
	        {
	        BigComplex	z2Big;

	        z2Big = zBig->CSqr();
	        temp1Big = z2Big * zBig->CDouble() + aBig;
	        temp2Big = z2Big * 3.0 + bBig;
	        *zBig = temp1Big / temp2Big;

	        vBig = *zBig - 1.0;
	        if (vBig.CSumSqr() <= 0.000001)
		    {
		    phaseflag = 0;				// first phase
		    return(TRUE);
		    }
	        // v_real = dz_real - lm5_real;
	        vBig = *zBig - lm5Big;			// v_imag = dz_imag - lm5_imag;
	        if (vBig.CSumSqr() <= 0.000001)
		    {
		    phaseflag = 1;				// second phase
		    return(TRUE);
		    }
	        // v_real = dz_real + lp5_real;
	        vBig = *zBig + lp5Big;			// v_imag = dz_imag + lp5_imag;
	        if (vBig.CSumSqr() <= 0.000001)
		    {
		    phaseflag = 2;				// third phase
		    return(TRUE);
		    }
	        return(FALSE);
	        }

	    case 2:					// Art Matriuc Matein fractal
	        {
	        double	epsilon = 0.01;
	        double	escape = 10.0E20;
	        BigComplex	t;

	        *zBig = *qBig * (*zBig + ozBig);		// Z = L*(Z + OZ)
	        ozBig = zBig->CInvert();			// OZ = 1/Z
	        t = -ozBig / *zBig;				// T = 1 - OZ/Z
	        t.x = t.x + 1.0;
	        // D = D*ABSL*(REAL(T)*REAL(T) + IMAG(T)*IMAG(T))
	        distance = distance * absolute * t.CSumSqr();

	        if (distance <= epsilon)
		    {
		    phaseflag = 0;				// first phase
		    return(TRUE);
		    }
	        if (distance > escape)
		    {
		    phaseflag = 1;				// second phase
		    return(TRUE);
		    }
	        return(FALSE);
	        }

	    }
    return 0;
    }
#endif

