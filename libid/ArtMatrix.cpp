/*
    ARTMATRIX.CPP a module for the per pixel calculations of Art Matrix fractals. 
    
    Written in Microsoft Visual 'C++' by Paul de Leeuw.

    This program is written in "standard" C. Hardware dependant code
    (console drivers & serial I/O) is in separate machine libraries.
*/

#include	<math.h>
#include	"fractype.h"
#include	"Complex.h"

#include	"cmplx.h"
#include    "bailout_formula.h"
#include    "pixel_grid.h"
#include    "id.h"

#define     ERROR   0

extern double g_params[];
extern bailouts g_bail_out_test;
extern double g_magnitude_limit;
extern DComplex g_old_z, g_new_z;
extern long g_max_iterations; // try this many iterations
extern long g_color_iter;
extern int g_color;
extern  bool juliaflag; // for the time being

static int type, subtype, special;
static Complex a, a2, aa3, alpha, b, cmcc, l2, lm5, lp5, oz, q, t2, t3, temp, temp1, temp3, temp4, v, z;
static double absolute, distance, der, epsilon, escape;
static int OldMaxIterations = 0; // check to see if max iterations changes

static int penpref[4] = {1, 2, 3, 4}; // can we do more in mappying these colours?
static int pennref[4] = {9, 10, 11, 12};
static int penp[4];
static int penn[4];

static long color;

unsigned char g_phaseflag;

/*
param 0 = type
param 1 = subtype
param 2 = special colour
param 3 = x pert
param 4 = y pert
*/


/**************************************************************************
	Initialise functions for each pixel
**************************************************************************/

int	InitArtMatrix(int type, Complex *z, Complex *q)
    {
    switch (type)
	    {
	    case 0:					        // Art Matrix Cubic
            {
	        switch ((int)g_params[1])
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

//	        period_level = FALSE;			            // no periodicity checking
	        if (subtype == 'B')				            // CBIN 
		        {
		        t3 = *q * 3;				            // T3 = 3*T
		        t2 = q->CSqr();				            // T2 = T*T
		        a = (t2 + 1) / t3;			            // A  = (T2 + 1)/T3
							                            // B  = 2*A*A*A + (T2 - 2)/T3    
		        temp = a.CCube() * 2;			        // 2*A*A*A
		        b = (t2 - 2) / t3 + temp;		        // B  = 2*A*A*A + (T2 - 2)/T3
		        }
	        else if (subtype == 'C' || subtype == 'F')	// CCIN or CFIN
		        {
		        a = *q;					                // A = T
							                            // find B = T + 2*T*T*T
		        temp = q->CCube();			            // B = T*T*T
		        if (subtype == 'C')
		            b = temp + temp + *q;		        // B = B * 2 + T
		        else
		            {
		            b = (temp - *q) * 2;		        // B = B * 2 - 2 * T
		            a2 = a + a;
		            }
		        }
	        else if (subtype == 'K')			        // CKIN 
		        {
		        a = 0;
		        v = 0;
		        b = *q;					                // B = T
		        }
	        aa3 = a.CSqr() * 3;				            // AA3 = A*A*3
	        if (!juliaflag)
		        *z = -a;				                // Z = -A
	        break;
            }
	    case 1:				                            // Art Matrix Newton
	        l2 = q->CSqr();				                // L2 = L*L
	        a = -l2 + 0.25;				                // A = ( .25,0) - L2
	        b = -l2 - 0.75;				                // B = (-.75,0) - L2 
	        lm5 = *q - 0.5;				                // LM5 = L - (.5,0)
	        lp5 = *q + 0.5;				                // LP5 = L + (.5,0)
	        break;

	    case 2:					                        // Art Matriuc Matein fractal
	        if ((absolute = q->CSumSqr()) > 1.0)
		        return(-1);				                // not inside set
	        if (!juliaflag)
		        *z = 1;

	        for (int i = 0; i < 100; ++i)		        // DO 300 I = 1,100 
		        {
		        temp = z->CInvert();			        // 300  Z = L*(Z + 1/Z)
		        *z = *q * (*z + temp);
		        }

	        distance = 1.0;				                // D = 1
	        oz = z->CInvert();				            // OZ = 1/Z
	        break;


	    case 3:				                            // Art Matrix Rational Map
            {
	        switch ((int)g_params[1])
		        {
		        case 0:
		            subtype = 'A';
		            break;
		        case 1:
		            subtype = 'B';
		            break;
		        default:
		            subtype = 'A';
		            break;
		        }
	        special = (int)g_params[2];
	        if (special < 0)
		        special = 0;
	        if (g_max_iterations != OldMaxIterations)
		        {
		        OldMaxIterations = g_max_iterations;
		        int gap = g_max_iterations / 16;			// split the colour map into 16 equal parts
		        for (int i = 0; i < 4; i++)
		            {
		            penp[i] = penpref[i] * gap;
		            penn[i] = pennref[i] * gap;
		            }
		        }

	        if (subtype == 'A')
		        {
		        cmcc = *q - q->CSqr();			        // CMCC = C - C*C
		        temp = -*q + 2;
		        a = temp / cmcc;			            // A = (2 - C)/CMCC
		        temp = cmcc + 1;
		        b = -temp / cmcc;			            // B = -(CMCC + 1)/CMCC

							                            // ALPHA = 1/(C*C * (B + B + B*B/A) * (2*A*C + B))
		        temp = a * *q * 2 + b;			        // 2*A*C + B
		        temp1 = b.CSqr() / a + b + b;		    // (B + B + B*B/A)
		        temp3 = q->CSqr()*temp1*temp;
		        alpha = temp3.CInvert();
		        }
	        else if (subtype == 'B')
		        {
		        a = *q;
		        b = a + 1;
		        temp = a.CSqr() - 1;
		        alpha = a / temp;			            // ALPHA = A/(A*A - 1)
		        }
	        else
		        return(ERROR);				            // unknown subtype

							                            // ESCAPE  =   4/ABS(ALPHA)
							                            // ESCAPE  =   ESCAPE*ESCAPE
	        if (alpha.x != 0.0 || alpha.y != 0.0)
		        escape = 16.0 / alpha.CSumSqr();
	        else
		        return(FALSE);				            // no naughty division
	        epsilon = 0.000001 / escape;		        // EPSILN = 0.000001/ESCAPE 

	        der = 1.0;					                // DER = 1.0 
	        if (!juliaflag)
		        {					                    // Z = -B/(A + A)
		        temp = -a * 2;
		        *z = b / temp;
		        }
	    // iterating Z  = 1/(A*Z*Z + B*Z + 1) has various proterties:
	    // 	   Z  = 1/(A*Z*Z + B*Z + 1)
	    // Julia   Z  = 1/(A*Z*Z + B*Z + 1)
	    // Julia   Z  = 1/(A*Z*Z + B*Z + 1)
	    // ?????   Z  = 1/(A*Z*Z + B*Z + 1)

	        int	zcount;
	        switch (subtype)
		        {
		        case 'A':
		            if (juliaflag)
			            zcount = 4;
		            else
			            zcount = 2;
		            break;
		        case 'B':
		            zcount = 3;
		            break;
		        }

	        for (int i = 0; i < zcount; ++i)
		        {
							                        // 1/(A*Z*Z + B*Z + 1)
		        temp = b * *z + 1;			        // B*Z + 1
		        temp1 = z->CSqr()*a + temp;		    // (A*Z*Z + B*Z + 1)
		        *z = temp1.CInvert();			    // Z = 1/(A*Z*Z + B*Z + 1)
		        }
	        break;
            }
        }
    return 0;
    }

/**************************************************************************
	Run functions for each iteration
**************************************************************************/

int RunArtMatrix(int type, Complex *z, Complex *q)
    {
    switch (type)
	    {
	    case 0:					            // Art Matrix Cubic
	        if (subtype == 'K')				// CKIN
		        {
		        *z = z->CCube() + b;			// Z = Z*Z*Z + B
		        z->x += g_params[3];
		        z->y += g_params[4];
		        }
	        else
		        {
		        temp = z->CCube() + b;			// Z = Z*Z*Z + B
		        *z = temp - aa3 * *z;			// Z = Z*Z*Z - AA3*Z + B
		        z->x += g_params[3];
		        z->y += g_params[4];
		        }
	        if (z->CSumSqr() > 100.0)
		        return (TRUE);
	        else
		        {
		        if (subtype == 'F')
		            {
		            if (q->CSumSqr() < 0.111111)
			            {
			            g_color_iter = special;
//        			    *SpecialFlag = TRUE;			// for decomp and biomorph
			            return (TRUE);
			            }
		            v = *z + a2;
		            }
		        else if (subtype == 'K')
		            v = *z - v;
		        else
		            v = *z - a;
		        if (v.CSumSqr() <= 0.000001)
		            {
		            g_color_iter = special;
    //		        *SpecialFlag = TRUE;			// for decomp and biomorph
		            return (TRUE);
		            }
		        return (FALSE);
		        }

	    case 1:				                        // Art Matrix Newton
	        {
		    special = (int)g_params[2];
	        if (special < 0)
		        special = 2;
	        Complex z2 = z->CSqr();			        // z2 = z*z
							                        // Z  =  (2*Z*Z2 + A)/(3*Z2 + B)
	        Complex top = z2 * *z * 2 + a;
	        Complex bottom = z2 * 3 + b;
	        *z = top / bottom;
	        z->x += g_params[3];
	        z->y += g_params[4];
	        v = *z - 1;
	        //    v.x += param[1];
	        //    v.y += param[2];
	        if (v.CSumSqr() <= 0.000001)
		        {
		        g_phaseflag = 0;				        // first phase
		        return(TRUE);
		        }
	        v = *z - lm5;
	        if (v.CSumSqr() <= 0.000001)
		        {
		        g_phaseflag = 1;				        // second phase
                g_color += special;
		        return(TRUE);
		        }
	        v = *z + lp5;
	        if (v.CSumSqr() <= 0.000001)
		        {
		        g_phaseflag = 2;				        // third phase
                g_color += special * 2;
		        return(TRUE);
		        }
	        return(FALSE);
	        }

	    case 2:					// Art Matriuc Matein fractal
	        epsilon = 0.01;
	        escape = 10.0E20;
	        *z = *q * (*z + oz);			// Z = L*(Z + OZ)
	        z->x += g_params[3];
	        z->y += g_params[4];
	        oz = z->CInvert();				// OZ = 1/Z
	        temp = -oz / *z + 1;			// T = 1 - OZ/Z
							                // D = D*ABSL*(REAL(T)*REAL(T) + IMAG(T)*IMAG(T))
	        distance = distance * absolute * temp.CSumSqr();

	        if (distance <= epsilon)
		        {
		        g_phaseflag = 0;			        // first phase
		        return(TRUE);
		        }
	        if (distance > escape)
		        {
		        g_phaseflag = 1;			        // second phase
		        return(TRUE);
		        }
	        return(FALSE);

	    case 3:				                // Art Matrix Rational Map 
	        {
	        Complex az = a * *z;			// AZ = A*Z
	        temp = az * *z + 1;
	        temp1 = b * *z + temp;
	        *z = temp1.CInvert();			// Z = 1/(AZ*Z + B*Z + 1)
	        z->x += g_params[3];
	        z->y += g_params[4];
	        if (z->CSumSqr() > escape)
		    {
		    if ((alpha.x * z->y + alpha.y * z->x) <= 0.0)
		        color = penp[g_color_iter % 4];
		    else
		        color = penn[g_color_iter % 4];
		    g_color_iter = color;
		    return (TRUE);
		    }
	        temp = az * 2 + b;
	        Complex d = z->CSqr()*temp;		// D = (2*AZ + B)*Z*Z
	        double dist = d.CSumSqr();
	        der *= dist;				    // DER  =   DER*(REAL(D)*REAL(D) + IMAG(D)*IMAG(D))
	        if (der < epsilon)
		        {
		        g_color_iter = g_max_iterations;
		        return(TRUE);
		        }
	        if (g_color_iter >= g_max_iterations)
		        {
		        g_color_iter = special;
		        return(TRUE);
		        }
            }
	    return(FALSE);				        // continue iterations
	    }

    return 0;
    }

/**************************************************************************
	Initialise functions for each pixel
**************************************************************************/

int init_art_matrix()
    {
    z.x = g_old_z.x;
    z.y = g_old_z.y;
    q.x = g_dx_pixel();
    q.y = g_dy_pixel();
    type = (int) g_params[0];
    subtype = (int) g_params[1];
//    g_bail_out_test = (bailouts) g_params[4];
    InitArtMatrix(type, &z, &q);
    return 0;
    }

/**************************************************************************
	Run functions for each orbit
**************************************************************************/

int run_art_matrix()
    {
    int ReturnMode;
    ReturnMode = RunArtMatrix(type, &z, &q);
    g_new_z.x = z.x;
    g_new_z.y = z.y;
    return ReturnMode;
    }

