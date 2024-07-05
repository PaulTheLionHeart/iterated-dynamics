/****************************************************
    FRAMECALCULATOR.CPP a module to explore Perturbation
    Thanks to Claude Heiland-Allen https://fractalforums.org/programming/11/perturbation-code-for-cubic-and-higher-order-polynomials/2783
    Written in Microsoft Visual C++ by Paul de Leeuw.
****************************************************/

#include <windows.h>
#include <time.h>
#include "Big.h"
//#include "Dib.h"
#include "filter.h"
#include "Big.h"
//#include "Manpwin.h"
//#include "mpfr.h"
#include "PertEngine.h"
#include "Complex.h"
#include "Drivers.h"
// #include "../../../../Program Files (x86)/Windows Kits/10/Include/10.0.19041.0/ucrt/stdio.h"





/*
extern	void	InitFilter(Complex q);
extern	int	DoTierazonFilter(Complex q);
extern	int	EndTierazonFilter(void);
extern	int	method;				// inside and outside filters
extern	long	iteration;			// globals for speed 
extern	Complex	z, q;
*/







//////////////////////////////////////////////////////////////////////
// Initialisation
//////////////////////////////////////////////////////////////////////

int CPertEngine::initialiseCalculateFrame(int WidthIn, int HeightIn, int threshold, BigDouble xZoomPointin, BigDouble yZoomPointin, double ZoomRadiusIn, int decimals /*, CTZfilter *TZfilter*/)
    {
    Complex q;

    width = WidthIn;
    height = HeightIn;
    MaxIteration = threshold;
    ZoomRadius = ZoomRadiusIn;
    
//    method = TZfilter->method;

    xZoomPt = xZoomPointin;
    yZoomPt = yZoomPointin;

/*
    if (method >= TIERAZONFILTERS)
	    {
	    q = { mpfr_get_d(xZoomPt, MPFR_RNDN), mpfr_get_d(yZoomPt, MPFR_RNDN) };
	    TZfilter->LoadFilterQ(q);		// initialise the constants used by Tierazon fractals
	    }
*/
    return 0;
    }

//////////////////////////////////////////////////////////////////////
// Full frame calculation
//////////////////////////////////////////////////////////////////////

int CPertEngine::calculateOneFrame(double bailout, char *StatusBarInfo, int powerin, int methodIn, int biomorphin, int subtypein, Complex rsrAin, bool rsrSignIn, int user_data(),
        void (*plot)(int, int, int) /*, int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol*/)

    {
    int i;
    BigComplex	C, ReferenceCoordinate;

    referencePoints = 0;
    GlitchPointCount = 0L;
    RemainingPointCount = 0L;

    rsrA = rsrAin;
    rsrSign = rsrSignIn;

    // get memory for all point arrays
    pointsRemaining = new Point[width * height];
    if (pointsRemaining == NULL)
	    return -1;
    glitchPoints = new Point[width * height];
    if (glitchPoints == NULL)
        {
	    if (pointsRemaining) { delete[] pointsRemaining; pointsRemaining = NULL; }
	    return -1;
        }
	// get memory for Perturbation Tolerance Check array
	PerturbationToleranceCheck = new double[MaxIteration * 2];
	if (PerturbationToleranceCheck == NULL)
	    {
	    if (pointsRemaining) { delete[] pointsRemaining; pointsRemaining = NULL; }
	    if (glitchPoints) { delete[] glitchPoints; glitchPoints = NULL; }
	    return -2;
	    }
	// get memory for Z array
	XSubN = new Complex[MaxIteration + 1];
	if (XSubN == NULL)
	    {
	    if (pointsRemaining) { delete[] pointsRemaining; pointsRemaining = NULL; }
	    if (glitchPoints) { delete[] glitchPoints; glitchPoints = NULL; }
	    if (PerturbationToleranceCheck) { delete[] PerturbationToleranceCheck; PerturbationToleranceCheck = NULL; }
	    return -2;
	    }
    biomorph = biomorphin;
    power = powerin;
    if (power < 2)
	    power = 2;
    if (power > MAXPOWER)
	    power = MAXPOWER;
/*                  keep this one for 'ron
    method = methodIn;
*/
    subtype = subtypein;

    // calculate the pascal's triangle coefficients for powers > 3
    LoadPascal(PascalArray, power);
    //Fill the list of points with all points in the image.
    for (long y = 0; y < height; y++) 
	    {
	    for (long x = 0; x < width; x++) 
	        {
	        Point pt(x, height - 1 - y);
	        pointsRemaining[y * width + x] = pt;
	        RemainingPointCount++;
	        }
	    }

    double magnifiedRadius = ZoomRadius;
    int window_radius = (width < height) ? width : height;

    while (RemainingPointCount > (width * height) * (percentGlitchTolerance / 100))
	    {
	    referencePoints++;
//	    if (user_data() == ID_KEY_ESC)
//	        {
//	        CloseTheDamnPointers();
//	        return -1;
//	        }

	    //Determine the reference point to calculate
	    //Check whether this is the first time running the loop. 
	    if (referencePoints == 1) 
	        {
            C.x = xZoomPt;
            C.y = yZoomPt;
            ReferenceCoordinate = C;
 
	        calculatedRealDelta = 0;
	        calculatedImaginaryDelta = 0;

	        if (ReferenceZoomPoint(ReferenceCoordinate, MaxIteration, user_data, StatusBarInfo) < 0)
		        {
		        CloseTheDamnPointers();
		        return -1;
		        }
	        }
	    else 
	        {
	        if (calculateGlitches == false) 
		        break;

	        int referencePointIndex = 0;
	        int Randomise;

	        srand((unsigned)time(NULL));		// Seed the random-number generator with current time 
	        Randomise = rand();
	        referencePointIndex = (int)((double)Randomise / (RAND_MAX + 1) * RemainingPointCount);

	        //Get the complex point at the chosen reference point
	        double deltaReal = ((magnifiedRadius * (2 * pointsRemaining[referencePointIndex].getX() - width)) / window_radius);
	        double deltaImaginary = ((-magnifiedRadius * (2 * pointsRemaining[referencePointIndex].getY() - height)) / window_radius);

	        // We need to store this offset because the formula we use to convert pixels into a complex point does so relative to the center of the image.
	        // We need to offset that calculation when our reference point isn't in the center. The actual offsetting is done in calculate point.

	        calculatedRealDelta = deltaReal;
	        calculatedImaginaryDelta = deltaImaginary;
	        ReferenceCoordinate.x = C.x + deltaReal;
	        ReferenceCoordinate.y = C.y + deltaImaginary;

	        if (ReferenceZoomPoint(ReferenceCoordinate, MaxIteration, user_data, StatusBarInfo) < 0)
		        {
		        CloseTheDamnPointers();
		        return -1;
		        }
	        }

	    int lastChecked = -1;
	    GlitchPointCount = 0;
	    for (i = 0; i < RemainingPointCount; i++)
	        {
            if (calculatePoint(pointsRemaining[i].getX(), pointsRemaining[i].getY(), magnifiedRadius, window_radius, bailout, glitchPoints, user_data, plot /*, potential, TZfilter, TrueCol*/) < 0)
		        return -1;
//            if (user_data() == ID_KEY_ESC)
//		        {
//		        CloseTheDamnPointers();
//		        return -1;
//		        }
	        //Everything else in this loop is just for updating the progress counter. 
	        double progress = (double)i / RemainingPointCount;
	        if (int(progress * 100) != lastChecked) 
		        {
		        lastChecked = int(progress * 100);
		        sprintf(StatusBarInfo, "Pass: %d, (%d%%)", referencePoints, int(progress * 100));
		        }
	        }

	    //These points are glitched, so we need to mark them for recalculation. We need to recalculate them using Pauldelbrot's glitch fixing method (see calculate point).
    //	memcpy(pointsRemaining, glitchPoints, sizeof(Point)*height*width);
	    memcpy(pointsRemaining, glitchPoints, sizeof(Point) * GlitchPointCount);
	    RemainingPointCount = GlitchPointCount;
	    }

    CloseTheDamnPointers();
    return 0;
    }

//////////////////////////////////////////////////////////////////////
// Cleanup
//////////////////////////////////////////////////////////////////////

void CPertEngine::CloseTheDamnPointers(void)
    {
    if (pointsRemaining) {delete[] pointsRemaining; pointsRemaining = NULL;}
    if (glitchPoints) {delete[] glitchPoints; glitchPoints = NULL;}
    if (XSubN) { delete[] XSubN; XSubN = NULL; }
	if (PerturbationToleranceCheck) { delete[] PerturbationToleranceCheck; PerturbationToleranceCheck = NULL; }
    }

//////////////////////////////////////////////////////////////////////
// Individual point calculation
//////////////////////////////////////////////////////////////////////

int CPertEngine::calculatePoint(int x, int y, double magnifiedRadius, int window_radius, double bailout, Point *glitchPoints, int user_data(), void (*plot)(int, int, int))
    //        int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol)
    {
// Get the complex number at this pixel.
// This calculates the number relative to the reference point, so we need to translate that to the center when the reference point isn't in the center.
// That's why for the first reference, calculatedRealDelta and calculatedImaginaryDelta are 0: it's calculating relative to the center.

    double deltaReal = ((magnifiedRadius * (2 * x - width)) / window_radius) - calculatedRealDelta;
    double deltaImaginary = ((-magnifiedRadius * (2 * y - height)) / window_radius) - calculatedImaginaryDelta;
    double	magnitude = 0.0;
    Complex DeltaSub0 {deltaReal, deltaImaginary};
    int iteration = skippedIterations;
    Complex DeltaSubN;
    DeltaSubN = DeltaSub0;
    iteration = 0;
    bool glitched = false;
//    int kbdchar;

    //if (method >= TIERAZONFILTERS)
	//TZfilter->LoadFilterQ(DeltaSub0);		// initialise the constants used by Tierazon filters

    //Iteration loop
    do
	    {
//        driver_wait_key_pressed(0);
//        kbdchar = driver_get_key();
//        if (kbdchar == ID_KEY_ESC)

//	    if (user_data() == ID_KEY_ESC)
//	        return -1;
        PertFunctions((XSubN + iteration), &DeltaSubN, &DeltaSub0);
	    iteration++;
        Complex CoordMag = *(XSubN + iteration) + DeltaSubN;
        ZCoordinateMagnitudeSquared = sqr(CoordMag.x) + sqr(CoordMag.y);

	    // This is Pauldelbrot's glitch detection method. You can see it here: http://www.fractalforums.com/announcements-and-news/pertubation-theory-glitches-improvement/.
	    // As for why it looks so weird, it's because I've squared both sides of his equation and moved the |ZsubN| to the other side to be precalculated.
	    // For more information, look at where the reference point is calculated.
	    // I also only want to store this point once.

    //	if (method >= TIERAZONFILTERS)
    //	    {
    //	    Complex z = XSubN[iteration] + DeltaSubN;
    //	    TZfilter->DoTierazonFilter(z, (long *)&iteration);
    //	    }

	    if (calculateGlitches == true && glitched == false && ZCoordinateMagnitudeSquared < PerturbationToleranceCheck[iteration])
	        {
	        Point pt(x, y, iteration);
	        glitchPoints[GlitchPointCount] = pt;
	        GlitchPointCount++;
	        glitched = true;
	        break;
	        }

	    //use bailout radius of 256 for smooth coloring.
	    } while (ZCoordinateMagnitudeSquared < bailout && iteration < MaxIteration);

        if (glitched == false) 
	        {
	        int	index;
	        double	rqlim2 = sqrt(bailout);
	        Complex	w = XSubN[iteration] + DeltaSubN;

	        if (biomorph >= 0)						// biomorph
	            {
	            if (iteration == MaxIteration)
		            index = MaxIteration;
	            else
		            {
		            if (fabs(w.x) < rqlim2 || fabs(w.y) < rqlim2)
		                index = biomorph;
		            else
		                index = iteration % 256;
		            }
	            }
	        else
	            {
    /*
	            switch (method)
		            {
		            case NONE:						// no filter
*/
		                if (iteration == MaxIteration)
			            index = MaxIteration;
		                else
			            index = iteration % 256;
/*
		            break;
		        case PERT1:						// something Shirom Makkad added
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = (int)((iteration - log2(log2(ZCoordinateMagnitudeSquared))) * 5) % 256; //Get the index of the color array that we are going to read from. 
		            break;
		        case PERT2:						// something Shirom Makkad added
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = (int)(iteration - (log(0.5*(ZCoordinateMagnitudeSquared)) - log(0.5*log(256))) / log(2)) % 256;
		            break;
		        case ZMAG:
		            if (iteration == MaxIteration)			// Zmag
			        index = (int)((CSumSqr(w)) * (MaxIteration >> 1) + 1);
		            else
			        index = iteration % 256;
		            break;
		        case REAL:						// "real"
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = iteration + (long)w.x + 7;
		            break;
		        case IMAG:	    					// "imag"
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = iteration + (long)w.y + 7;
		            break;
		        case MULT:						// "mult"
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else if (w.y)
			        index = (long)((double)iteration * (w.x / w.y));
		            break;
		        case SUM:						// "sum"
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = iteration + (long)(w.x + w.y);
		            break;
		        case ATAN:						// "atan"
		            if (iteration == MaxIteration)
			        index = MaxIteration;
		            else
			        index = (long)fabs(atan2(w.y, w.x)*180.0 / PI);
		            break;
		        case POTENTIAL:
		            magnitude = sqr(w.x) + sqr(w.y);
		            index = potential(magnitude, iteration);
		            break;
		        default:
		            if (method >= TIERAZONFILTERS)			// suite of Tierazon filters and colouring schemes
			        {
			        TZfilter->EndTierazonFilter(w, (long *)&iteration, TrueCol);
			        index = iteration;
			        }
		            else						// no filter
			        {
			        if (iteration == MaxIteration)
			            index = MaxIteration;
			        else
			            index = iteration % 256;
			        }
		            break;
		        }
*/
	        }
	    plot(x, height - 1 - y, index);
	    }
    return 0;
    }

//////////////////////////////////////////////////////////////////////
// Reference Zoom Point
//////////////////////////////////////////////////////////////////////

int CPertEngine::ReferenceZoomPoint(BigComplex& centre, int maxIteration, int user_data(), char* StatusBarInfo)
    {
    // Raising this number makes more calculations, but less variation between each calculation (less chance of mis-identifying a glitched point).
    BigComplex	ZTimes2, Z;
    double	glitchTolerancy = 1e-6;
    mpfr_t	TempReal, TempImag;
    BigDouble	zisqr, zrsqr, realimag;

    mpfr_init(TempReal);
    mpfr_init(TempImag);
    Z = centre;

    for (int i = 0; i <= maxIteration; i++)
	    {
	    // pre multiply by two
	    mpfr_mul_ui(ZTimes2.x.x, Z.x.x, 2, MPFR_RNDN);
	    mpfr_mul_ui(ZTimes2.y.x, Z.y.x, 2, MPFR_RNDN);
	    Complex c {Z.x.BigDoubleToDouble(), Z.y.BigDoubleToDouble()};
	    Complex TwoC {ZTimes2.x.BigDoubleToDouble(), ZTimes2.y.BigDoubleToDouble()};
	    // The reason we are storing the same value times two is that we can precalculate this value here because multiplying this value by two is needed many times in the program.
	    // Also, for some reason, we can't multiply complex numbers by anything greater than 1 using std::complex, so we have to multiply the individual terms each time.
	    // This is expensive to do above, so we are just doing it here.

	    XSubN[i] = c; 
	    // Norm is the squared version of abs and 0.000001 is 10^-3 squared.
	    // The reason we are storing this into an array is that we need to check the magnitude against this value to see if the value is glitched. 
	    // We are leaving it squared because otherwise we'd need to do a square root operation, which is expensive, so we'll just compare this to the squared magnitude.
	
	    if (user_data() < 0)
	        {
	        mpfr_clear(TempReal);
	        mpfr_clear(TempImag);
	        return -1;
	        }
	    //Everything else in this loop is just for updating the progress counter. 
	    int lastChecked = -1;
	    double progress = (double)i / maxIteration;
	    if (int(progress * 100) != lastChecked)
	        {
	        lastChecked = int(progress * 100);
	        sprintf(StatusBarInfo, "Pass: %d, Ref (%d%%)", referencePoints, int(progress * 100));
	        }

	    mpfr_mul_d(TempReal, Z.x.x, glitchTolerancy, MPFR_RNDN);
	    mpfr_mul_d(TempImag, Z.y.x, glitchTolerancy, MPFR_RNDN);
	    Complex tolerancy { mpfr_get_d(TempReal, MPFR_RNDN), mpfr_get_d(TempImag, MPFR_RNDN) };
//	    PerturbationToleranceCheck[i] = (CSumSqr(tolerancy));
        PerturbationToleranceCheck[i] = sqr(tolerancy.x) + sqr(tolerancy.y);

	    // Calculate the set
	    RefFunctions(&centre, &Z, &ZTimes2);
	    }
    mpfr_clear(TempReal);
    mpfr_clear(TempImag);
    return 0;
    }
    