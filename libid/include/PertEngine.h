#pragma once
#include <windows.h>
#include <stdio.h>
//#include "Big.h"
#include "Point.h"
#include "filter.h"
#include "Complex.h"
#include "id_keys.h"
/*
#include "Dib.h"
#include "colour.h"
*/
#include "BigDouble.h"
#include "BigComplex.h"

#define	MAXPOWER    28
#define MAXFILTER   9
#define ALLOW_MPFR

class CPertEngine 
    {
    public:
#ifdef ALLOW_MPFR
    int initialiseCalculateFrame(int WidthIn, int HeightIn, int threshold, BigDouble xZoomPointin, BigDouble yZoomPointin, double ZoomRadiusIn, int decimals /*, CTZfilter *TZfilter*/);
#else
    int initialiseCalculateFrame(int WidthIn, int HeightIn, int threshold, double xZoomPointin, double yZoomPointin, double ZoomRadiusIn, int decimals /*, CTZfilter *TZfilter*/);
#endif // ALLOW_MPFR
    int calculateOneFrame(double bailout, char *StatusBarInfo, int powerin, int FilterTypeIn, int biomorph, int subtype, Complex RSRA, bool RSRsign, int user_data(),
                            void (*plot)(int, int, int) /*, int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol*/);
    private:
    int calculatePoint(int x, int y, double tempRadius, int window_radius, double bailout,
            Point *glitchPoints, int user_data(), void (*plot)(int, int, int) /*, int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol*/);
#ifdef ALLOW_MPFR
    int     ReferenceZoomPoint(BigComplex *centre, int maxIteration, int user_data(), char *StatusBarInfo);
#else
    int     ReferenceZoomPoint(Complex *centre, int maxIteration, int user_data(), char *StatusBarInfo);
#endif // ALLOW_MPFR
	void	LoadPascal(long PascalArray[], int n);
	double	DiffAbs(const double c, const double d);
	void	PertFunctions(Complex *XRef, Complex *DeltaSubN, Complex *DeltaSub0);
#ifdef ALLOW_MPFR
    void RefFunctions(BigComplex *centre, BigComplex *Z, BigComplex *ZTimes2);
#else
    void RefFunctions(Complex *centre, Complex *Z, Complex *ZTimes2);
#endif // ALLOW_MPFR
	void	CloseTheDamnPointers(void);

	Complex *XSubN = NULL;
	double	*PerturbationToleranceCheck = NULL;
	double	calculatedRealDelta, calculatedImaginaryDelta;
	double	ZCoordinateMagnitudeSquared;
	int	    skippedIterations = 0;

	long	PascalArray[MAXPOWER];
	Point	*pointsRemaining = NULL;
	Point	*glitchPoints = NULL;

//	Complex	q;			// location of current pixel
	Complex	rsrA;			// TheRedshiftRider value of a
	bool	rsrSign;		// TheRedshiftRider sign true if positive

	int	    width, height;
	int	    MaxIteration;
	int	    power, subtype, method, biomorph; 
	long	GlitchPointCount;
	long	RemainingPointCount;
#ifdef ALLOW_MPFR
    BigDouble    xZoomPt, yZoomPt;
#else
    double    xZoomPt, yZoomPt;
#endif // ALLOW_MPFR
	double	ZoomRadius;
	bool	calculateGlitches = true;
	bool	seriesApproximation = false;
	unsigned int numCoefficients = 5;
	//What percentage of the image is okay to be glitched. 
	double	percentGlitchTolerance = 0.1;
	int	referencePoints = 0;
    int saved;              // saved value of bignum stack

    };


