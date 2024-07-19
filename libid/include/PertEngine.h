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

class CPertEngine 
    {
    public:
    int     initialiseCalculateFrame(int WidthIn, int HeightIn, int threshold, BigDouble xZoomPointin, BigDouble yZoomPointin, double ZoomRadiusIn, bool IsPotentialIn /*, CTZfilter *TZfilter*/);
    int     calculateOneFrame(double bailout, char *StatusBarInfo, int powerin, int InsideFilterIn, int OutsideFilterIn, int biomorph, int subtype, Complex RSRA, bool RSRsign, int user_data(),
                            void (*plot)(int, int, int), int potential(double, long)/*, CTZfilter *TZfilter, CTrueCol *TrueCol*/);
    private:
    int     calculatePoint(int x, int y, double tempRadius, int window_radius, double bailout,
            Point *glitchPoints, int user_data(), void (*plot)(int, int, int), int potential(double, long)/*, CTZfilter *TZfilter, CTrueCol *TrueCol*/);
    int     ReferenceZoomPoint(BigComplex *centre, int maxIteration, int user_data(), char *StatusBarInfo);
	void	LoadPascal(long PascalArray[], int n);
	double	DiffAbs(const double c, const double d);
	void	PertFunctions(Complex *XRef, Complex *DeltaSubN, Complex *DeltaSub0);
    void    RefFunctions(BigComplex *centre, BigComplex *Z, BigComplex *ZTimes2);
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
    bool    IsPotential = false;

	int	    width, height;
	int	    MaxIteration;
	int	    power, subtype, biomorph; 
	int	    InsideMethod;			// the number of the inside filter
	int	    OutsideMethod;			// the number of the outside filter
	long	GlitchPointCount;
	long	RemainingPointCount;
    mpfr_t	xZoomPt, yZoomPt;
	double	ZoomRadius;
	bool	calculateGlitches = true;
	bool	seriesApproximation = false;
	unsigned int numCoefficients = 5;
	//What percentage of the image is okay to be glitched. 
	double	percentGlitchTolerance = 0.1;
	int	    referencePoints = 0;
    int     kbdchar;            // keyboard character

    };


