#pragma once
#include <windows.h>
#include <stdio.h>
#include "Big.h"
#include "Point.h"
#include "filter.h"
#include "Complex.h"
#include "id_keys.h"
/*
#include "Dib.h"
#include "colour.h"
#include "BigDouble.h"
#include "BigComplex.h"
*/

#define	MAXPOWER    28
#define MAXFILTER 9



// a few norty externals for the moment

extern int save_stack();                   // Returns stack pointer offset so it can be saved.
extern void restore_stack(int old_offset); // Restores stack pointer, effectively freeing local variables allocated since save_stack()
extern bn_t alloc_stack(size_t size);      // Allocates a bn_t variable on stack



class CPertEngine 
    {
    public:
        int initialiseCalculateFrame(int WidthIn, int HeightIn, int threshold, bf_t xZoomPointin, bf_t yZoomPointin, double ZoomRadiusIn, int decimals /*, CTZfilter *TZfilter */);
	int	calculateOneFrame(double bailout, char* StatusBarInfo, int powerin, int FilterTypeIn, int biomorph, int subtype, Complex RSRA, bool RSRsign, int user_data(),
                            void (*plot)(int, int, int) /*, int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol*/);
    private:
        int calculatePoint(int x, int y, double tempRadius, int window_radius, double bailout,
            Point *glitchPoints, int user_data(), void (*plot)(int, int, int) /*, int potential(double, int), CTZfilter *TZfilter, CTrueCol *TrueCol*/);
        int ReferenceZoomPoint(bf_t *centreX, bf_t *centreY, int MAX_ITERATIONS, int user_data(), char *StatusBarInfo);
	void	LoadPascal(long PascalArray[], int n);
	double	DiffAbs(const double c, const double d);
	void	PertFunctions(Complex *XRef, Complex *DeltaSubN, Complex *DeltaSub0);
    void    RefFunctions(bf_t *centreX, bf_t *centreY, bf_t *zX, bf_t *zY, bf_t *ZTimes2X, bf_t *ZTimes2Y);
	void	CloseTheDamnPointers(void);

	Complex	*XSubN;
	double	*PerturbationToleranceCheck;
	double	calculatedRealDelta, calculatedImaginaryDelta;
	double	ZCoordinateMagnitudeSquared;
	int	    skippedIterations = 0;

	long	PascalArray[MAXPOWER];
	Point	*pointsRemaining;
	Point	*glitchPoints;

//	Complex	q;			// location of current pixel
	Complex	rsrA;			// TheRedshiftRider value of a
	bool	rsrSign;		// TheRedshiftRider sign true if positive

	int	    width, height;
	int	    MaxIteration;
	int	    power, subtype, method, biomorph; 
	long	GlitchPointCount;
	long	RemainingPointCount;
    bf_t    xZoomPt, yZoomPt;
	double	ZoomRadius;
	bool	calculateGlitches = true;
	bool	seriesApproximation = false;
	unsigned int numCoefficients = 5;
	//What percentage of the image is okay to be glitched. 
	double	percentGlitchTolerance = 0.1;
	int	referencePoints = 0;
    int saved;              // saved value of bignum stack

    };


