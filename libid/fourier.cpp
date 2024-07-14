/*
    FOURIER.CPP a module for Fourier Synthesis of signals
    
    Written in Microsoft Visual C++ by Paul de Leeuw.
*/

#include	<conio.h>
#include	<string.h>
#include	<stdio.h>
#include	<math.h>
#include	<windows.h>
#include	<windowsx.h>
#include	"resource.h"
#include	"fractype.h"
#include	"id.h"
#include    "drivers.h"
#include    "read_ticker.h"

#define	PREVIEW_HEIGHT	168
#define	PREVIEW_WIDTH	180
#define	HOR_OFFSET	145
#define	VERT_OFFSET	142

#define	DEG2RAD		57.29577951
#define	FOURIERMAX	1000
#define	MAXSTEPS	2500
#define	LEVELS		30		// number of harmonic sliders

#define	SQUARE		0
#define	TRIANGLE	1
#define	SAWTOOTH	2
#define	FULLWAVE	3
#define	SINEWAVE	4
#define	COSINEWAVE	5
#define	IMPULSE		6
#define	USER		'U'     // can be implemented when we get a better GUI
#define HALF_PI     PI/2
#define TWO_PI      PI*2
#define BYTE        UINT8
#define WORD        UINT16
#define DWORD       UINT32;
// Macro to determine to round off the given value to the closest byte
#define WIDTHBYTES(i) ((i + 31) / 32 * 4)

struct	fourierstruct
    {    
    double	x;				        // starting point for harmonic
    double	y;
    double	magsin;				    // magnitude of the sine component of the harmonic
    double	magcos;				    // magnitude of the cosine component of the harmonic
    unsigned short	c;				// colour of harmonic
    }	FourierArray[FOURIERMAX];   

extern	long	g_max_iterations;
extern	int	g_row, g_col;
extern  int     g_screen_x_dots, g_screen_y_dots;
extern  double  g_params[];

int	WaveformArray[MAXSTEPS];   

static  int	    WavePtr = 0;	    // point to the currently updated location
static	short	level[LEVELS];
static  bool	MovingWave = true;
static  int     time_to_break;      // time to break out of animation?
static bool     first = true;
static int      NumHarmonics = 60, steps = 400, delay = 20;
static double   MagX, MagY;         // magnitude of the fundamental expressed in screen size

extern int      driver_key_pressed();

/**************************************************************************
	Plot Harmonics
***************************************************************************/

int	Harmonics(bool clear, int TotalFrames)

    {
    int	i, x1, y1, x2, y2, centrex, centrey;
    double	xold, yold;
    double	xnew, ynew;
    int	test;

    xnew = ynew = xold = yold = 0.0;
    centrey = g_screen_y_dots / 2;
    centrex = g_screen_x_dots / 4;

    for (i = 0; i < NumHarmonics; ++i)
	    {
	    if (FourierArray[i].magsin == 0.0 && FourierArray[i].magcos == 0.0)		// nothing to do
	        continue;
	    xnew += FourierArray[i].x;
	    ynew += FourierArray[i].y;
	    test = ((int)(MagX * xold) + centrex);
	    x1 = (test < 0) ? 0 : test;
	    test = ((int)(MagY * yold) + centrey);
	    y1 = (test < 0) ? 0 : test;
	    test = ((int)(MagX * xnew) + centrex);
	    x2 = (test < 0) ? 0 : test;
	    test = ((int)(MagY * ynew) + centrey);
	    y2 = (test < 0) ? 0 : test;
	    if (x1 > g_screen_x_dots / 2 - 10)					// make sure we don't write a harmonic into the history area
	        x1 = g_screen_x_dots / 2 - 10;
	    if (x2 > g_screen_x_dots / 2 - 10)
	        x2 = g_screen_x_dots / 2 - 10;
	    if (x1 < 0)	
	        x1 = 0;
	    if (x2 < 0)
	        x2 = 0;
	    if (y1 > g_screen_y_dots - 1)						// check boundaries
	        y1 = g_screen_y_dots - 1;
	    if (y2 > g_screen_y_dots - 1)
	        y2 = g_screen_y_dots - 1;
	    if (y1 < 0)	
	        y1 = 0;
	    if (y2 < 0)
	        y2 = 0;
	    driver_draw_line(x1, (WORD)(g_screen_y_dots - y1), x2, (WORD)(g_screen_y_dots - y2), clear ? 0 : FourierArray[i].c);
	    xold = xnew;
	    yold = ynew;
	    }
    // load waveform array with the current
    if (!clear)
	    {
	    WaveformArray[WavePtr++] = y2;
	    if (WavePtr >= TotalFrames)
	        WavePtr = 0;			// cycle through again
	    }
    return 0;
    }

/**************************************************************************
	Set up table of values for the starting point of each harmonic
***************************************************************************/

void	CalculateFourier(double angle)

    {
    int	i;

    for (i = 1; i < NumHarmonics; i++)
	    {
	    if (FourierArray[i].magsin == 0.0 && FourierArray[i].magcos == 0.0)		// nothing to do
	        continue;
	    FourierArray[i].x = -sin(HALF_PI + i * angle) * FourierArray[i].magsin + sin(i * angle) * FourierArray[i].magcos;
	    FourierArray[i].y = -cos(HALF_PI + i * angle) * FourierArray[i].magsin + cos(i * angle) * FourierArray[i].magcos;
	    FourierArray[i].c = i * (int)g_max_iterations / NumHarmonics + 1;
	    }
    }

/**************************************************************************
	Set up screen scaling and other inits
***************************************************************************/

void	InitFourier(void)

    {
    int i, WaveType;
    double	sign;

    WaveType = (int)g_params[0];
    if (WaveType < 0 || WaveType > 6)
        WaveType = 0;
    delay = (int) g_params[1];
    if (delay < 0 || delay > 5000)
        delay = 100;
    NumHarmonics = (int) g_params[2];
    if (NumHarmonics < 1 || NumHarmonics >= FOURIERMAX)
        NumHarmonics = 60;
    steps = (int) g_params[3];
    if (steps < 1 || steps >= MAXSTEPS)
        steps = 400;

    WavePtr = 0;							// point to the currently updated location
    MagX = (double) g_screen_x_dots / 16.0;
    MagY = (double) g_screen_y_dots / 8.0;
    if (WaveType != USER)						// we don't want to do this if we are using sliders
	    {
	    for (i = 0; i < FOURIERMAX; i++)
	        {
	        FourierArray[i].magsin = 0.0;
	        FourierArray[i].magcos = 0.0;
	        }

	    if (WaveType == IMPULSE || WaveType == SAWTOOTH || WaveType == SQUARE || WaveType == SINEWAVE)
	        FourierArray[1].magsin = 1.0;
	    else
	        FourierArray[1].magcos = 1.0;

	    sign = (WaveType == IMPULSE || WaveType == SAWTOOTH) ? -1.0 : 1.0;

	    if (WaveType == COSINEWAVE)
	        FourierArray[0].magcos = 1.0;				// only fundamental
	    if (WaveType == SINEWAVE)
	        FourierArray[0].magsin = 1.0;				// only fundamental
	    }

    for (i = 1; i < NumHarmonics; i++)
	    {
	    switch (WaveType)
	        {
	        case IMPULSE:
		        if (i % 2 == 0)					// reject even values
		            continue;
		        sign *= -1.0;					// change sign
		        FourierArray[i].magsin = sign / (PI * 5.0);	// all the same
		        break;
	        case SAWTOOTH:
		        sign *= -1.0;					// change sign
		        FourierArray[i].magsin = sign / (double)i;
		        break;
	        case SQUARE:
		        if (i % 2 == 0)					// reject even values
		            continue;
		        FourierArray[i].magsin = 1.0 / (double)i;
		        break;
	        case TRIANGLE:
		        if (i % 2 == 0)					// reject even values
		            continue;
		        FourierArray[i].magcos = 1.0 / (double)(i * i);	// n squared
		        break;
	        case FULLWAVE:
		        if (i == 1)					    // reject first harmonic  - no divide by 0
		            {
		            FourierArray[i].magcos = 0.0;
		            continue;
		            }
		        if (i % 2 == 1)					// reject odd values
		            continue;
		        sign *= -1.0;					// change sign
		            FourierArray[i].magcos = sign * 2.0 / (1.0 - (double)(i * i));
		        break;
	        }
	    }
    for (i = 0; i < LEVELS; i++)
	    level[i] = (i < LEVELS / 2) ? (short)(FourierArray[i + 1].magsin * 255.0) : (short)(FourierArray[i - LEVELS / 2 + 1].magcos * 255.0);
    first = TRUE;						// don't draw artificats when starting up
    }

/**************************************************************************
	Fourier Images
***************************************************************************/

int	FourierStep(int TotalFrames, int ThisStep)

    {
    static	double	angle;
    int		k, m, index1, index2;
    static	int	HorPos, OldHorPos;
    BYTE	*buffer = nullptr;
    long	address;
    static	int	colour = 0;
    bool    MovingWave = true;

    UINT32  tick;

    if (driver_key_pressed())
       return -1;
                     
    buffer = new BYTE[g_screen_x_dots * 3];     // ready for 24 bit true colour
    if (buffer == nullptr)
        return -1;

    tick = readticker();
    while (readticker() < tick + delay);

    MovingWave = (g_params[4] != 0.0);
    OldHorPos = HorPos;
    HorPos = (UINT32) g_screen_x_dots * (UINT32) ThisStep / (UINT32) (TotalFrames * 2) + g_screen_x_dots / 2;
    Harmonics(TRUE, TotalFrames);		// clear previous vectors
    angle = (TWO_PI / (double)TotalFrames) * (double) ThisStep;
    CalculateFourier(angle);
    Harmonics(FALSE, TotalFrames);		// write new vectors
    driver_flush();
    index1 = WavePtr - 1;
    index2 = WavePtr - 2;
    if (index1 < 0)
	    index1 += TotalFrames;
    if (index2 < 0)
	    index2 += TotalFrames;
    
    if (first)
	    first = FALSE;
    else						// don't draw artificats when starting up
	    {

	    if (MovingWave)
	        {
	        for (k = 0; k < g_screen_y_dots; k++)
		        {
//		        address = WIDTHBYTES((UINT32)g_screen_x_dots * (UINT32)bits_per_pixel) * k + (g_screen_x_dots / 2 - 3) * 3;
//		        memcpy(buffer, Dib.DibPixels + address, (g_screen_x_dots / 2) * 3);
//		        memcpy(Dib.DibPixels + address + 3, buffer, (g_screen_x_dots / 2 - 1) * 3);
                driver_read_span(k, g_screen_x_dots / 2 - 1, g_screen_x_dots - 1, buffer);
                driver_write_span(k, g_screen_x_dots / 2, g_screen_x_dots - 1, buffer);
                }

	        driver_draw_line((WORD)(g_screen_x_dots / 2), (WORD)(g_screen_y_dots - WaveformArray[index2]), 
					        (WORD)(g_screen_x_dots / 2 + 1), (WORD)(g_screen_y_dots - WaveformArray[index1]), colour);
	        }
	    else
	        driver_draw_line((WORD) (OldHorPos - 2), (WORD) (g_screen_y_dots - WaveformArray[index2]), 
		    (WORD)((HorPos > OldHorPos) ? (HorPos - 2) : OldHorPos),	// remove retrace
		    (WORD)(g_screen_y_dots - WaveformArray[index1]), colour);
	    }

    if (colour++ >= g_max_iterations)                   // rotate through palette
        colour = 1;
    if (buffer)  {delete[] buffer; buffer = nullptr;}

    return 0;
    }                                                 

/**************************************************************************
	Fourier Images
***************************************************************************/

bool	Fourier(void)

    {
    int	j;

//    DisplayAnimation = TRUE;
//    time_to_break = FALSE;
    InitFourier();

    while (time_to_break == FALSE)
	    {
	    for (j = 0; j < steps; ++j)
	        {
	        if (time_to_break)
		        break;
	        if (FourierStep(steps, j) < 0)
		        return 0;
	        }
	    }

    return 0;
    }

