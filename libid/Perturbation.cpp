/****************************************************
    PERTURBATION.CPP a module to explore Perturbation
    Thanks to Shirom Makkad fractaltodesktop@gmail.com
    Written in Microsoft Visual C++ by Paul de Leeuw.
****************************************************/

#include "fractalp.h"
#include "PertEngine.h"
//#include "Big.h"

/*
#include "manpwin.h"
#include "BigDouble.h"
#include "colour.h"
#include "resource.h"
#include "dib.h"
#include "filter.h"
#include "Fract1.h"
*/

char* xZoomPointString = "-1.25736802846652839265383159384773654166836713857126000896912753375688559878664765114255696457015368246531973104439755978333044015506759938503739206829441575363669402497147343368904702066174408250247081855416385744218741909521990441308969603994513271641284298225323509381146075334937409491188";
char* yZoomPointString = "0.37873083102862491151257052392242106932532193327534173605649141946411779667848532042309666819671311329800095959763206221251386819369531602358854394169140220049675504811341950346171196600590463661845947574424944950533273158558278821586333530950155398785389980291835095955110139525825547062070";
//const char* xZoomPointString = "-1.9428";
//const char* yZoomPointString = "0.0";
//const char* xZoomPointString = "-0.75";
//const char* yZoomPointString = "0.01";
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

//extern	void	ConvertBignum2String(char *s, mpfr_t num);
extern int getakeynohelp();
extern void (*g_plot)(int, int, int); // function pointer
//extern	int	potential(double, int);
extern  int g_screen_x_dots, g_screen_y_dots; // # of dots on the physical screen

extern	long	g_max_iterations;
extern  double  g_x_max, g_x_min, g_y_max, g_y_min;

	    double	hor;				// horizontal address
	    double	vert;				// vertical address
	    double	mandel_width;		// width of display
//extern	BigDouble   BigHor, BigVert, BigWidth;
extern	BYTE	BigNumFlag;			// True if bignum used
//extern	DWORD 	*wpixels;			// the device-independent bitmap  pixels
extern	char	szStatus[];			// status bar text
extern double g_magnitude_limit; // bailout level
//extern	int	decimals;
extern	double	param[];
//extern	int	method;				// inside and outside filters
//extern	BYTE	degree;				// power
//extern	int	biomorph;			// biomorph colour
//extern	RGBTRIPLE FilterRGB;			// for Tierazon filters
//extern	double	dStrands;
//extern	BOOL	UseCurrentPalette;		// do we use the ManpWIN palette? If false, generate internal filter palette

//////////////////////////////////////////////////////////////////////////////////////////////////////
extern	HWND	GlobalHwnd;			// This is the main windows handle
//////////////////////////////////////////////////////////////////////////////////////////////////////

/*
extern	CTrueCol    TrueCol;			// palette info
extern	CDib	Dib;
extern	CFract	Fractal;			// current fractal stuff
static	CTZfilter	TZfilter;		// Tierazon filters
*/
CPertEngine PertEngine;
char		PertStatus[200];

extern void init_bf_dec(int dec);
    //static	int	PerturbationPtr = 0, PerturbationNum = 0;

/**************************************************************************
	Initialise Perturbation engine
**************************************************************************/

bool	InitPerturbation(void)
    {
    char	*HorString = NULL;
    char	*VertString = NULL;
    int saved;
    saved = save_stack();

    bf_math = bf_math_type::BIGFLT;

    int     Sizeof_BF_Vars = strlen_needed_bf();
    bf_t    Big_centrex, Big_centrey;

    HorString = new char[Sizeof_BF_Vars];
    VertString = new char[Sizeof_BF_Vars];

    init_bf_dec(24);                                        // note this is only for testing
    mandel_width = g_x_max - g_x_min;
//    Height = g_y_max - g_y_min;

    /*
        if (!BigNumFlag)
	    {
	    BigHor = hor;
	    BigVert = vert;
	    BigWidth = mandel_width;
	    }
    Big_centrex = BigHor + (BigWidth * ((double)Dib.DibWidth / (double)(2 * Dib.DibHeight)));
    Big_centrey = -(BigVert + (BigWidth / 2.0));						// negative because of screen layout
    ConvertBignum2String(HorString, Big_centrex.x);
    ConvertBignum2String(VertString, Big_centrey.x);
    if (method >= TIERAZONFILTERS)
	TZfilter.InitFilter(method, threshold, dStrands, &FilterRGB, UseCurrentPalette);		// initialise the constants used by Tierazon fractals
    if (BigNumFlag)
	mandel_width = mpfr_get_d(BigWidth.x, MPFR_RNDN);
*/
    PertEngine.initialiseCalculateFrame(g_screen_x_dots, g_screen_x_dots, g_max_iterations, /*HorString, VertString*/xZoomPointString, yZoomPointString, mandel_width / 2, bflength/*, &TZfilter*/);
    if (HorString) { delete[] HorString; HorString = NULL; }
    if (VertString) { delete[] VertString; VertString = NULL; }
    return false;
    }

/**************************************************************************
	The Perturbation engine
**************************************************************************/

int	DoPerturbation(void)
    {
    int (*UserData)() = getakeynohelp;
    Complex a = {0, 0};
    bool    IsPositive = false;

/*    
    if (PerturbationSpecific[subtype].RedShiftRider == 1)
	{
	a.x = param[0];
	a.y = param[1];
	degree = (int)param[2];
	IsPositive = (param[3] == 1.0);
	}
    else if (param[0] != 0.0)
	degree = (int)param[0];
*/
    if (PertEngine.calculateOneFrame(g_magnitude_limit, PertStatus, /*degree, method, biomorph, subtype*/2, 0, 0, 0, a, IsPositive, UserData, g_plot/*, potential, &TZfilter, &TrueCol*/) < 0)
        //    if (frameCalculator.calculateOneFrame(rqlim, PertStatus, degree, method, biomorph, subtype, a, IsPositive, UserData, plot, potential) < 0)
	    return -1;
    return 0;
    }

/**************************************************************************
    Some simple processing
**************************************************************************/

/*
int	setup_Perturbation(void)
    {
    if (!PerturbationNum)	    // we'd better count how many records we have
	{
	while (PerturbationSpecific[PerturbationNum].name)
	    PerturbationNum++;
	}
    return 0;
    }

void	LoadPerturbationParams(void)
    {
    // here is where we can do some specific updates to individual Perturbation fractals
    int	i;

    Fractal.NumParam = PerturbationSpecific[subtype].numparams;

    for (i = 0; i < NUMPERTPARAM; i++)
	{
	param[i] = PerturbationSpecific[subtype].paramvalue[i];
	Fractal.ParamName[i] = PerturbationSpecific[subtype].paramname[i];
	}
    rqlim = PerturbationSpecific[subtype].rqlim;
    }
*/

/**************************************************************************
    Show Perturbation Fractal
**************************************************************************/

/*
DLGPROC FAR PASCAL SelectPertDlg(HWND hDlg, unsigned message, WPARAM wParam, LPARAM lParam)

    {
    static	int	i;
    static	int	index = 1;

    switch (message)
	{
	case WM_INITDIALOG:
	    SetDlgItemText(hDlg, ID_LISTTITLE, "Perturbation");
	    for (i = 0; i < PerturbationNum; i++)
		SendDlgItemMessage(hDlg, IDM_LSYSTEM, LB_ADDSTRING, (WPARAM)NULL, (LPARAM)(LPSTR)PerturbationSpecific[i].name);
	    SendDlgItemMessage(hDlg, IDM_LSYSTEM, LB_SETCURSEL, (WPARAM)PerturbationPtr, 0L);
	    return ((DLGPROC)TRUE);

	case WM_COMMAND:
	    switch ((int)LOWORD(wParam))
		{
		case IDOK:
		    index = (int)SendDlgItemMessage(hDlg, IDM_LSYSTEM, LB_GETCURSEL, 0, 0L);
		    if (index == LB_ERR)
			{
			MessageBox(hDlg, "No Choice selected",
			    "Select From a List", MB_OK | MB_ICONEXCLAMATION);
			break;
			}
		    PerturbationPtr = index;
		    subtype = PerturbationPtr;
		    EndDialog(hDlg, TRUE);
		    return ((DLGPROC)TRUE);

		case IDCANCEL:
		    PerturbationPtr = 1;
		    EndDialog(hDlg, FALSE);
		    return ((DLGPROC)FALSE);

		case IDM_LSYSTEM:
		    switch (HIWORD(wParam) & 0x0003)
			{
			case LBN_SELCHANGE:
			    index = (int)SendDlgItemMessage(hDlg, IDM_LSYSTEM, LB_GETCURSEL, 0, 0L);
			    if (index == LB_ERR)
				break;
			    break;

			case LBN_DBLCLK:
			    index = (int)SendDlgItemMessage(hDlg, IDM_LSYSTEM, LB_GETCURSEL, 0, 0L);
			    if (index == LB_ERR)
				{
				MessageBox(hDlg, "No Choice selected",
				    "Select From a List", MB_OK | MB_ICONEXCLAMATION);
				break;
				}
			    subtype = PerturbationPtr = index;
			    EndDialog(hDlg, TRUE);
			    return ((DLGPROC)TRUE);
			}
		    return ((DLGPROC)TRUE);
		}
	}
    return ((DLGPROC)FALSE);
    }
*/
