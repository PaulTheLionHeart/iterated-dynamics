#include <float.h>

#include "port.h"
#include "prototyp.h"
#include "helpdefs.h"
#include "fractype.h"
#include "drivers.h"

// these need to be accessed elsewhere for saving data
double g_julibrot_x_min = -.83;
double g_julibrot_y_min = -.25;
double g_julibrot_x_max = -.83;
double g_julibrot_y_max =  .25;

static long mxmin, mymin;
static long x_per_inch, y_per_inch, inch_per_xdot, inch_per_ydot;
static double x_per_inchfp, y_per_inchfp, inch_per_xdotfp, inch_per_ydotfp;
static int bbase;
static long xpixel, ypixel;
static double xpixelfp, ypixelfp;
static long initz, djx, djy, dmx, dmy;
static double initzfp, djxfp, djyfp, dmxfp, dmyfp;
static long jx, jy, mx, my, xoffset, yoffset;
static double jxfp, jyfp, mxfp, myfp, xoffsetfp, yoffsetfp;

struct Perspective
{
    long x, y, zx, zy;
};

struct Perspectivefp
{
    double x, y, zx, zy;
};

Perspective LeftEye, RightEye, *Per;
Perspectivefp LeftEyefp, RightEyefp, *Perfp;

LComplex jbc;
DComplex jbcfp;

#ifndef XFRACT
static double fg, fg16;
#endif
int zdots = 128;

float g_julibrot_origin_fp  = 8.0F;
float g_julibrot_height_fp  = 7.0F;
float widthfp   = 10.0F;
float g_dist_fp    = 24.0F;
float g_eyes_fp    = 2.5F;
float g_julibrot_depth_fp   = 8.0F;
float brratiofp = 1.0F;
static long width, dist, depth, brratio;
#ifndef XFRACT
static long eyes;
#endif
int g_julibrot_3d_mode = 0;

fractal_type g_new_orbit_type = fractal_type::JULIA;

bool
JulibrotSetup()
{
    int r = 0;
    char const *mapname;

#ifndef XFRACT
    if (g_colors < 255)
    {
        stopmsg(STOPMSG_NONE, "Sorry, but Julibrots require a 256-color video mode");
        return false;
    }
#endif

    xoffsetfp = (xxmax + xxmin) / 2;     // Calculate average
    yoffsetfp = (yymax + yymin) / 2;     // Calculate average
    dmxfp = (g_julibrot_x_max - g_julibrot_x_min) / zdots;
    dmyfp = (g_julibrot_y_max - g_julibrot_y_min) / zdots;
    g_float_param = &jbcfp;
    x_per_inchfp = (xxmin - xxmax) / widthfp;
    y_per_inchfp = (yymax - yymin) / g_julibrot_height_fp;
    inch_per_xdotfp = widthfp / xdots;
    inch_per_ydotfp = g_julibrot_height_fp / ydots;
    initzfp = g_julibrot_origin_fp - (g_julibrot_depth_fp / 2);
    if (g_julibrot_3d_mode == 0)
    {
        RightEyefp.x = 0.0;
    }
    else
    {
        RightEyefp.x = g_eyes_fp / 2;
    }
    LeftEyefp.x = -RightEyefp.x;
    RightEyefp.y = 0;
    LeftEyefp.y = RightEyefp.y;
    RightEyefp.zx = g_dist_fp;
    LeftEyefp.zx = RightEyefp.zx;
    RightEyefp.zy = g_dist_fp;
    LeftEyefp.zy = RightEyefp.zy;
    bbase = 128;

#ifndef XFRACT
    if (fractalspecific[static_cast<int>(fractype)].isinteger > 0)
    {
        long jxmin, jxmax, jymin, jymax, mxmax, mymax;
        if (fractalspecific[static_cast<int>(g_new_orbit_type)].isinteger == 0)
        {
            stopmsg(STOPMSG_NONE, "Julibrot orbit type isinteger mismatch");
        }
        if (fractalspecific[static_cast<int>(g_new_orbit_type)].isinteger > 1)
        {
            bitshift = fractalspecific[static_cast<int>(g_new_orbit_type)].isinteger;
        }
        fg = (double)(1L << bitshift);
        fg16 = (double)(1L << 16);
        jxmin = (long)(xxmin * fg);
        jxmax = (long)(xxmax * fg);
        xoffset = (jxmax + jxmin) / 2;    // Calculate average
        jymin = (long)(yymin * fg);
        jymax = (long)(yymax * fg);
        yoffset = (jymax + jymin) / 2;    // Calculate average
        mxmin = (long)(g_julibrot_x_min * fg);
        mxmax = (long)(g_julibrot_x_max * fg);
        mymin = (long)(g_julibrot_y_min * fg);
        mymax = (long)(g_julibrot_y_max * fg);
        long origin = (long)(g_julibrot_origin_fp * fg16);
        depth = (long)(g_julibrot_depth_fp * fg16);
        width = (long)(widthfp * fg16);
        dist = (long)(g_dist_fp * fg16);
        eyes = (long)(g_eyes_fp * fg16);
        brratio = (long)(brratiofp * fg16);
        dmx = (mxmax - mxmin) / zdots;
        dmy = (mymax - mymin) / zdots;
        g_long_param = &jbc;

        x_per_inch = (long)((xxmin - xxmax) / widthfp * fg);
        y_per_inch = (long)((yymax - yymin) / g_julibrot_height_fp * fg);
        inch_per_xdot = (long)((widthfp / xdots) * fg16);
        inch_per_ydot = (long)((g_julibrot_height_fp / ydots) * fg16);
        initz = origin - (depth / 2);
        if (g_julibrot_3d_mode == 0)
        {
            RightEye.x = 0L;
        }
        else
        {
            RightEye.x = eyes / 2;
        }
        LeftEye.x = -RightEye.x;
        RightEye.y = 0L;
        LeftEye.y = RightEye.y;
        RightEye.zx = dist;
        LeftEye.zx = RightEye.zx;
        RightEye.zy = dist;
        LeftEye.zy = RightEye.zy;
        bbase = (int)(128.0 * brratiofp);
    }
#endif

    if (g_julibrot_3d_mode == 3)
    {
        savedac = 0;
        mapname = g_glasses1_map.c_str();
    }
    else
    {
        mapname = g_gray_map_file.c_str();
    }
    if (savedac != 1)
    {
        if (ValidateLuts(mapname))
        {
            return false;
        }
        spindac(0, 1);               // load it, but don't spin
        if (savedac == 2)
        {
            savedac = 1;
        }
    }
    return r >= 0;
}


int
jb_per_pixel()
{
    jx = multiply(Per->x - xpixel, initz, 16);
    jx = divide(jx, dist, 16) - xpixel;
    jx = multiply(jx << (bitshift - 16), x_per_inch, bitshift);
    jx += xoffset;
    djx = divide(depth, dist, 16);
    djx = multiply(djx, Per->x - xpixel, 16) << (bitshift - 16);
    djx = multiply(djx, x_per_inch, bitshift) / zdots;

    jy = multiply(Per->y - ypixel, initz, 16);
    jy = divide(jy, dist, 16) - ypixel;
    jy = multiply(jy << (bitshift - 16), y_per_inch, bitshift);
    jy += yoffset;
    djy = divide(depth, dist, 16);
    djy = multiply(djy, Per->y - ypixel, 16) << (bitshift - 16);
    djy = multiply(djy, y_per_inch, bitshift) / zdots;

    return (1);
}

int
jbfp_per_pixel()
{
    jxfp = ((Perfp->x - xpixelfp) * initzfp / g_dist_fp - xpixelfp) * x_per_inchfp;
    jxfp += xoffsetfp;
    djxfp = (g_julibrot_depth_fp / g_dist_fp) * (Perfp->x - xpixelfp) * x_per_inchfp / zdots;

    jyfp = ((Perfp->y - ypixelfp) * initzfp / g_dist_fp - ypixelfp) * y_per_inchfp;
    jyfp += yoffsetfp;
    djyfp = g_julibrot_depth_fp / g_dist_fp * (Perfp->y - ypixelfp) * y_per_inchfp / zdots;

    return (1);
}

static int zpixel, plotted;
static long n;

int
zline(long x, long y)
{
    xpixel = x;
    ypixel = y;
    mx = mxmin;
    my = mymin;
    switch (g_julibrot_3d_mode)
    {
    case 0:
    case 1:
        Per = &LeftEye;
        break;
    case 2:
        Per = &RightEye;
        break;
    case 3:
        if ((row + col) & 1)
        {
            Per = &LeftEye;
        }
        else
        {
            Per = &RightEye;
        }
        break;
    }
    jb_per_pixel();
    for (zpixel = 0; zpixel < zdots; zpixel++)
    {
        g_l_old_z.x = jx;
        g_l_old_z.y = jy;
        jbc.x = mx;
        jbc.y = my;
        if (driver_key_pressed())
        {
            return (-1);
        }
        g_l_temp_sqr_x = multiply(g_l_old_z.x, g_l_old_z.x, bitshift);
        g_l_temp_sqr_y = multiply(g_l_old_z.y, g_l_old_z.y, bitshift);
        for (n = 0; n < g_max_iterations; n++)
        {
            if (fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc())
            {
                break;
            }
        }
        if (n == g_max_iterations)
        {
            if (g_julibrot_3d_mode == 3)
            {
                g_color = (int)(128l * zpixel / zdots);
                if ((row + col) & 1)
                {

                    (*plot)(col, row, 127 - g_color);
                }
                else
                {
                    g_color = (int)(multiply((long) g_color << 16, brratio, 16) >> 16);
                    if (g_color < 1)
                    {
                        g_color = 1;
                    }
                    if (g_color > 127)
                    {
                        g_color = 127;
                    }
                    (*plot)(col, row, 127 + bbase - g_color);
                }
            }
            else
            {
                g_color = (int)(254l * zpixel / zdots);
                (*plot)(col, row, g_color + 1);
            }
            plotted = 1;
            break;
        }
        mx += dmx;
        my += dmy;
        jx += djx;
        jy += djy;
    }
    return (0);
}

int
zlinefp(double x, double y)
{
#ifdef XFRACT
    static int keychk = 0;
#endif
    xpixelfp = x;
    ypixelfp = y;
    mxfp = g_julibrot_x_min;
    myfp = g_julibrot_y_min;
    switch (g_julibrot_3d_mode)
    {
    case 0:
    case 1:
        Perfp = &LeftEyefp;
        break;
    case 2:
        Perfp = &RightEyefp;
        break;
    case 3:
        if ((row + col) & 1)
        {
            Perfp = &LeftEyefp;
        }
        else
        {
            Perfp = &RightEyefp;
        }
        break;
    }
    jbfp_per_pixel();
    for (zpixel = 0; zpixel < zdots; zpixel++)
    {
        // Special initialization for Mandelbrot types
        if ((g_new_orbit_type == fractal_type::QUATFP || g_new_orbit_type == fractal_type::HYPERCMPLXFP)
                && save_release > 2002)
        {
            g_old_z.x = 0.0;
            g_old_z.y = 0.0;
            jbcfp.x = 0.0;
            jbcfp.y = 0.0;
            qc = jxfp;
            qci = jyfp;
            qcj = mxfp;
            qck = myfp;
        }
        else
        {
            g_old_z.x = jxfp;
            g_old_z.y = jyfp;
            jbcfp.x = mxfp;
            jbcfp.y = myfp;
            qc = param[0];
            qci = param[1];
            qcj = param[2];
            qck = param[3];
        }
#ifdef XFRACT
        if (keychk++ > 500)
        {
            keychk = 0;
            if (driver_key_pressed())
            {
                return (-1);
            }
        }
#else
        if (driver_key_pressed())
        {
            return (-1);
        }
#endif
        tempsqrx = sqr(g_old_z.x);
        tempsqry = sqr(g_old_z.y);

        for (n = 0; n < g_max_iterations; n++)
        {
            if (fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc())
            {
                break;
            }
        }
        if (n == g_max_iterations)
        {
            if (g_julibrot_3d_mode == 3)
            {
                g_color = (int)(128l * zpixel / zdots);
                if ((row + col) & 1)
                {
                    (*plot)(col, row, 127 - g_color);
                }
                else
                {
                    g_color = (int)(g_color * brratiofp);
                    if (g_color < 1)
                    {
                        g_color = 1;
                    }
                    if (g_color > 127)
                    {
                        g_color = 127;
                    }
                    (*plot)(col, row, 127 + bbase - g_color);
                }
            }
            else
            {
                g_color = (int)(254l * zpixel / zdots);
                (*plot)(col, row, g_color + 1);
            }
            plotted = 1;
            break;
        }
        mxfp += dmxfp;
        myfp += dmyfp;
        jxfp += djxfp;
        jyfp += djyfp;
    }
    return (0);
}

int
Std4dFractal()
{
    long x;
    g_c_exponent = (int)param[2];
    if (g_new_orbit_type == fractal_type::LJULIAZPOWER)
    {
        if (g_c_exponent < 1)
        {
            g_c_exponent = 1;
        }
        if (param[3] == 0.0 && g_debug_flag != debug_flags::force_complex_power && (double)g_c_exponent == param[2])
        {
            fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc = longZpowerFractal;
        }
        else
        {
            fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc = longCmplxZpowerFractal;
        }
    }

    long y = 0;
    for (int ydot = (ydots >> 1) - 1; ydot >= 0; ydot--, y -= inch_per_ydot)
    {
        plotted = 0;
        x = -(width >> 1);
        for (int xdot = 0; xdot < xdots; xdot++, x += inch_per_xdot)
        {
            col = xdot;
            row = ydot;
            if (zline(x, y) < 0)
            {
                return (-1);
            }
            col = xdots - col - 1;
            row = ydots - row - 1;
            if (zline(-x, -y) < 0)
            {
                return (-1);
            }
        }
        if (plotted == 0)
        {
            if (y == 0)
            {
                plotted = -1;  // no points first pass; don't give up
            }
            else
            {
                break;
            }
        }
    }
    return (0);
}
int
Std4dfpFractal()
{
    double x;
    g_c_exponent = (int)param[2];

    if (g_new_orbit_type == fractal_type::FPJULIAZPOWER)
    {
        if (param[3] == 0.0 && g_debug_flag != debug_flags::force_complex_power && (double)g_c_exponent == param[2])
        {
            fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc = floatZpowerFractal;
        }
        else
        {
            fractalspecific[static_cast<int>(g_new_orbit_type)].orbitcalc = floatCmplxZpowerFractal;
        }
        get_julia_attractor(param[0], param[1]);  // another attractor?
    }

    double y = 0.0;
    for (int ydot = (ydots >> 1) - 1; ydot >= 0; ydot--, y -= inch_per_ydotfp)
    {
        plotted = 0;
        x = -widthfp / 2;
        for (int xdot = 0; xdot < xdots; xdot++, x += inch_per_xdotfp)
        {
            col = xdot;
            row = ydot;
            if (zlinefp(x, y) < 0)
            {
                return (-1);
            }
            col = xdots - col - 1;
            row = ydots - row - 1;
            if (zlinefp(-x, -y) < 0)
            {
                return (-1);
            }
        }
        if (plotted == 0)
        {
            if (y == 0)
            {
                plotted = -1;  // no points first pass; don't give up
            }
            else
            {
                break;
            }
        }
    }
    return (0);
}
