#include "newton.h"

#include "port.h"
#include "prototyp.h"

#include "calcfrac.h"
#include "cmplx.h"
#include "fractals.h"
#include "fractype.h"
#include "pixel_grid.h"

#include <cmath>

#define distance(z1, z2)  (sqr((z1).x-(z2).x)+sqr((z1).y-(z2).y))

static double t2{};

// this code translated to asm - lives in newton.asm
// transform points with reciprocal function
void invertz2(DComplex *z)
{
    z->x = g_dx_pixel();
    z->y = g_dy_pixel();
    z->x -= g_f_x_center;
    z->y -= g_f_y_center;  // Normalize values to center of circle

    g_temp_sqr_x = sqr(z->x) + sqr(z->y);  // Get old radius
    if (std::fabs(g_temp_sqr_x) > FLT_MIN)
    {
        g_temp_sqr_x = g_f_radius / g_temp_sqr_x;
    }
    else
    {
        g_temp_sqr_x = FLT_MAX;   // a big number, but not TOO big
    }
    z->x *= g_temp_sqr_x;
    z->y *= g_temp_sqr_x;      // Perform inversion
    z->x += g_f_x_center;
    z->y += g_f_y_center; // Renormalize
}

// Distance of complex z from unit circle
#define DIST1(z) (((z).x-1.0)*((z).x-1.0)+((z).y)*((z).y))
#define LDIST1(z) (lsqr((((z).x)-g_fudge_factor)) + lsqr(((z).y)))

int NewtonFractal2()
{
    static char start = 1;
    if (start)
    {
        start = 0;
    }
    cpower(&g_old_z, g_degree-1, &g_tmp_z);
    complex_mult(g_tmp_z, g_old_z, &g_new_z);

    if (DIST1(g_new_z) < g_threshold)
    {
        if (g_fractal_type == fractal_type::NEWTBASIN || g_fractal_type == fractal_type::MPNEWTBASIN)
        {
            long tmpcolor;
            tmpcolor = -1;
            /* this code determines which degree-th root of root the
               Newton formula converges to. The roots of a 1 are
               distributed on a circle of radius 1 about the origin. */
            for (int i = 0; i < g_degree; i++)
            {
                /* color in alternating shades with iteration according to
                   which root of 1 it converged to */
                if (distance(g_roots[i], g_old_z) < g_threshold)
                {
                    if (g_basin == 2)
                    {
                        tmpcolor = 1+(i&7)+((g_color_iter&1) << 3);
                    }
                    else
                    {
                        tmpcolor = 1+i;
                    }
                    break;
                }
            }
            if (tmpcolor == -1)
            {
                g_color_iter = g_max_color;
            }
            else
            {
                g_color_iter = tmpcolor;
            }
        }
        return 1;
    }
    g_new_z.x = g_degree_minus_1_over_degree * g_new_z.x + g_newton_r_over_d;
    g_new_z.y *= g_degree_minus_1_over_degree;

    // Watch for divide underflow
    t2 = g_tmp_z.x*g_tmp_z.x + g_tmp_z.y*g_tmp_z.y;
    if (t2 < FLT_MIN)
    {
        return 1;
    }
    else
    {
        t2 = 1.0 / t2;
        g_old_z.x = t2 * (g_new_z.x * g_tmp_z.x + g_new_z.y * g_tmp_z.y);
        g_old_z.y = t2 * (g_new_z.y * g_tmp_z.x - g_new_z.x * g_tmp_z.y);
    }
    return 0;
}
