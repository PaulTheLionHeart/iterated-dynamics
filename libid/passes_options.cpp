// SPDX-License-Identifier: GPL-3.0-only
//
#include "passes_options.h"

#include "ValueSaver.h"
#include "cmdfiles.h"
#include "fractals/lorenz.h"
#include "full_screen_prompt.h"
#include "get_corners.h"
#include "helpdefs.h"
#include "id_data.h"
#include "id_keys.h"
#include "sticky_orbits.h"

#include <algorithm>

/*
     passes_options invoked by <p> key
*/

int passes_options()
{
    char const *choices[20];
    char const *pass_calc_modes[] = {"rect", "line"};

    FullScreenValues values[25];
    int i;

    bool const old_keep_screen_coords = g_keep_screen_coords;

    int ret = 0;

pass_option_restart:
    // fill up the choices (and previous values) arrays
    int k = -1;

    choices[++k] = "Periodicity (0=off, <0=show, >0=on, -255..+255)";
    values[k].type = 'i';
    int old_periodicity = g_user_periodicity_value;
    values[k].uval.ival = old_periodicity;

    choices[++k] = "Orbit delay (0 = none)";
    values[k].type = 'i';
    int old_orbit_delay = g_orbit_delay;
    values[k].uval.ival = old_orbit_delay;

    choices[++k] = "Orbit interval (1 ... 255)";
    values[k].type = 'i';
    int old_orbit_interval = (int) g_orbit_interval;
    values[k].uval.ival = old_orbit_interval;

    choices[++k] = "Maintain screen coordinates";
    values[k].type = 'y';
    values[k].uval.ch.val = g_keep_screen_coords ? 1 : 0;

    choices[++k] = "Orbit pass shape (rect, line)";
    //   choices[++k] = "Orbit pass shape (rect,line,func)";
    values[k].type = 'l';
    values[k].uval.ch.vlen = 5;
    values[k].uval.ch.list_len = sizeof(pass_calc_modes)/sizeof(*pass_calc_modes);
    values[k].uval.ch.list = pass_calc_modes;
    values[k].uval.ch.val = (g_draw_mode == 'r') ? 0
                             : (g_draw_mode == 'l') ? 1
                             :   /* function */    2;
    char old_draw_mode = g_draw_mode;

    {
        ValueSaver saved_help_mode{g_help_mode, HelpLabels::HELP_PASSES_OPTIONS};
        i = full_screen_prompt("Passes Options\n"
                              "(not all combinations make sense)\n"
                              "(Press F2 for corner parameters)\n"
                              "(Press F6 for calculation parameters)", k+1, choices, values, 64 | 4, nullptr);
    }
    if (i < 0)
    {
        return -1;
    }

    // now check out the results (*hopefully* in the same order <grin>)
    k = -1;
    int j = 0;   // return code

    g_user_periodicity_value = values[++k].uval.ival;
    g_user_periodicity_value = std::min(g_user_periodicity_value, 255);
    g_user_periodicity_value = std::max(g_user_periodicity_value, -255);
    if (g_user_periodicity_value != old_periodicity)
    {
        j = 1;
    }

    g_orbit_delay = values[++k].uval.ival;
    if (g_orbit_delay != old_orbit_delay)
    {
        j = 1;
    }

    g_orbit_interval = values[++k].uval.ival;
    g_orbit_interval = std::min(g_orbit_interval, 255L);
    g_orbit_interval = std::max(g_orbit_interval, 1L);
    if (g_orbit_interval != old_orbit_interval)
    {
        j = 1;
    }

    g_keep_screen_coords = values[++k].uval.ch.val != 0;
    if (g_keep_screen_coords != old_keep_screen_coords)
    {
        j = 1;
    }
    if (!g_keep_screen_coords)
    {
        g_set_orbit_corners = false;
    }

    {
        switch (values[++k].uval.ch.val)
        {
        default:
        case 0:
            g_draw_mode = 'r';
            break;
        case 1:
            g_draw_mode = 'l';
            break;
        case 2:
            g_draw_mode = 'f';
            break;
        }
    }
    if (g_draw_mode != old_draw_mode)
    {
        j = 1;
    }

    if (i == ID_KEY_F2)
    {
        if (get_screen_corners() > 0)
        {
            ret = 1;
        }
        if (j)
        {
            ret = 1;
        }
        goto pass_option_restart;
    }

    if (i == ID_KEY_F6)
    {
        if (get_corners() > 0)
        {
            ret = 1;
        }
        if (j)
        {
            ret = 1;
        }
        goto pass_option_restart;
    }

    return j + ret;
}
