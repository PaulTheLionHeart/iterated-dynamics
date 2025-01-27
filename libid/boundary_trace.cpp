#include "boundary_trace.h"

#include "calcfrac.h"
#include "cmdfiles.h"
#include "get_color.h"
#include "id_data.h"
#include "stop_msg.h"
#include "work_list.h"

#include <cstring>

enum class direction
{
    North,
    East,
    South,
    West
};
inline int operator+(direction value)
{
    return static_cast<int>(value);
}

static direction going_to;
static int trail_row = 0;
static int trail_col = 0;
static BYTE dstack[4096] = { 0 };              // common temp, two put_line calls

// boundary trace method
constexpr int BK_COLOR{0};

inline direction advance(direction dir, int increment)
{
    return static_cast<direction>((+dir + increment) & 0x03);
}

// take one step in the direction of going_to
static void step_col_row()
{
    switch (going_to)
    {
    case direction::North:
        g_col = trail_col;
        g_row = trail_row - 1;
        break;
    case direction::East:
        g_col = trail_col + 1;
        g_row = trail_row;
        break;
    case direction::South:
        g_col = trail_col;
        g_row = trail_row + 1;
        break;
    case direction::West:
        g_col = trail_col - 1;
        g_row = trail_row;
        break;
    }
}

int boundary_trace()
{
    if (g_inside_color == COLOR_BLACK || g_outside_color == COLOR_BLACK)
    {
        stopmsg("Boundary tracing cannot be used with inside=0 or outside=0");
        return -1;
    }
    if (g_colors < 16)
    {
        stopmsg("Boundary tracing cannot be used with < 16 colors");
        return -1;
    }

    int last_fillcolor_used = -1;
    g_got_status = 2;
    int max_putline_length = 0; // reset max_putline_length
    for (int cur_row = g_i_y_start; cur_row <= g_i_y_stop; cur_row++)
    {
        g_reset_periodicity = true; // reset for a new row
        g_color = BK_COLOR;
        for (int cur_col = g_i_x_start; cur_col <= g_i_x_stop; cur_col++)
        {
            if (getcolor(cur_col, cur_row) != BK_COLOR)
            {
                continue;
            }

            int trail_color = g_color;
            g_row = cur_row;
            g_col = cur_col;
            if ((*g_calc_type)() == -1) // g_color, g_row, g_col are global
            {
                if (g_show_dot != BK_COLOR)   // remove show dot pixel
                {
                    (*g_plot)(g_col, g_row, BK_COLOR);
                }
                if (g_i_y_stop != g_yy_stop)
                {
                    g_i_y_stop = g_yy_stop - (cur_row - g_yy_start); // allow for sym
                }
                add_worklist(g_xx_start, g_xx_stop, cur_col, cur_row, g_i_y_stop, cur_row, 0, g_work_symmetry);
                return -1;
            }
            g_reset_periodicity = false; // normal periodicity checking

            // This next line may cause a few more pixels to be calculated,
            // but at the savings of quite a bit of overhead
            if (g_color != trail_color)
            {
                continue;
            }

            // sweep clockwise to trace outline
            trail_row = cur_row;
            trail_col = cur_col;
            trail_color = g_color;
            const int fill_color_used = g_fill_color > 0 ? g_fill_color : trail_color;
            direction coming_from = direction::West;
            going_to = direction::East;
            unsigned int matches_found = 0;
            bool continue_loop = true;
            do
            {
                step_col_row();
                if (g_row >= cur_row        //
                    && g_col >= g_i_x_start //
                    && g_col <= g_i_x_stop  //
                    && g_row <= g_i_y_stop)
                {
                    g_color = getcolor(g_col, g_row);
                    // g_color, g_row, g_col are global for (*g_calc_type)()
                    if (g_color == BK_COLOR && (*g_calc_type)()== -1)
                    {
                        if (g_show_dot != BK_COLOR)   // remove show dot pixel
                        {
                            (*g_plot)(g_col, g_row, BK_COLOR);
                        }
                        if (g_i_y_stop != g_yy_stop)
                        {
                            g_i_y_stop = g_yy_stop - (cur_row - g_yy_start); // allow for sym
                        }
                        add_worklist(g_xx_start, g_xx_stop, cur_col, cur_row, g_i_y_stop, cur_row, 0, g_work_symmetry);
                        return -1;
                    }
                    if (g_color == trail_color)
                    {
                        if (matches_found < 4)   // to keep it from overflowing
                        {
                            matches_found++;
                        }
                        trail_row = g_row;
                        trail_col = g_col;
                        going_to = advance(going_to, -1);
                        coming_from = advance(going_to, -1);
                    }
                    else
                    {
                        going_to = advance(going_to, 1);
                        continue_loop = going_to != coming_from || matches_found > 0;
                    }
                }
                else
                {
                    going_to = advance(going_to, 1);
                    continue_loop = going_to != coming_from || matches_found > 0;
                }
            } while (continue_loop && (g_col != cur_col || g_row != cur_row));

            if (matches_found <= 3)
            {
                // no hole
                g_color = BK_COLOR;
                g_reset_periodicity = true;
                continue;
            }

            // Fill in region by looping around again, filling lines to the left
            // whenever going_to is South or West
            trail_row = cur_row;
            trail_col = cur_col;
            coming_from = direction::West;
            going_to = direction::East;
            do
            {
                bool match_found = false;
                do
                {
                    step_col_row();
                    if (g_row >= cur_row                          //
                        && g_col >= g_i_x_start                   //
                        && g_col <= g_i_x_stop                    //
                        && g_row <= g_i_y_stop                    //
                        && getcolor(g_col, g_row) == trail_color) // getcolor() must be last
                    {
                        if (going_to == direction::South ||
                            (going_to == direction::West && coming_from != direction::East))
                        {
                            // fill a row, but only once
                            int right = g_col;
                            while (--right >= g_i_x_start && (g_color = getcolor(right,g_row)) == trail_color)
                            {
                                // do nothing
                            }
                            if (g_color == BK_COLOR) // check last color
                            {
                                int left = right;
                                while (getcolor(--left,g_row) == BK_COLOR)
                                {
                                    // Should NOT be possible for left < g_i_x_start
                                    // do nothing
                                }
                                left++; // one pixel too far
                                if (right == left)   // only one hole
                                {
                                    (*g_plot)(left,g_row,fill_color_used);
                                }
                                else
                                {
                                    // fill the line to the left
                                    const int length = right - left + 1;
                                    if (fill_color_used != last_fillcolor_used || length > max_putline_length)
                                    {
                                        // only reset dstack if necessary
                                        std::memset(dstack, fill_color_used, length);
                                        last_fillcolor_used = fill_color_used;
                                        max_putline_length = length;
                                    }
                                    sym_fill_line(g_row, left, right, dstack);
                                }
                            } // end of fill line
                        }
                        trail_row = g_row;
                        trail_col = g_col;
                        going_to = advance(going_to, -1);
                        coming_from = advance(going_to, -1);
                        match_found = true;
                    }
                    else
                    {
                        going_to = advance(going_to, 1);
                    }
                } while (!match_found && going_to != coming_from);

                if (!match_found)
                {
                    // next one has to be a match
                    step_col_row();
                    trail_row = g_row;
                    trail_col = g_col;
                    going_to = advance(going_to, -1);
                    coming_from = advance(going_to, -1);
                }
            } while (trail_col != cur_col || trail_row != cur_row);
            g_reset_periodicity = true; // reset after a trace/fill
            g_color = BK_COLOR;
        }
    }
    return 0;
}
