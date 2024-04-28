/* d_win32_disk.c
 *
 * Routines for a Win32 disk video mode driver for id.
 */
#include "port.h"

#include "calcfrac.h"
#include "cmdfiles.h"
#include "diskvid.h"
#include "drivers.h"
#include "get_color.h"
#include "id_data.h"
#include "os.h"
#include "plot3d.h"
#include "put_color_a.h"
#include "rotate.h"
#include "spindac.h"
#include "video_mode.h"

#define WIN32_LEAN_AND_MEAN
#define STRICT
#include <Windows.h>

#include <cassert>

#include "WinText.h"
#include "frame.h"
#include "d_win32.h"
#include "instance.h"
#include "ods.h"

// read/write-a-dot/line routines
using t_dotwriter = void(int, int, int);
using t_dotreader = int(int, int);
using t_linewriter = void(int y, int x, int lastx, BYTE *pixels);
using t_linereader = void(int y, int x, int lastx, BYTE *pixels);

#define DRAW_INTERVAL 6
#define TIMER_ID 1

class Win32DiskDriver : public Win32BaseDriver
{
public:
    Win32DiskDriver() :
        Win32BaseDriver("disk", "Windows Disk")
    {
    }

    bool init(int *argc, char **argv) override;
    bool resize() override;
    int read_palette() override;
    int write_palette() override;
    void schedule_alarm(int secs) override;
    void write_pixel(int x, int y, int color) override;
    int read_pixel(int x, int y) override;
    void write_span(int y, int x, int lastx, BYTE *pixels) override;
    void read_span(int y, int x, int lastx, BYTE *pixels) override;
    void set_line_mode(int mode) override;
    void draw_line(int x1, int y1, int x2, int y2, int color) override;
    void redraw() override;
    void unget_key(int key) override;
    void window() override;
    void set_video_mode(VIDEOINFO *mode) override;
    void set_clear() override;
    void display_string(int x, int y, int fg, int bg, char const *text) override;
    void hide_text_cursor() override;
    bool is_text() override;
    void set_for_text() override;
    void set_for_graphics() override;
    bool diskp() override;
    bool validate_mode(VIDEOINFO *mode) override;
    void pause() override;
    void resume() override;
    void save_graphics() override;
    void restore_graphics() override;
    void get_max_screen(int &xmax, int &ymax) override;
    void flush() override;

private:
    int width{};
    int height{};
    unsigned char clut[256][3]{};
};

#define DRIVER_MODE(width_, height_ ) \
    { 0, width_, height_, 256, nullptr, "                        " }
static VIDEOINFO modes[] =
{
    DRIVER_MODE( 800,  600),
    DRIVER_MODE(1024,  768),
    DRIVER_MODE(1200,  900),
    DRIVER_MODE(1280,  960),
    DRIVER_MODE(1400, 1050),
    DRIVER_MODE(1500, 1125),
    DRIVER_MODE(1600, 1200)
};
#undef DRIVER_MODE

/*----------------------------------------------------------------------
*
* initdacbox --
*
* Put something nice in the dac.
*
* The conditions are:
*   Colors 1 and 2 should be bright so ifs fractals show up.
*   Color 15 should be bright for lsystem.
*   Color 1 should be bright for bifurcation.
*   Colors 1, 2, 3 should be distinct for periodicity.
*   The color map should look good for mandelbrot.
*
* Results:
*   None.
*
* Side effects:
*   Loads the dac.
*
*----------------------------------------------------------------------
*/
static void initdacbox()
{
    for (int i = 0; i < 256; i++)
    {
        g_dac_box[i][0] = (i >> 5)*8+7;
        g_dac_box[i][1] = (((i+16) & 28) >> 2)*8+7;
        g_dac_box[i][2] = (((i+2) & 3))*16+15;
    }
    g_dac_box[0][2] = 0;
    g_dac_box[0][1] = g_dac_box[0][2];
    g_dac_box[0][0] = g_dac_box[0][1];
    g_dac_box[1][2] = 255;
    g_dac_box[1][1] = g_dac_box[1][2];
    g_dac_box[1][0] = g_dac_box[1][1];
    g_dac_box[2][0] = 190;
    g_dac_box[2][2] = 255;
    g_dac_box[2][1] = g_dac_box[2][2];
}

/***********************************************************************
////////////////////////////////////////////////////////////////////////
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
***********************************************************************/

bool Win32DiskDriver::init(int *argc, char **argv)
{
    const bool base_init = Win32BaseDriver::init(argc, argv);
    if (!base_init)
    {
        return false;
    }

    initdacbox();

    // add default list of video modes
    for (VIDEOINFO &mode : modes)
    {
        add_video_mode(this, &mode);
    }

    return true;
}

/* resize
 *
 * Check if we need resizing.  If no, return 0.
 * If yes, resize internal buffers and return 1.
 */
bool Win32DiskDriver::resize()
{
    frame_resize(wintext.max_width, wintext.max_height);
    if ((g_video_table[g_adapter].xdots == width)
        && (g_video_table[g_adapter].ydots == height))
    {
        return false;
    }

    if (g_disk_flag)
    {
        enddisk();
    }
    startdisk();

    return true;
}


/*----------------------------------------------------------------------
* read_palette
*
*   Reads the current video palette into g_dac_box.
*
*
* Results:
*   None.
*
* Side effects:
*   Fills in g_dac_box.
*
*----------------------------------------------------------------------
*/
int Win32DiskDriver::read_palette()
{
    ODS("disk_read_palette");
    if (!g_got_real_dac)
    {
        return -1;
    }
    for (int i = 0; i < 256; i++)
    {
        g_dac_box[i][0] = clut[i][0];
        g_dac_box[i][2] = clut[i][2];
    }
    return 0;
}

/*
*----------------------------------------------------------------------
*
* write_palette --
*   Writes g_dac_box into the video palette.
*
*
* Results:
*   None.
*
* Side effects:
*   Changes the displayed colors.
*
*----------------------------------------------------------------------
*/
int Win32DiskDriver::write_palette()
{
    ODS("disk_write_palette");
    for (int i = 0; i < 256; i++)
    {
        clut[i][0] = g_dac_box[i][0];
        clut[i][1] = g_dac_box[i][1];
        clut[i][2] = g_dac_box[i][2];
    }

    return 0;
}

/*
*----------------------------------------------------------------------
*
* schedule_alarm --
*
*   Start the refresh alarm
*
* Results:
*   None.
*
* Side effects:
*   Starts the alarm.
*
*----------------------------------------------------------------------
*/
void Win32DiskDriver::schedule_alarm(int secs)
{
    wintext_schedule_alarm(&wintext, (secs ? 1 : DRAW_INTERVAL)*1000);
}

/*
*----------------------------------------------------------------------
*
* write_pixel --
*
*   Write a point to the screen
*
* Results:
*   None.
*
* Side effects:
*   Draws point.
*
*----------------------------------------------------------------------
*/
void Win32DiskDriver::write_pixel(int x, int y, int color)
{
    putcolor_a(x, y, color);
}

/*
*----------------------------------------------------------------------
*
* read_pixel --
*
*   Read a point from the screen
*
* Results:
*   Value of point.
*
* Side effects:
*   None.
*
*----------------------------------------------------------------------
*/
int Win32DiskDriver::read_pixel(int x, int y)
{
    return getcolor(x, y);
}

/*
*----------------------------------------------------------------------
*
* write_span --
*
*   Write a line of pixels to the screen.
*
* Results:
*   None.
*
* Side effects:
*   Draws pixels.
*
*----------------------------------------------------------------------
*/
void Win32DiskDriver::write_span(int y, int x, int lastx, BYTE *pixels)
{
    int width = lastx-x+1;
    ODS3("disk_write_span (%d,%d,%d)", y, x, lastx);

    for (int i = 0; i < width; i++)
    {
        write_pixel(x+i, y, pixels[i]);
    }
}

/*
*----------------------------------------------------------------------
*
* read_span --
*
*   Reads a line of pixels from the screen.
*
* Results:
*   None.
*
* Side effects:
*   Gets pixels
*
*----------------------------------------------------------------------
*/
void Win32DiskDriver::read_span(int y, int x, int lastx, BYTE *pixels)
{
    ODS3("disk_read_span (%d,%d,%d)", y, x, lastx);
    int width = lastx-x+1;
    for (int i = 0; i < width; i++)
    {
        pixels[i] = read_pixel(x+i, y);
    }
}

void Win32DiskDriver::set_line_mode(int mode)
{
    ODS1("disk_set_line_mode %d", mode);
}

void Win32DiskDriver::draw_line(int x1, int y1, int x2, int y2, int color)
{
    ODS5("disk_draw_line (%d,%d) (%d,%d) %d", x1, y1, x2, y2, color);
    ::draw_line(x1, y1, x2, y2, color);
}

/*
*----------------------------------------------------------------------
*
* redraw --
*
*   Refresh the screen.
*
* Results:
*   None.
*
* Side effects:
*   Redraws the screen.
*
*----------------------------------------------------------------------
*/
void Win32DiskDriver::redraw()
{
    ODS("disk_redraw");
    wintext_paintscreen(&wintext, 0, 80, 0, 25);
}

/* unget_key
 *
 * Unread a key!  The key buffer is only one character deep, so we
 * assert if its already full.  This should never happen in real life :-).
 */
void Win32DiskDriver::unget_key(int key)
{
    _ASSERTE(0 == key_buffer);
    key_buffer = key;
}

void Win32DiskDriver::window()
{
    frame_window(wintext.max_width, wintext.max_height);
    wintext.hWndParent = g_frame.window;
    wintext_texton(&wintext);
}

extern void set_disk_dot();
extern void set_normal_line();
void Win32DiskDriver::set_video_mode(VIDEOINFO *mode)
{
    // initially, set the virtual line to be the scan line length
    g_is_true_color = false;            // assume not truecolor
    g_vesa_x_res = 0;                   // reset indicators used for
    g_vesa_y_res = 0;                   // virtual screen limits estimation
    g_good_mode = true;
    if (g_dot_mode != 0)
    {
        g_and_color = g_colors-1;
        g_box_count = 0;
        g_dac_learn = true;
        g_dac_count = g_cycle_limit;
        g_got_real_dac = true;

        read_palette();
    }

    resize();

    set_disk_dot();
    set_normal_line();
}

void Win32DiskDriver::set_clear()
{
    wintext_clear(&wintext);
}

void Win32DiskDriver::display_string(int x, int y, int fg, int bg, char const *text)
{
}

void Win32DiskDriver::hide_text_cursor()
{
    if (cursor_shown)
    {
        cursor_shown = false;
        wintext_hide_cursor(&wintext);
    }
    ODS("disk_hide_text_cursor");
}

bool Win32DiskDriver::is_text()
{
    return true;
}

void Win32DiskDriver::set_for_text()
{
}

void Win32DiskDriver::set_for_graphics()
{
    hide_text_cursor();
}

bool Win32DiskDriver::diskp()
{
    return true;
}

bool Win32DiskDriver::validate_mode(VIDEOINFO *mode)
{
    /* allow modes of any size with 256 colors */
    return mode->colors == 256;
}

void Win32DiskDriver::pause()
{
    if (wintext.hWndCopy)
    {
        ShowWindow(wintext.hWndCopy, SW_HIDE);
    }
}

void Win32DiskDriver::resume()
{
    if (!wintext.hWndCopy)
    {
        window();
    }

    if (wintext.hWndCopy)
    {
        ShowWindow(wintext.hWndCopy, SW_NORMAL);
    }
    wintext_resume(&wintext);
}

void Win32DiskDriver::save_graphics()
{
}

void Win32DiskDriver::restore_graphics()
{
}

void Win32DiskDriver::get_max_screen(int &xmax, int &ymax)
{
    xmax = -1;
    ymax = -1;
}

void Win32DiskDriver::flush()
{
}

Win32DiskDriver disk_driver_info{};

Driver *disk_driver = &disk_driver_info;
