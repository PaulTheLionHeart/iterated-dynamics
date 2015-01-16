/*
        encoder.c - GIF Encoder and associated routines
*/
#include <algorithm>
#include <string>

#include <limits.h>
#include <string.h>
#if defined(XFRACT)
#include <unistd.h>
#else
#include <io.h>
#endif

#include "port.h"
#include "prototyp.h"
#include "fractype.h"
#include "drivers.h"

static bool compress(int rowlimit);
static int shftwrite(BYTE const *color, int numcolors);
static int extend_blk_len(int datalen);
static int put_extend_blk(int block_id, int block_len, char const *block_data);
static int store_item_name(char const *name);
static void setup_save_info(FRACTAL_INFO *save_info);

/*
                        Save-To-Disk Routines (GIF)

The following routines perform the GIF encoding when the 's' key is pressed.

The compression logic in this file has been replaced by the classic
UNIX compress code. We have extensively modified the sources to fit
Fractint's needs, but have left the original credits where they
appear. Thanks to the original authors for making available these
classic and reliable sources. Of course, they are not responsible for
all the changes we have made to integrate their sources into Fractint.

MEMORY ALLOCATION

There are two large arrays:

   long htab[HSIZE]              (5003*4 = 20012 bytes)
   unsigned short codetab[HSIZE] (5003*2 = 10006 bytes)

At the moment these arrays reuse extraseg and strlocn, respectively.

*/

static int numsaves = 0;        // For adjusting 'save-to-disk' filenames
static FILE *g_outfile;
static int last_colorbar;
static bool save16bit;
static int outcolor1s, outcolor2s;
static int startbits;

static BYTE paletteBW[] =
{   // B&W palette
    0, 0, 0, 63, 63, 63,
};

#ifndef XFRACT
static BYTE paletteCGA[] =
{   // 4-color (CGA) palette
    0, 0, 0, 21, 63, 63, 63, 21, 63, 63, 63, 63,
};
#endif

static BYTE paletteEGA[] =
{   // 16-color (EGA/CGA) pal
    0, 0, 0, 0, 0, 42, 0, 42, 0, 0, 42, 42,
    42, 0, 0, 42, 0, 42, 42, 21, 0, 42, 42, 42,
    21, 21, 21, 21, 21, 63, 21, 63, 21, 21, 63, 63,
    63, 21, 21, 63, 21, 63, 63, 63, 21, 63, 63, 63,
};

static int gif_savetodisk(char *filename)      // save-to-disk routine
{
    char tmpmsg[41];                 // before openfile in case of overrun
    char openfile[FILE_MAX_PATH], openfiletype[10];
    char tmpfile[FILE_MAX_PATH];
    char const *period;
    bool newfile = false;
    int interrupted;

restart:
    save16bit = disk16bit;
    if (gif87a_flag)             // not storing non-standard fractal info
        save16bit = 0;

    strcpy(openfile, filename);  // decode and open the filename
    strcpy(openfiletype, DEFAULTFRACTALTYPE);    // determine the file extension
    if (save16bit)
        strcpy(openfiletype, ".pot");

    period = has_ext(openfile);
    if (period != nullptr)
    {
        strcpy(openfiletype, period);
        openfile[period - openfile] = 0;
    }
    if (resave_flag != 1)
        updatesavename(filename); // for next time

    strcat(openfile, openfiletype);

    strcpy(tmpfile, openfile);
    if (access(openfile, 0) != 0)// file doesn't exist
        newfile = true;
    else
    {   // file already exists
        if (!fract_overwrite)
        {
            if (resave_flag == 0)
                goto restart;
            if (!started_resaves)
            {   // first save of a savetime set
                updatesavename(filename);
                goto restart;
            }
        }
        if (access(openfile, 2) != 0)
        {
            sprintf(tmpmsg, "Can't write %s", openfile);
            stopmsg(STOPMSG_NONE, tmpmsg);
            return -1;
        }
        newfile = false;
        int i = (int) strlen(tmpfile);
        while (--i >= 0 && tmpfile[i] != SLASHC)
            tmpfile[i] = 0;
        strcat(tmpfile, "fractint.tmp");
    }

    started_resaves = (resave_flag == 1);
    if (resave_flag == 2)        // final save of savetime set?
        resave_flag = 0;

    g_outfile = fopen(tmpfile, "wb");
    if (g_outfile == nullptr)
    {
        sprintf(tmpmsg, "Can't create %s", tmpfile);
        stopmsg(STOPMSG_NONE, tmpmsg);
        return -1;
    }

    if (driver_diskp())
    {   // disk-video
        char buf[61];
        extract_filename(tmpmsg, openfile);

        sprintf(buf, "Saving %s", tmpmsg);
        dvid_status(1, buf);
    }
#ifdef XFRACT
    else
    {
        driver_put_string(3, 0, 0, "Saving to:");
        driver_put_string(4, 0, 0, openfile);
        driver_put_string(5, 0, 0, "               ");
    }
#endif

    busy = true;

    if (debugflag != debug_flags::benchmark_encoder)
        interrupted = encoder() ? 1 : 0;
    else
        interrupted = timer(2, nullptr);     // invoke encoder() via timer

    busy = false;

    fclose(g_outfile);

    if (interrupted)
    {
        char buf[200];
        sprintf(buf, "Save of %s interrupted.\nCancel to ", openfile);
        if (newfile)
            strcat(buf, "delete the file,\ncontinue to keep the partial image.");
        else
            strcat(buf, "retain the original file,\ncontinue to replace original with new partial image.");
        interrupted = 1;
        if (stopmsg(STOPMSG_CANCEL, buf) < 0)
        {
            interrupted = -1;
            unlink(tmpfile);
        }
    }

    if (!newfile && interrupted >= 0)
    {   // replace the real file
        unlink(openfile);         // success assumed since we checked
        rename(tmpfile, openfile);// earlier with access
    }

    if (!driver_diskp())
    {   // supress this on disk-video
        int outcolor1 = outcolor1s;
        int outcolor2 = outcolor2s;
        for (int j = 0; j <= last_colorbar; j++)
        {
            if ((j & 4) == 0)
            {
                if (++outcolor1 >= colors)
                    outcolor1 = 0;
                if (++outcolor2 >= colors)
                    outcolor2 = 0;
            }
            for (int i = 0; 250*i < xdots; i++)
            {   // clear vert status bars
                putcolor(i, j, getcolor(i, j) ^ outcolor1);
                putcolor(xdots - 1 - i, j,
                         getcolor(xdots - 1 - i, j) ^ outcolor2);
            }
        }
    }
    else                         // disk-video
        dvid_status(1, "");

    if (interrupted)
    {
        texttempmsg(" *interrupted* save ");
        if (init_batch >= batch_modes::NORMAL)
            init_batch = batch_modes::BAILOUT_ERROR_NO_SAVE;         // if batch mode, set error level
        return -1;
    }
    if (timedsave == 0)
    {
        driver_buzzer(buzzer_codes::COMPLETE);
        if (init_batch == batch_modes::NONE)
        {
            extract_filename(tmpfile, openfile);
            sprintf(tmpmsg, " File saved as %s ", tmpfile);
            texttempmsg(tmpmsg);
        }
    }
    if (initsavetime < 0)
        goodbye();
    return 0;
}

enum e_save_format
{
    SAVEFORMAT_GIF = 0,
    SAVEFORMAT_PNG,
    SAVEFORMAT_JPEG
};

int savetodisk(char *filename)
{
    e_save_format format = SAVEFORMAT_GIF;

    switch (format)
    {
    case SAVEFORMAT_GIF:
        return gif_savetodisk(filename);

    default:
        return -1;
    }
}

int savetodisk(std::string &filename)
{
    char buff[FILE_MAX_PATH];
    strncpy(buff, filename.c_str(), FILE_MAX_PATH);
    int const result = savetodisk(buff);
    filename = buff;
    return result;
}

bool encoder()
{
    bool interrupted;
    int width, rowlimit;
    BYTE bitsperpixel, x;
    FRACTAL_INFO save_info;

    if (init_batch != batch_modes::NONE)               // flush any impending keystrokes
        while (driver_key_pressed())
            driver_get_key();

    setup_save_info(&save_info);

#ifndef XFRACT
    bitsperpixel = 0;            // calculate bits / pixel
    for (int i = colors; i >= 2; i /= 2)
        bitsperpixel++;

    startbits = bitsperpixel + 1;// start coding with this many bits
    if (colors == 2)
        startbits++;              // B&W Klooge
#else
    if (colors == 2)
    {
        bitsperpixel = 1;
        startbits = 3;
    }
    else
    {
        bitsperpixel = 8;
        startbits = 9;
    }
#endif

    int i = 0;
    if (gif87a_flag)
    {
        if (fwrite("GIF87a", 6, 1, g_outfile) != 1)
            goto oops;             // old GIF Signature
    }
    else
    {
        if (fwrite("GIF89a", 6, 1, g_outfile) != 1)
            goto oops;             // new GIF Signature
    }

    width = xdots;
    rowlimit = ydots;
    if (save16bit)
    {
        /* pot16bit info is stored as: file:    double width rows, right side
         * of row is low 8 bits diskvid: ydots rows of colors followed by ydots
         * rows of low 8 bits decoder: returns (row of color info then row of
         * low 8 bits) * ydots */
        rowlimit <<= 1;
        width <<= 1;
    }
    if (write2(&width, 2, 1, g_outfile) != 1)
        goto oops;                // screen descriptor
    if (write2(&ydots, 2, 1, g_outfile) != 1)
        goto oops;
    x = (BYTE)(128 + ((6 - 1) << 4) + (bitsperpixel - 1));      // color resolution == 6 bits worth
    if (write1(&x, 1, 1, g_outfile) != 1)
        goto oops;
    if (fputc(0, g_outfile) != 0)
        goto oops;                // background color

    // TODO: pixel aspect ratio should be 1:1?
    if (viewwindow                               // less than full screen?
            && (viewxdots == 0 || viewydots == 0))   // and we picked the dots?
        i = (int)(((double) sydots / (double) sxdots) * 64.0 / screenaspect - 14.5);
    else   // must risk loss of precision if numbers low
        i = (int)((((double) ydots / (double) xdots) / finalaspectratio) * 64 - 14.5);
    if (i < 1)
        i = 1;
    if (i > 255)
        i = 255;
    if (gif87a_flag)
        i = 0;                    // for some decoders which can't handle aspect
    if (fputc(i, g_outfile) != i)
        goto oops;                // pixel aspect ratio

#ifndef XFRACT
    if (colors == 256)
    {   // write out the 256-color palette
        if (g_got_real_dac)
        {   // got a DAC - must be a VGA
            if (!shftwrite((BYTE *) g_dac_box, colors))
                goto oops;
#else
    if (colors > 2)
    {
        if (g_got_real_dac || fake_lut)
        {   // got a DAC - must be a VGA
            if (!shftwrite((BYTE *) g_dac_box, 256))
                goto oops;
#endif
        }
        else
        {   // uh oh - better fake it
            for (int i = 0; i < 256; i += 16)
                if (!shftwrite((BYTE *)paletteEGA, 16))
                    goto oops;
        }
    }
    if (colors == 2)
    {   // write out the B&W palette
        if (!shftwrite((BYTE *)paletteBW, colors))
            goto oops;
    }
#ifndef XFRACT
    if (colors == 4)
    {   // write out the CGA palette
        if (!shftwrite((BYTE *)paletteCGA, colors))
            goto oops;
    }
    if (colors == 16)
    {   // Either EGA or VGA
        if (g_got_real_dac)
        {
            if (!shftwrite((BYTE *) g_dac_box, colors))
                goto oops;
        }
        else
        {   // no DAC - must be an EGA
            if (!shftwrite((BYTE *)paletteEGA, colors))
                goto oops;
        }
    }
#endif

    if (fwrite(",", 1, 1, g_outfile) != 1)
        goto oops;                // Image Descriptor
    i = 0;
    if (write2(&i, 2, 1, g_outfile) != 1)
        goto oops;
    if (write2(&i, 2, 1, g_outfile) != 1)
        goto oops;
    if (write2(&width, 2, 1, g_outfile) != 1)
        goto oops;
    if (write2(&ydots, 2, 1, g_outfile) != 1)
        goto oops;
    if (write1(&i, 1, 1, g_outfile) != 1)
        goto oops;

    bitsperpixel = (BYTE)(startbits - 1);

    if (write1(&bitsperpixel, 1, 1, g_outfile) != 1)
        goto oops;

    interrupted = compress(rowlimit);

    if (ferror(g_outfile))
        goto oops;

    if (fputc(0, g_outfile) != 0)
        goto oops;

    if (!gif87a_flag)
    {   // store non-standard fractal info
        // loadfile.c has notes about extension block structure
        if (interrupted)
            save_info.calc_status = static_cast<short>(calc_status_value::PARAMS_CHANGED);     // partial save is not resumable
        save_info.tot_extend_len = 0;
        if (!resume_data.empty() && save_info.calc_status == static_cast<short>(calc_status_value::RESUMABLE))
        {
            // resume info block, 002
            save_info.tot_extend_len += extend_blk_len(resume_len);
            std::copy(&resume_data[0], &resume_data[resume_len], &block[0]);
            if (!put_extend_blk(2, resume_len, (char *)block))
                goto oops;
        }
        // save_info.fractal_type gets modified in setup_save_info() in float only version, so we need to use fractype.
        //    if (save_info.fractal_type == FORMULA || save_info.fractal_type == FFORMULA)
        if (fractype == fractal_type::FORMULA || fractype == fractal_type::FFORMULA)
            save_info.tot_extend_len += store_item_name(FormName.c_str());
        //    if (save_info.fractal_type == LSYSTEM)
        if (fractype == fractal_type::LSYSTEM)
            save_info.tot_extend_len += store_item_name(LName.c_str());
        //    if (save_info.fractal_type == IFS || save_info.fractal_type == IFS3D)
        if (fractype == fractal_type::IFS || fractype == fractal_type::IFS3D)
            save_info.tot_extend_len += store_item_name(IFSName.c_str());
        if (display_3d <= display_3d_modes::NONE && rangeslen)
        {
            // ranges block, 004
            int const num_bytes = rangeslen*2;
            save_info.tot_extend_len += extend_blk_len(num_bytes);
            std::vector<char> buffer;
            for (int range : ranges)
            {
                // ranges are stored as 16-bit ints in little-endian byte order
                buffer.push_back(range & 0xFF);
                buffer.push_back(static_cast<unsigned>(range & 0xFFFF) >> 8);
            }
            if (!put_extend_blk(4, num_bytes, &buffer[0]))
                goto oops;

        }
        // Extended parameters block 005
        if (bf_math != bf_math_type::NONE)
        {
            save_info.tot_extend_len += extend_blk_len(22 * (bflength + 2));
            /* note: this assumes variables allocated in order starting with
             * bfxmin in init_bf2() in BIGNUM.C */
            if (!put_extend_blk(5, 22 * (bflength + 2), (char *) bfxmin))
                goto oops;
        }

        // Extended parameters block 006
        if (evolving & 1)
        {
            EVOLUTION_INFO esave_info;
            if (!have_evolve_info || calc_status == calc_status_value::COMPLETED)
            {
                esave_info.x_parameter_range = evolve_x_parameter_range;
                esave_info.y_parameter_range = evolve_y_parameter_range;
                esave_info.x_parameter_offset = evolve_x_parameter_offset;
                esave_info.y_parameter_offset = evolve_y_parameter_offset;
                esave_info.discrete_x_parameter_offset = (short) evolve_discrete_x_parameter_offset;
                esave_info.discrete_y_paramter_offset = (short) evolve_discrete_y_parameter_offset;
                esave_info.px              = (short)px;
                esave_info.py              = (short)py;
                esave_info.sxoffs          = (short)sxoffs;
                esave_info.syoffs          = (short)syoffs;
                esave_info.xdots           = (short)xdots;
                esave_info.ydots           = (short)ydots;
                esave_info.image_grid_size = (short) evolve_image_grid_size;
                esave_info.evolving        = (short)evolving;
                esave_info.this_generation_random_seed = (unsigned short) evolve_this_generation_random_seed;
                esave_info.max_random_mutation = evolve_max_random_mutation;
                esave_info.ecount          = (short)(evolve_image_grid_size * evolve_image_grid_size);  // flag for done
            }
            else
            {
                // we will need the resuming information
                esave_info.x_parameter_range = evolve_info.x_parameter_range;
                esave_info.y_parameter_range = evolve_info.y_parameter_range;
                esave_info.x_parameter_offset = evolve_info.x_parameter_offset;
                esave_info.y_parameter_offset = evolve_info.y_parameter_offset;
                esave_info.discrete_x_parameter_offset = (short)evolve_info.discrete_x_parameter_offset;
                esave_info.discrete_y_paramter_offset = (short)evolve_info.discrete_y_paramter_offset;
                esave_info.px              = (short)evolve_info.px;
                esave_info.py              = (short)evolve_info.py;
                esave_info.sxoffs          = (short)evolve_info.sxoffs;
                esave_info.syoffs          = (short)evolve_info.syoffs;
                esave_info.xdots           = (short)evolve_info.xdots;
                esave_info.ydots           = (short)evolve_info.ydots;
                esave_info.image_grid_size = (short)evolve_info.image_grid_size;
                esave_info.evolving        = (short)evolve_info.evolving;
                esave_info.this_generation_random_seed = (unsigned short)evolve_info.this_generation_random_seed;
                esave_info.max_random_mutation = evolve_info.max_random_mutation;
                esave_info.ecount          = evolve_info.ecount;
            }
            for (int i = 0; i < NUMGENES; i++)
                esave_info.mutate[i] = (short)gene_bank[i].mutate;

            for (int i = 0; i < sizeof(esave_info.future) / sizeof(short); i++)
                esave_info.future[i] = 0;

            // some XFRACT logic for the doubles needed here
#ifdef XFRACT
            decode_evolver_info(&esave_info, 0);
#endif
            // evolution info block, 006
            save_info.tot_extend_len += extend_blk_len(sizeof(esave_info));
            if (!put_extend_blk(6, sizeof(esave_info), (char *) &esave_info))
                goto oops;
        }

        // Extended parameters block 007
        if (stdcalcmode == 'o')
        {
            ORBITS_INFO osave_info;
            osave_info.oxmin     = oxmin;
            osave_info.oxmax     = oxmax;
            osave_info.oymin     = oymin;
            osave_info.oymax     = oymax;
            osave_info.ox3rd     = ox3rd;
            osave_info.oy3rd     = oy3rd;
            osave_info.keep_scrn_coords = (short) (keep_scrn_coords ? 1 : 0);
            osave_info.drawmode  = drawmode;
            for (int i = 0; i < sizeof(osave_info.future) / sizeof(short); i++)
                osave_info.future[i] = 0;

            // some XFRACT logic for the doubles needed here
#ifdef XFRACT
            decode_orbits_info(&osave_info, 0);
#endif
            // orbits info block, 007
            save_info.tot_extend_len += extend_blk_len(sizeof(osave_info));
            if (!put_extend_blk(7, sizeof(osave_info), (char *) &osave_info))
                goto oops;
        }

        // main and last block, 001
        save_info.tot_extend_len += extend_blk_len(FRACTAL_INFO_SIZE);
#ifdef XFRACT
        decode_fractal_info(&save_info, 0);
#endif
        if (!put_extend_blk(1, FRACTAL_INFO_SIZE, (char *) &save_info))
        {
            goto oops;
        }
    }

    if (fwrite(";", 1, 1, g_outfile) != 1)
        goto oops;                // GIF Terminator

    return interrupted;

oops:
    {
        fflush(g_outfile);
        stopmsg(STOPMSG_NONE, "Error Writing to disk (Disk full?)");
        return true;
    }
}

// TODO: should we be doing this?  We need to store full colors, not the VGA truncated business.
// shift IBM colors to GIF
static int shftwrite(BYTE const *color, int numcolors)
{
    for (int i = 0; i < numcolors; i++)
        for (int j = 0; j < 3; j++)
        {
            BYTE thiscolor = color[3 * i + j];
            thiscolor = (BYTE)(thiscolor << 2);
            thiscolor = (BYTE)(thiscolor + (BYTE)(thiscolor >> 6));
            if (fputc(thiscolor, g_outfile) != (int) thiscolor)
                return (0);
        }
    return (1);
}

static int extend_blk_len(int datalen)
{
    return (datalen + (datalen + 254) / 255 + 15);
    // data   +     1.per.block   + 14 for id + 1 for null at end
}

static int put_extend_blk(int block_id, int block_len, char const *block_data)
{
    int i, j;
    char header[15];
    strcpy(header, "!\377\013fractint");
    sprintf(&header[11], "%03u", block_id);
    if (fwrite(header, 14, 1, g_outfile) != 1)
        return (0);
    i = (block_len + 254) / 255;
    while (--i >= 0)
    {
        block_len -= (j = std::min(block_len, 255));
        if (fputc(j, g_outfile) != j)
            return (0);
        while (--j >= 0)
            fputc(*(block_data++), g_outfile);
    }
    if (fputc(0, g_outfile) != 0)
        return (0);
    return (1);
}

static int store_item_name(char const *nameptr)
{
    formula_info fsave_info;
    for (int i = 0; i < 40; i++)
        fsave_info.form_name[i] = 0;      // initialize string
    strcpy(fsave_info.form_name, nameptr);
    if (fractype == fractal_type::FORMULA || fractype == fractal_type::FFORMULA)
    {
        fsave_info.uses_p1 = (short) (uses_p1 ? 1 : 0);
        fsave_info.uses_p2 = (short) (uses_p2 ? 1 : 0);
        fsave_info.uses_p3 = (short) (uses_p3 ? 1 : 0);
        fsave_info.uses_ismand = (short) (uses_ismand ? 1 : 0);
        fsave_info.ismand = (short) (ismand ? 1 : 0);
        fsave_info.uses_p4 = (short) (uses_p4 ? 1 : 0);
        fsave_info.uses_p5 = (short) (uses_p5 ? 1 : 0);
    }
    else
    {
        fsave_info.uses_p1 = 0;
        fsave_info.uses_p2 = 0;
        fsave_info.uses_p3 = 0;
        fsave_info.uses_ismand = 0;
        fsave_info.ismand = 0;
        fsave_info.uses_p4 = 0;
        fsave_info.uses_p5 = 0;
    }
    for (int i = 0; i < sizeof(fsave_info.future) / sizeof(short); i++)
        fsave_info.future[i] = 0;
    // formula/lsys/ifs info block, 003
    put_extend_blk(3, sizeof(fsave_info), (char *) &fsave_info);
    return (extend_blk_len(sizeof(fsave_info)));
}

static void setup_save_info(FRACTAL_INFO *save_info)
{
    if (fractype != fractal_type::FORMULA && fractype != fractal_type::FFORMULA)
        maxfn = 0;
    // set save parameters in save structure
    strcpy(save_info->info_id, INFO_ID);
    save_info->version = FRACTAL_INFO_VERSION;

    if (maxit <= SHRT_MAX)
        save_info->iterationsold = (short) maxit;
    else
        save_info->iterationsold = (short) SHRT_MAX;

    save_info->fractal_type = (short) fractype;
    save_info->xmin = xxmin;
    save_info->xmax = xxmax;
    save_info->ymin = yymin;
    save_info->ymax = yymax;
    save_info->creal = param[0];
    save_info->cimag = param[1];
    save_info->videomodeax = (short) g_video_entry.videomodeax;
    save_info->videomodebx = (short) g_video_entry.videomodebx;
    save_info->videomodecx = (short) g_video_entry.videomodecx;
    save_info->videomodedx = (short) g_video_entry.videomodedx;
    save_info->dotmode = (short)(g_video_entry.dotmode % 100);
    save_info->xdots = (short) g_video_entry.xdots;
    save_info->ydots = (short) g_video_entry.ydots;
    save_info->colors = (short) g_video_entry.colors;
    save_info->parm3 = 0;        // pre version==7 fields
    save_info->parm4 = 0;
    save_info->dparm3 = param[2];
    save_info->dparm4 = param[3];
    save_info->dparm5 = param[4];
    save_info->dparm6 = param[5];
    save_info->dparm7 = param[6];
    save_info->dparm8 = param[7];
    save_info->dparm9 = param[8];
    save_info->dparm10 = param[9];
    save_info->fillcolor = (short) fillcolor;
    save_info->potential[0] = (float) potparam[0];
    save_info->potential[1] = (float) potparam[1];
    save_info->potential[2] = (float) potparam[2];
    save_info->rflag = (short) (rflag ? 1 : 0);
    save_info->rseed = (short) rseed;
    save_info->inside = (short) inside;
    if (LogFlag <= SHRT_MAX)
        save_info->logmapold = (short) LogFlag;
    else
        save_info->logmapold = (short) SHRT_MAX;
    save_info->invert[0] = (float) inversion[0];
    save_info->invert[1] = (float) inversion[1];
    save_info->invert[2] = (float) inversion[2];
    save_info->decomp[0] = (short) decomp[0];
    save_info->biomorph = (short) usr_biomorph;
    save_info->symmetry = (short) forcesymmetry;
    for (int i = 0; i < 16; i++)
        save_info->init3d[i] = (short) init3d[i];
    save_info->previewfactor = (short) previewfactor;
    save_info->xtrans = (short) xtrans;
    save_info->ytrans = (short) ytrans;
    save_info->red_crop_left = (short) red_crop_left;
    save_info->red_crop_right = (short) red_crop_right;
    save_info->blue_crop_left = (short) blue_crop_left;
    save_info->blue_crop_right = (short) blue_crop_right;
    save_info->red_bright = (short) red_bright;
    save_info->blue_bright = (short) blue_bright;
    save_info->xadjust = (short) xadjust;
    save_info->yadjust = (short) yadjust;
    save_info->eyeseparation = (short) g_eye_separation;
    save_info->glassestype = (short) g_glasses_type;
    save_info->outside = (short) outside;
    save_info->x3rd = xx3rd;
    save_info->y3rd = yy3rd;
    save_info->calc_status = (short) calc_status;
    save_info->stdcalcmode = (char)((three_pass && stdcalcmode == '3') ? 127 : stdcalcmode);
    if (distest <= 32000)
        save_info->distestold = (short) distest;
    else
        save_info->distestold = 32000;
    save_info->floatflag = floatflag ? 1 : 0;
    if (bailout >= 4 && bailout <= 32000)
        save_info->bailoutold = (short) bailout;
    else
        save_info->bailoutold = 0;

    save_info->calctime = calctime;
    save_info->trigndx[0] = static_cast<BYTE>(trigndx[0]);
    save_info->trigndx[1] = static_cast<BYTE>(trigndx[1]);
    save_info->trigndx[2] = static_cast<BYTE>(trigndx[2]);
    save_info->trigndx[3] = static_cast<BYTE>(trigndx[3]);
    save_info->finattract = (short) (finattract ? 1 : 0);
    save_info->initorbit[0] = initorbit.x;
    save_info->initorbit[1] = initorbit.y;
    save_info->useinitorbit = useinitorbit;
    save_info->periodicity = (short) periodicitycheck;
    save_info->pot16bit = (short) (disk16bit ? 1 : 0);
    save_info->faspectratio = finalaspectratio;
    save_info->system = (short) save_system;

    if (check_back())
        save_info->release = (short) std::min(save_release, g_release);
    else
        save_info->release = (short) g_release;

    save_info->flag3d = (short) display_3d;
    save_info->ambient = (short) Ambient;
    save_info->randomize = (short) RANDOMIZE;
    save_info->haze = (short) haze;
    save_info->transparent[0] = (short) transparent[0];
    save_info->transparent[1] = (short) transparent[1];
    save_info->rotate_lo = (short) rotate_lo;
    save_info->rotate_hi = (short) rotate_hi;
    save_info->distestwidth = (short) distestwidth;
    save_info->mxmaxfp = mxmaxfp;
    save_info->mxminfp = mxminfp;
    save_info->mymaxfp = mymaxfp;
    save_info->myminfp = myminfp;
    save_info->zdots = (short) zdots;
    save_info->originfp = originfp;
    save_info->depthfp = depthfp;
    save_info->heightfp = heightfp;
    save_info->widthfp = widthfp;
    save_info->distfp = distfp;
    save_info->eyesfp = eyesfp;
    save_info->orbittype = (short) neworbittype;
    save_info->juli3Dmode = (short) juli3Dmode;
    save_info->maxfn = maxfn;
    save_info->inversejulia = (short)((static_cast<int>(major_method) << 8) + static_cast<int>(minor_method));
    save_info->bailout = bailout;
    save_info->bailoutest = (short) bailoutest;
    save_info->iterations = maxit;
    save_info->bflength = (short) bnlength;
    save_info->bf_math = (short) bf_math;
    save_info->old_demm_colors = (short) (old_demm_colors ? 1 : 0);
    save_info->logmap = LogFlag;
    save_info->distest = distest;
    save_info->dinvert[0] = inversion[0];
    save_info->dinvert[1] = inversion[1];
    save_info->dinvert[2] = inversion[2];
    save_info->logcalc = (short) Log_Fly_Calc;
    save_info->stoppass = (short) stoppass;
    save_info->quick_calc = (short) (quick_calc ? 1 : 0);
    save_info->closeprox = closeprox;
    save_info->nobof = (short) (nobof ? 1 : 0);
    save_info->orbit_interval = orbit_interval;
    save_info->orbit_delay = (short) orbit_delay;
    save_info->math_tol[0] = math_tol[0];
    save_info->math_tol[1] = math_tol[1];
    for (int i = 0; i < sizeof(save_info->future)/sizeof(short); i++)
        save_info->future[i] = 0;

}

/***************************************************************************
 *
 *  GIFENCOD.C       - GIF Image compression routines
 *
 *  Lempel-Ziv compression based on 'compress'.  GIF modifications by
 *  David Rowley (mgardi@watdcsu.waterloo.edu).
 *  Thoroughly massaged by the Stone Soup team for Fractint's purposes.
 *
 ***************************************************************************/

#define BITSF   12
#define HSIZE  5003            // 80% occupancy

/*
 *
 * GIF Image compression - modified 'compress'
 *
 * Based on: compress.c - File compression ala IEEE Computer, June 1984.
 *
 * By Authors:  Spencer W. Thomas       (decvax!harpo!utah-cs!utah-gr!thomas)
 *              Jim McKie               (decvax!mcvax!jim)
 *              Steve Davies            (decvax!vax135!petsd!peora!srd)
 *              Ken Turkowski           (decvax!decwrl!turtlevax!ken)
 *              James A. Woods          (decvax!ihnp4!ames!jaw)
 *              Joe Orost               (decvax!vax135!petsd!joe)
 *
 */

// prototypes

static void output(int code);
static void char_out(int c);
static void flush_char();
static void cl_block();

static int n_bits;                        // number of bits/code
static int maxbits = BITSF;                // user settable max # bits/code
static int maxcode;                  // maximum code, given n_bits
static int maxmaxcode = (int)1 << BITSF; // should NEVER generate this code
# define MAXCODE(n_bits)        (((int) 1 << (n_bits)) - 1)

BYTE block[4096] = { 0 };

static long htab[HSIZE];
static unsigned short codetab[10240] = { 0 };

/*
 * To save much memory, we overlay the table used by compress() with those
 * used by decompress().  The tab_prefix table is the same size and type
 * as the codetab.  The tab_suffix table needs 2**BITSF characters.  We
 * get this from the beginning of htab.  The output stack uses the rest
 * of htab, and contains characters.  There is plenty of room for any
 * possible stack (stack used to be 8000 characters).
 */

#define tab_prefixof(i)   codetab[i]
#define tab_suffixof(i)   ((char_type *)(htab))[i]
#define de_stack          ((char_type *)&tab_suffixof((int)1 << BITSF))

static int free_ent;                  // first unused entry

/*
 * block compression parameters -- after all codes are used up,
 * and compression rate changes, start over.
 */
static bool clear_flg = false;

/*
 * compress stdin to stdout
 *
 * Algorithm:  use open addressing double hashing (no chaining) on the
 * prefix code / next character combination.  We do a variant of Knuth's
 * algorithm D (vol. 3, sec. 6.4) along with G. Knott's relatively-prime
 * secondary probe.  Here, the modular division first probe is gives way
 * to a faster exclusive-or manipulation.  Also do block compression with
 * an adaptive reset, whereby the code table is cleared when the compression
 * ratio decreases, but after the table fills.  The variable-length output
 * codes are re-sized at this point, and a special CLEAR code is generated
 * for the decompressor.  Late addition:  construct the table according to
 * file size for noticeable speed improvement on small files.  Please direct
 * questions about this implementation to ames!jaw.
 */

static int ClearCode;
static int EOFCode;
static int a_count; // Number of characters so far in this 'packet'
static unsigned long cur_accum = 0;
static int  cur_bits = 0;

/*
 * Define the storage for the packet accumulator
 */
static char accum[256];

static bool compress(int rowlimit)
{
    int outcolor1, outcolor2;
    int ent;
    int disp;
    int hsize_reg;
    int hshift;
    int color;
    int in_count = 0;
    bool interrupted = false;
    int tempkey;

    outcolor1 = 0;               // use these colors to show progress
    outcolor2 = 1;               // (this has nothing to do with GIF)

    if (colors > 2)
    {
        outcolor1 = 2;
        outcolor2 = 3;
    }
    if (((++numsaves) & 1) == 0)
    {   // reverse the colors on alt saves
        int i = outcolor1;
        outcolor1 = outcolor2;
        outcolor2 = i;
    }
    outcolor1s = outcolor1;
    outcolor2s = outcolor2;

    // Set up the necessary values
    cur_accum = 0;
    cur_bits = 0;
    clear_flg = false;
    ent = 0;
    n_bits = startbits;
    maxcode = MAXCODE(n_bits);

    ClearCode = (1 << (startbits - 1));
    EOFCode = ClearCode + 1;
    free_ent = ClearCode + 2;

    a_count = 0;
    hshift = 0;
    for (long fcode = (long) HSIZE;  fcode < 65536L; fcode *= 2L)
        hshift++;
    hshift = 8 - hshift;                // set hash code range bound

    memset(htab, 0xff, (unsigned)HSIZE*sizeof(long));
    hsize_reg = HSIZE;

    output((int)ClearCode);

    for (int rownum = 0; rownum < ydots; rownum++)
    {   // scan through the dots
        for (int ydot = rownum; ydot < rowlimit; ydot += ydots)
        {
            for (int xdot = 0; xdot < xdots; xdot++)
            {
                if (save16bit == 0 || ydot < ydots)
                    color = getcolor(xdot, ydot);
                else
                    color = readdisk(xdot + sxoffs, ydot + syoffs);
                if (in_count == 0)
                {
                    in_count = 1;
                    ent = color;
                    continue;
                }
                long fcode = (long)(((long) color << maxbits) + ent);
                int i = (((int)color << hshift) ^ ent);    // xor hashing

                if (htab[i] == fcode)
                {
                    ent = codetab[i];
                    continue;
                }
                else if ((long)htab[i] < 0)      // empty slot
                    goto nomatch;
                disp = hsize_reg - i;           // secondary hash (after G. Knott)
                if (i == 0)
                    disp = 1;
probe:
                if ((i -= disp) < 0)
                    i += hsize_reg;

                if (htab[i] == fcode)
                {
                    ent = codetab[i];
                    continue;
                }
                if ((long)htab[i] > 0)
                    goto probe;
nomatch:
                output((int) ent);
                ent = color;
                if (free_ent < maxmaxcode)
                {
                    // code -> hashtable
                    codetab[i] = (unsigned short)free_ent++;
                    htab[i] = fcode;
                }
                else
                    cl_block();
            } // end for xdot
            if (! driver_diskp()       // supress this on disk-video
                    && ydot == rownum)
            {
                if ((ydot & 4) == 0)
                {
                    if (++outcolor1 >= colors)
                        outcolor1 = 0;
                    if (++outcolor2 >= colors)
                        outcolor2 = 0;
                }
                for (int i = 0; 250*i < xdots; i++)
                {   // display vert status bars
                    // (this is NOT GIF-related)
                    putcolor(i, ydot, getcolor(i, ydot) ^ outcolor1);
                    putcolor(xdots - 1 - i, ydot,
                             getcolor(xdots - 1 - i, ydot) ^ outcolor2);
                }
                last_colorbar = ydot;
            } // end if !driver_diskp()
            tempkey = driver_key_pressed();
            if (tempkey && (tempkey != 's'))  // keyboard hit - bail out
            {
                interrupted = true;
                rownum = ydots;
                break;
            }
            if (tempkey == 's')
                driver_get_key();   // eat the keystroke
        } // end for ydot
    } // end for rownum

    // Put out the final code.
    output((int)ent);
    output((int) EOFCode);
    return interrupted;
}

/*****************************************************************
 * TAG(output)
 *
 * Output the given code.
 * Inputs:
 *      code:   A n_bits-bit integer.  If == -1, then EOF.  This assumes
 *              that n_bits =< (long)wordsize - 1.
 * Outputs:
 *      Outputs code to the file.
 * Assumptions:
 *      Chars are 8 bits long.
 * Algorithm:
 *      Maintain a BITSF character long buffer (so that 8 codes will
 * fit in it exactly).  Use the VAX insv instruction to insert each
 * code in turn.  When the buffer fills up empty it and start over.
 */


static void output(int code)
{
    static unsigned long masks[] =
    {   0x0000, 0x0001, 0x0003, 0x0007, 0x000F,
        0x001F, 0x003F, 0x007F, 0x00FF,
        0x01FF, 0x03FF, 0x07FF, 0x0FFF,
        0x1FFF, 0x3FFF, 0x7FFF, 0xFFFF
    };

    cur_accum &= masks[ cur_bits ];

    if (cur_bits > 0)
        cur_accum |= ((long)code << cur_bits);
    else
        cur_accum = code;

    cur_bits += n_bits;

    while (cur_bits >= 8)
    {
        char_out((unsigned int)(cur_accum & 0xff));
        cur_accum >>= 8;
        cur_bits -= 8;
    }

    /*
     * If the next entry is going to be too big for the code size,
     * then increase it, if possible.
     */
    if (free_ent > maxcode || clear_flg)
    {
        if (clear_flg)
        {
            n_bits = startbits;
            maxcode = MAXCODE(n_bits);
            clear_flg = false;
        }
        else
        {
            n_bits++;
            if (n_bits == maxbits)
                maxcode = maxmaxcode;
            else
                maxcode = MAXCODE(n_bits);
        }
    }

    if (code == EOFCode)
    {
        // At EOF, write the rest of the buffer.
        while (cur_bits > 0)
        {
            char_out((unsigned int)(cur_accum & 0xff));
            cur_accum >>= 8;
            cur_bits -= 8;
        }

        flush_char();

        fflush(g_outfile);
    }
}

/*
 * Clear out the hash table
 */
static void cl_block()             // table clear for block compress
{
    memset(htab, 0xff, (unsigned)HSIZE*sizeof(long));
    free_ent = ClearCode + 2;
    clear_flg = true;
    output((int)ClearCode);
}

/*
 * Add a character to the end of the current packet, and if it is 254
 * characters, flush the packet to disk.
 */
static void char_out(int c)
{
    accum[ a_count++ ] = (char)c;
    if (a_count >= 254)
        flush_char();
}

/*
 * Flush the packet to disk, and reset the accumulator
 */
static void flush_char()
{
    if (a_count > 0)
    {
        fputc(a_count, g_outfile);
        fwrite(accum, 1, a_count, g_outfile);
        a_count = 0;
    }
}
