// SPDX-License-Identifier: GPL-3.0-only
//
// Thanks to Shirom Makkad fractaltodesktop@gmail.com

#include "engine/perturbation.h"

#include "engine/PertEngine.h"
#include "engine/convert_center_mag.h"
#include "engine/id_data.h"
#include "math/biginit.h"

#include <config/port.h>

#include <stdexcept>
#include <string>

static PertEngine s_pert_engine;

extern  DComplex g_old_z;
extern int g_row;
extern int g_col;
extern long g_color_iter;
extern double g_magnitude_limit;

bool perturbation()
{
    BigStackSaver saved;
    double mandel_width; // width of display
    DComplex center{};
    BFComplex center_bf{};
    double x_mag_factor{};
    double rotation{};
    double skew{};
    LDouble magnification{};
    if (g_bf_math != BFMathType::NONE)
    {
        center_bf.x = alloc_stack(g_bf_length + 2);
        center_bf.y = alloc_stack(g_bf_length + 2);
        cvt_center_mag_bf(center_bf.x, center_bf.y, magnification, x_mag_factor, rotation, skew);
        neg_bf(center_bf.y, center_bf.y);
    }
    else
    {
        LDouble magnification_ld;
        cvt_center_mag(center.x, center.y, magnification_ld, x_mag_factor, rotation, skew);
        center.y = -center.y;
    }

    if (g_bf_math == BFMathType::NONE)
    {
        mandel_width = g_y_max - g_y_min;
    }
    else
    {
        BigFloat tmp_bf{alloc_stack(g_bf_length + 2)};
        sub_bf(tmp_bf, g_bf_y_max, g_bf_y_min);
        mandel_width = bf_to_float(tmp_bf);
    }

    s_pert_engine.initialize_frame(center_bf, {center.x, center.y}, mandel_width / 2.0);
    if (const int result = s_pert_engine.calculate_one_frame(); result < 0)
    {
        throw std::runtime_error("Failed to initialize perturbation engine (" + std::to_string(result) + ")");
    }
    g_calc_status = CalcStatus::COMPLETED;
    return false;
}

int perturbation_per_orbit()
{
//    juliaflag = false;
    std::complex<double> z;
    z = { g_old_z.x, g_old_z.y};

    int status = s_pert_engine.calculate_orbit(g_col, g_row, g_color_iter, &z);
    return status;
}

int perturbation_per_pixel()
{
    int result;
    if (result = s_pert_engine.perturbation_per_pixel(g_col, g_row, g_magnitude_limit); result < 0)
    {
        throw std::runtime_error("Failed to run perturbation pixel (" + std::to_string(result) + ")");
    }
    return result;
}

int perturbation_per_image()
{
    if (const int result = s_pert_engine.calculate_one_frame(); result < 0)
    {
        throw std::runtime_error("Failed to initialize perturbation engine (" + std::to_string(result) + ")");
    }

    s_pert_engine.set_glitch_points_count(0);
    return 0;
}
