// SPDX-License-Identifier: GPL-3.0-only
//
#pragma once

bool perturbation();
extern bool mandel_perturbation_setup();
extern bool perturbation_per_image();
extern int perturbation_per_pixel();
extern int perturbation_per_orbit();
extern bool is_pixel_finished(int x, int y);
extern long get_glitch_point_count();
extern int calculate_reference();
extern void cleanup_perturbation();

enum class PerturbationMode
{
    AUTO = 0, // the default
    YES = 1,
    NO = 2
};
extern PerturbationMode g_perturbation;
// Raising this number makes more calculations, but less variation between each calculation (less chance
// of mis-identifying a glitched point).
extern double g_perturbation_tolerance;
