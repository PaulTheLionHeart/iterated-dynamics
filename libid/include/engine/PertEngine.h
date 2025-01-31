// SPDX-License-Identifier: GPL-3.0-only
//
#pragma once

#include "engine/Point.h"
#include "math/big.h"

#include <complex>
#include <string>
#include <vector>

class PertEngine
{
public:
    void initialize_frame(const BFComplex &center_bf, const std::complex<double> &center, double zoom_radius);
    int calculate_one_frame();
    int perturbation_per_pixel(int x, int y, double bailout);
    int calculate_orbit(int x, int y, long iteration, std::complex<double> *z);
    int calculate_reference();
    bool is_glitched();
    bool is_pixel_complete(int x, int y);
    void set_glitched(bool status);
    void PertEngine::decrement_remaining_point_count();
    void cleanup();


private:
    void reference_zoom_point(const BFComplex &center, int max_iteration);
    void reference_zoom_point(const std::complex<double> &center, int max_iteration);

    std::string m_status;
    std::vector<std::complex<double>> m_xn;
    std::vector<double> m_perturbation_tolerance_check;
    double m_delta_real{};
    double m_delta_imag{};
    std::vector<Point> m_points_remaining;
    long m_remaining_point_count{};
    std::complex<double> m_center{};
    BFComplex m_center_bf{};
    std::complex<double> m_c{};
    BFComplex m_c_bf{};
    double m_zoom_radius{};
    bool m_calculate_glitches{true};
    double m_percent_glitch_tolerance{0.1}; // What percentage of the image is okay to be glitched.
    int m_reference_points{};
    int m_saved_stack{};
    std::vector<int> m_glitches{};
    bool m_glitched{};

    std::complex<double> m_old_reference_coordinate{};

    std::complex<double> m_delta_sub_0{};
    std::complex<double> m_delta_sub_n{};
};
