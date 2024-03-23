#pragma once

#include "cmplx.h"

#include <vector>

struct MP;
struct MPC;

extern int                   g_basin;
extern int                   g_bit_shift_less_1;
extern int                   g_c_exponent;
extern double                g_cos_x;
extern int                   g_degree;
extern double                g_degree_minus_1_over_degree;
extern DComplex *            g_float_param;
extern long                  g_fudge_half;
extern long                  g_fudge_one;
extern long                  g_fudge_two;
extern LComplex              g_l_coefficient;
extern LComplex              g_l_init;
extern LComplex              g_l_old_z;
extern LComplex              g_l_new_z;
extern LComplex              g_l_param;
extern LComplex              g_l_param2;
extern LComplex              g_l_temp;
extern long                  g_l_temp_sqr_x;
extern long                  g_l_temp_sqr_y;
extern LComplex *            g_long_param;
extern DComplex              g_marks_coefficient;
extern int                   g_max_color;
extern MP                    g_mp_degree_minus_1_over_degree;
extern MP                    g_mp_one;
extern MP                    g_mp_temp_param2_x;
extern MP                    g_mp_temp2;
extern MP                    g_mp_threshold;
extern std::vector<MPC>      g_mpc_roots;
extern MPC                   g_mpc_temp_param;
extern MP                    g_newton_mp_r_over_d;
extern double                g_newton_r_over_d;
extern bool                  g_overflow;
extern DComplex              g_param_z1;
extern DComplex              g_param_z2;
extern DComplex              g_power_z;
extern double                g_quaternion_c;
extern double                g_quaternion_ci;
extern double                g_quaternion_cj;
extern double                g_quaternino_ck;
extern std::vector<DComplex> g_roots;
extern double                g_sin_x;
extern double                g_temp_sqr_x;
extern double                g_temp_sqr_y;
extern double                g_threshold;
extern long                  g_xx_one;

void FloatPreCalcMagnet2();
void cpower(DComplex *, int, DComplex *);
int lcpower(LComplex *, int, LComplex *, int);
int lcomplex_mult(LComplex, LComplex, LComplex *, int);
int MPCNewtonFractal();
int Barnsley1Fractal();
int Barnsley1FPFractal();
int Barnsley2Fractal();
int Barnsley2FPFractal();
int JuliaFractal();
int JuliafpFractal();
int LambdaFPFractal();
int LambdaFractal();
int SierpinskiFractal();
int SierpinskiFPFractal();
int LambdaexponentFractal();
int LongLambdaexponentFractal();
int FloatTrigPlusExponentFractal();
int LongTrigPlusExponentFractal();
int UnityFractal();
int UnityfpFractal();
int Mandel4Fractal();
int Mandel4fpFractal();
int floatZtozPluszpwrFractal();
int longZpowerFractal();
int longCmplxZpowerFractal();
int floatZpowerFractal();
int floatCmplxZpowerFractal();
int Barnsley3Fractal();
int Barnsley3FPFractal();
int TrigPlusZsquaredFractal();
int TrigPlusZsquaredfpFractal();
int Richard8fpFractal();
int Richard8Fractal();
int PopcornFractal();
int LPopcornFractal();
int PopcornFractal_Old();
int LPopcornFractal_Old();
int PopcornFractalFn();
int LPopcornFractalFn();
int SpiderfpFractal();
int SpiderFractal();
int TetratefpFractal();
int ZXTrigPlusZFractal();
int ScottZXTrigPlusZFractal();
int SkinnerZXTrigSubZFractal();
int ZXTrigPlusZfpFractal();
int ScottZXTrigPlusZfpFractal();
int SkinnerZXTrigSubZfpFractal();
int Sqr1overTrigFractal();
int Sqr1overTrigfpFractal();
int TrigPlusTrigFractal();
int TrigPlusTrigfpFractal();
int ScottTrigPlusTrigFractal();
int ScottTrigPlusTrigfpFractal();
int SkinnerTrigSubTrigFractal();
int SkinnerTrigSubTrigfpFractal();
int TrigXTrigfpFractal();
int TrigXTrigFractal();
int TrigPlusSqrFractal();
int TrigPlusSqrfpFractal();
int ScottTrigPlusSqrFractal();
int ScottTrigPlusSqrfpFractal();
int SkinnerTrigSubSqrFractal();
int SkinnerTrigSubSqrfpFractal();
int TrigZsqrdfpFractal();
int TrigZsqrdFractal();
int SqrTrigFractal();
int SqrTrigfpFractal();
int Magnet1Fractal();
int Magnet2Fractal();
int LambdaTrigFractal();
int LambdaTrigfpFractal();
int LambdaTrigFractal1();
int LambdaTrigfpFractal1();
int LambdaTrigFractal2();
int LambdaTrigfpFractal2();
int ManOWarFractal();
int ManOWarfpFractal();
int CirclefpFractal();
int long_julia_per_pixel();
int long_richard8_per_pixel();
int long_mandel_per_pixel();
int julia_per_pixel();
int mandel_per_pixel();
int mandelfp_per_pixel();
int juliafp_per_pixel();
int MPCjulia_per_pixel();
int otherrichard8fp_per_pixel();
int othermandelfp_per_pixel();
int otherjuliafp_per_pixel();
int LambdaTrigOrTrigFractal();
int LambdaTrigOrTrigfpFractal();
int JuliaTrigOrTrigFractal();
int JuliaTrigOrTrigfpFractal();
int dynamfloat(double *, double *, double*);
int mandelcloudfloat(double *, double *, double*);
int dynam2dfloat();
int fpMODbailout();
int fpREALbailout();
int fpIMAGbailout();
int fpORbailout();
int fpANDbailout();
int fpMANHbailout();
int fpMANRbailout();
int bnMODbailout();
int bnREALbailout();
int bnIMAGbailout();
int bnORbailout();
int bnANDbailout();
int bnMANHbailout();
int bnMANRbailout();
int bfMODbailout();
int bfREALbailout();
int bfIMAGbailout();
int bfORbailout();
int bfANDbailout();
int bfMANHbailout();
int bfMANRbailout();
void set_pixel_calc_functions();
