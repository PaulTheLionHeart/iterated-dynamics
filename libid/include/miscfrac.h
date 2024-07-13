#pragma once

int Bifurcation();
int BifurcAddTrigPi();
int BifurcLambda();
int BifurcLambdaTrig();
int BifurcMay();
bool BifurcMaySetup();
int BifurcSetTrigPi();
int BifurcStewartTrig();
int BifurcVerhulstTrig();

int LongBifurcAddTrigPi();
int LongBifurcLambdaTrig();
int LongBifurcMay();
int LongBifurcSetTrigPi();
int LongBifurcStewartTrig();
int LongBifurcVerhulstTrig();

int diffusion();
bool lya_setup();
int lyapunov();
int plasma();
int popcorn();
int test();

// PHD 240702
bool InitPerturbation();
int DoPerturbation();
// PHD 240709
int init_mand_derivatives();
int run_mand_derivatives();
int init_big_mand_derivatives();
int run_big_mand_derivatives();
// PHD 240710
int init_tierazon();
int run_tierazon();
int init_big_tierazon();
int run_big_tierazon();
int init_art_matrix();
int run_art_matrix();
// PHD 240712
int init_big_art_matrix();
int run_big_art_matrix();
// PHD 240713
bool Fourier();
