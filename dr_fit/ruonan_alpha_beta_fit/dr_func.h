#include "constant.h"
extern "C"{
   double v_compton_diff_cross_driver_dr_(int *, int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

extern "C"{
   double p_lt_tot_driver_(double *, double *,int *);
}

extern "C"{
   double p_ll_tot_driver_(double *, double *,int *);
}

extern "C"{
   double p_tt_driver_(double *,int *);
}

extern "C"{
   double beta_driver_dr_(double *, double *);
}

extern "C"{
   double alpha_driver_dr_(double *, double *);
}

//calc E' and th for the new W
void Calckinvars(double W, double Q2, double ebeam, double &eth, double &ep, double &g_p_lab, double &g_p, double &Q20)
{
    ep=ebeam-(W*W+Q2-Mp*Mp)/(2.*Mp);
    eth=2.*asin(sqrt(Q2/(4.*ebeam*ep)));
    double ge = (W*W+Q2-Mp*Mp)/(2.*Mp);
    g_p_lab =sqrt(Q2+ge*ge);
    g_p =Mp/W*g_p_lab;
    double q0=Mp-sqrt(Mp*Mp+g_p*g_p);
    Q20=-2.*Mp*q0;
}
//calc E' and th for the new W
void Calc_ep_eth(double W, double Q2, double ebeam, double &eth, double &ep)
{
    ep=ebeam-(W*W+Q2-Mp*Mp)/(2.*Mp);
    eth=2.*asin(sqrt(Q2/(4.*ebeam*ep)));
}
