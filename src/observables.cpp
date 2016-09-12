#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_OBSERVABLES
 #include "observables.hpp"

#include "conf.hpp"
#include "theory.hpp"

#include <iostream>

using namespace std;

void obs_pars_t::measure_all(double t,theory_t &theory,conf_t &conf)
{
  int it=(t+meas_each/10)/meas_each;
  
  //CRASH("asking to measure it=%d but only %lu possible",it,kin_ener.size());
  
  // kin_ener[it].add(conf.kinetic_energy()-conf.kinetic_energy_trace());
  // common_pot[it].add(theory.common_potential(conf.X));
  // mass_pot[it].add(theory.mass_potential(conf.X,t));
  ener[it].add(theory.hamiltonian(conf,t));
  // constraint[it].add(theory.constraint(conf));
  // trace[it].add(conf.X[0].trace().real());
  // double temp=conf.sq_X_trace();
  // sq_X_trace[it].add(temp);
  // sq_X_trace_sub[it].add(temp-sq_X_trace_ref);
  // sq_Y_trace[it].add(conf.sq_Y_trace());
  // sq_Ymom_trace[it].add(conf.sq_Ymom_trace());
  // sq_Y_trace_ch1[it].add(conf.sq_Y_trace_ch1());
  sq_Y_trace_ch2[it].add(conf.sq_Y_trace_ch2());
  sq_Y_trace_ch_extra[it].add(conf.sq_Y_trace_ch_extra());
  sq_Y_trace_ch_modulo[it].add(conf.sq_Y_trace_ch_modulo());
  // sq_Ymom_trace_ch_modulo[it].add(conf.sq_Ymom_trace_ch_modulo());
  
  // auto ei0=eigenvalues(conf.X[0]);
  // for(int i=0;i<NCOL;i++) eig_x0[it][i].add(ei0(i));
  // auto ei1=eigenvalues(conf.X[1]);
  // for(int i=0;i<NCOL;i++) eig_x1[it][i].add(ei1(i));
  // for(int i=0;i<NCOL;i++) eig_x01[it][i].add(ei0(i)*ei1(i));
  
  // //angular momentum
  // int nset=2,imin[2]={0,nX},imax[2]={nX,glb_N},ipair=0;
  // for(int iset=0;iset<nset;iset++)
  //   for(int i=imin[iset];i<imax[iset];i++)
  //     for(int j=i+1;j<imax[iset];j++)
  // 	L[it][ipair++].add(pow((conf.X[i]*conf.P[j]-conf.X[j]*conf.P[i]).trace().real(),2));
}
