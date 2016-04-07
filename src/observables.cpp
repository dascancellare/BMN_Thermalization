#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "conf.hpp"
#include "observables.hpp"
#include "theory.hpp"

#include <iostream>

using namespace std;

void obs_pars_t::measure_all(double t,theory_t &theory,conf_t &conf)
{
  int it=(t+meas_each/10)/meas_each;
  
  //CRASH("asking to measure it=%d but only %lu possible",it,kin_ener.size());
  
  kin_ener[it].add(conf.kinetic_energy()-conf.kinetic_energy_trace());
  common_pot[it].add(theory.common_potential(conf.X));
  mass_pot[it].add(theory.mass_potential(conf.X,t));
  ener[it].add(theory.hamiltonian(conf,t));
  constraint[it].add(theory.constraint(conf));
  trace[it].add(conf.X[0].trace().real());
  sq_X_trace[it].add(conf.sq_X_trace());
  sq_Y_trace[it].add(conf.sq_Y_trace());
  sq_Y_trace_ch1[it].add(conf.sq_Y_trace_ch1());
  sq_Y_trace_ch2[it].add(conf.sq_Y_trace_ch2());
  
  auto ei=eigenvalues(conf.X[0]);
  for(int i=0;i<N;i++) eig_x0[it][i].add(ei(i));
  //eig_x0<<t<<" "<<<<endl;
  // eig_x1<<t<<" "<<eigenvalues(conf.X[1])<<endl;
  // eig_y0<<t<<" "<<eigenvalues(conf.X[nX])<<endl;
}
