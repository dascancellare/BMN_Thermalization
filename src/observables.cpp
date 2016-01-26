#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "conf.hpp"
#include "observables.hpp"
#include "theory.hpp"

void obs_file::measure_all(double t,theory_t &theory,conf_t &conf)
{
  kin_ener<<t<<" "<<conf.kinetic_energy()-conf.kinetic_energy_trace()<<endl;
  common_pot<<t<<" "<<theory.common_potential(conf.X)<<endl;
  mass_pot<<t<<" "<<theory.mass_potential(conf.X,t)<<endl;
  ener<<t<<" "<<theory.hamiltonian(conf,t)<<endl;
  constraint<<t<<" "<<theory.constraint(conf)<<endl;
  trace<<t<<" "<<conf.X[0].trace().real()<<endl;
  eig_x0<<t<<" "<<eigenvalues(conf.X[0])<<endl;
  eig_x1<<t<<" "<<eigenvalues(conf.X[1])<<endl;
  eig_y0<<t<<" "<<eigenvalues(conf.X[3])<<endl;
}
