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
  constraint[it].add(theory.constraint(conf));
  // trace[it].add(conf.X[0].trace().real());
  // double temp=conf.sq_X_trace();
  // sq_X_trace[it].add(temp);
  // sq_X_trace_sub[it].add(temp-sq_X_trace_ref);
  sq_Y_trace[it].add(conf.sq_Y_trace());
  // sq_Ymom_trace[it].add(conf.sq_Ymom_trace());
  // sq_Y_trace_ch1[it].add(conf.sq_Y_trace_ch1());
  sq_Y_trace_ch2[it].add(conf.sq_Y_trace_ch2());
  //sq_Y_trace_ch_extra[it].add(conf.sq_Y_trace_ch_extra());
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

double theory_t::constraint(conf_t &conf)
{
  matr_t C;
  C.setZero();
  
  for(int i=0;i<glb_N;i++) C+=comm(conf.X[i],conf.P[i]);
  
  return trace_square(C);
}

double conf_t::kinetic_energy_trace()
{
  double K=0;
  for(int i=0;i<glb_N;i++) K+=square(P[i].trace().real())/NCOL;
  K/=2;
  
  return K;
}

double conf_t::kinetic_energy()
{
  double K=0;
  for(int i=0;i<glb_N;i++) K+=trace_square(P[i]);
  K/=2;
  
  return K;
}

double conf_t::sq_X_trace()
{
  double S=0;
  for(int i=0;i<nX;i++) S+=trace_square(X[i]);
  
  return S;
}

double conf_t::sq_Y_trace_weighted(double *coef)
{
  double S=0;
  for(int i=nX;i<glb_N;i++) S+=coef[i-nX]*trace_square(X[i]);
  
  return S;
}

double conf_t::sq_Y_trace()
{
  double coef[6]={1,1,1,1,1,1};
  return sq_Y_trace_weighted(coef);
}

double conf_t::sq_Ymom_trace()
{
  double out=0;
  for(int i=nX;i<glb_N;i++) out+=trace_square(P[i]);
  
  return out;
}

double conf_t::sq_Y_trace_ch1()
{
  double coef[6]={1,1,1,1,-2,-2};
  return sq_Y_trace_weighted(coef);
}

double conf_t::sq_Y_trace_ch2()
{
  double coef[6]={1,1,-1,-1,0,0};
  return sq_Y_trace_weighted(coef);
}

double conf_t::sq_Y_trace_ch_extra()
{return (X[3]*X[4]).trace().real();}

double conf_t::sq_Y_trace_ch_modulo()
{return std::norm(((X[3]+I*X[4])*(X[3]+I*X[4])).trace());}

double conf_t::sq_Ymom_trace_ch_modulo()
{return std::norm(((P[3]+I*P[4])*(P[3]+I*P[4])).trace());}
