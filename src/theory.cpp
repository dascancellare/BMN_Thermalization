#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#include "conf.hpp"
#include "matr.hpp"
#include "theory.hpp"

double theory_t::constraint(conf_t &conf)
{
  matr_t C;
  C.setZero();
  
  for(int i=0;i<glb_N;i++) C+=comm(conf.X[i],conf.P[i]);
  
  return trace_square(C);
}

double theory_t::mass_potential(vector<matr_t> &X,double t)
{
  double V=0;
  
  //potential X
  for(int i=0;i<nX;i++)
    V+=
      mass2*trace_square(X[i])/2+
      mass*(I*X[i]*comm(X[(i+1)%nX],X[(i+2)%nX])).trace().real());
  
  //potential Y
  for(int a=nX;a<glb_N;a++) V+=mass2*trace_square(X[a])/8;
  
  //see ref arXiv:0306054 their -m/3 in (21) is our m in agreement with (1) of 1104.5469
  return V;
}

double theory_t::common_potential(vector<matr_t> &X)
{
  double V=0;
  
  //common part
  for(int a=0;a<glb_N;a++)
    for(int b=0;b<glb_N;b++)
      V+=-trace_square(comm(X[a],X[b]))/2;
  
  if(pert)
    {
      matr_t XA2;
      XA2.Zero();
      for(int a=0;a<glb_N;a++) XA2+=X[a]*X[a];
      
      V+=-c1*XA2.trace().real()-c2*trace_square(XA2);
    }
  
  return V/2;
}

double theory_t::potential(vector<matr_t> &X, double t)
{
  double V=0;
  
  //piece coming from mass
  V+=mass_potential(X,t);
  
  //add the common part (commutator)
  V+=common_potential(X);
  
  return V;
}

double theory_t::hamiltonian(conf_t &conf,double t)
{
  double K=conf.kinetic_energy();
  double V=potential(conf.X,t);
  
  return K+V;
}

void theory_t::get_force(vector<matr_t> &F,conf_t &conf,double t)
{
  // calcola la forza
  for(auto &Fi : F) Fi.setZero();
  
  //primo pezzo diverso per X ed Y
  for(int i=0;i<nX;i++) F[i]=sqm(t)*(-conf.X[i]-3.0*I*comm(conf.X[(i+1)%nX],conf.X[(i+2)%nX]));
  for(int a=nX;a<glb_N;a++) F[a]=-0.25*sqm(t)*conf.X[a];
  
  //secondo pezzo uguale per tutti
  for(int i=0;i<glb_N;i++) for(int a=0;a<glb_N;a++) F[i]+=comm(comm(conf.X[a],conf.X[i]),conf.X[a]);
  
  //perturbation
  if(pert)
    {
      //common part of the anticommutator
      matr_t XA2;
      XA2.Zero();
      for(int a=0;a<glb_N;a++) XA2+=conf.X[a]*conf.X[a];
      
      for(int i=0;i<glb_N;i++) F[i]+=2*c1*conf.X[i]+2*c2*anti_comm(conf.X[i],XA2);
    }
  
  //#define DEBUG
  
#ifdef DEBUG
  int icheck=1;
  double e_bef=hamiltonian();
  
  double eps=1e-6;
  X[icheck](1,0)+=eps/2;
  X[icheck](0,1)+=eps/2;
  
  double e_aft=hamiltonian();
  double f_num=-(e_aft-e_bef)/eps;
  double f_exa=(F[icheck](1,0)+F[icheck](0,1)).real()/2.0;
  cout<<"Exa: "<<f_exa<<", num: "<<f_num<<endl;
  
  X[icheck](1,0)-=eps;
  X[icheck](0,1)-=eps;
#endif
  //print the force
  //double n=0;
  //for(auto &Fi : F) n+=Fi.norm();
  //cout<<n<<endl;
}
