#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_CONF
#include "conf.hpp"

#include "random.hpp"
#include "tools.hpp"

#include <iostream>

void::conf_t::generate_static(double v)
{
  //first is zero but for L
  X[0].setZero();
  
  //fill with L
  auto J=generate_L(n);
  for(int i=0;i<nX;i++) X[i].block(0,0,n,n)=J[i];
  
  //original value for bottom-right corner of X[0]
  X[0](N-1,N-1)=X[0](N-2,N-2)-1.0;
  
  //metti a 0 le P tranne p[0]
  for(int i=0;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

void::conf_t::generate_static_traceless(double v)
{
  //first is zero but for L
  X[0].setZero();
  
  //fill with L
  auto J=generate_L(N);
  for(int i=0;i<nX;i++) X[i]=J[i];
  
  //metti a 0 le P tranne p[0]
  for(int i=0;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

void conf_t::generate_angular(string &path)
{
  //put everything to zero
  for(int i=0;i<nX;i++) X[i].setZero();
  
  //open file
  ifstream in(path);
  if(!in.good()) CRASH("unable to open path: %s",path.c_str());
  
  //momentum
  double ang_J;
  if(!(in>>ang_J)) CRASH("unable to read J");
  
  //read A
  vector<double> ang_A(N);
  for(int i=0;i<N;i++) if(!(in>>ang_A[i])) CRASH("unable to read A[%d]",i);
  
  //read Z
  vector<double> ang_Z(N);
  ang_Z[N-1]=0;
  for(int i=0;i<N-1;i++)
    {
      if(!(in>>ang_Z[i])) CRASH("unable to read Z[%d]",i);
      ang_Z[N-1]-=ang_Z[i];
    }
  
  //fill X[0], X+ and P+
  matr_t XPlus,PPlus;
  for(int ic=0;ic<N;ic++)
    {
      X[0](ic,ic)=ang_Z[ic];
      XPlus(ic,(ic+1)%N)=ang_A[ic]; //there is a mismatch of a factor 2 because of the way X[1] and X[2] are defined in terms of XPlus
      PPlus(ic,(ic+1)%N)=ang_A[ic]*I*ang_J/(4*N*square(ang_A[ic]));
    }
  
  //define X[1] and X[2]
  X[1]=    XPlus+XPlus.transpose();
  X[2]=-I*(XPlus-XPlus.transpose());
  
  //define P
  P[0].setZero();
  P[1]=    PPlus+PPlus.adjoint();
  P[2]=-I*(PPlus-PPlus.adjoint());
}

void conf_t::generate(init_setup_pars_t &pars)
{
  for(auto &Xi : X) Xi.setZero();
  
  //fluctuations
  for(int i=0;i<glb_N;i++) //first X is 0 apart from L
    {
      double h;
      if(i<nX) h=pars.hx;
      else     h=pars.hy;
      for(int ir=0;ir<N;ir++)
	for(int ic=ir+1;ic<N;ic++)
	  {
	    complex<double> delta_y=get_rand_gauss(h);
	    X[i](ir,ic)=delta_y;
	    X[i](ic,ir)=conj(delta_y);
	  }
    }
  
  //discriminate
  switch(pars.kind)
    {
    case init_static:generate_static(pars.v);break;
    case init_static_traceless:generate_static_traceless(pars.v);break;
    case init_angular:generate_angular(pars.path);break;
    default: CRASH("unknown init format %d!",pars.kind);break;
    }
}

double conf_t::kinetic_energy_trace()
{
  double K=0;
  for(int i=0;i<glb_N;i++) K+=square(P[i].trace().real())/N;
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

void conf_t::gauge_transf(matr_t fix)
{
  for(int i=0;i<glb_N;i++)
    {
      X[i]=fix.adjoint()*X[i]*fix;
      P[i]=fix.adjoint()*P[i]*fix;
    }
}

conf_t conf_t::get_gauge_transformed(matr_t fix)
{
  conf_t out=(*this);
  out.gauge_transf(fix);
  
  return out;
}

double conf_t::get_norm_with(conf_t oth)
{
  double norm=0;
  for(int i=0;i<glb_N;i++)
    {
      norm+=(X[i].adjoint()*oth.X[i]).trace().real();
      norm+=(P[i].adjoint()*oth.P[i]).trace().real();
    }
  return norm;
  
}

double conf_t::squared_norm()
{
  double out=0;
  for(auto &x : X) out+=x.squaredNorm();
  for(auto &p : P) out+=p.squaredNorm();
  
  return out;
}
