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
  auto J=generate_L(NCOL-1);
  for(int i=0;i<nX;i++) X[i].block(0,0,NCOL-1,NCOL-1)=J[i];
  
  //original value for bottom-right corner of X[0]
  X[0](NCOL-1,NCOL-1)=X[0](NCOL-2,NCOL-2)-1.0;
  
  //metti a 0 le P tranne p[0]
  for(int i=0;i<glb_N;i++) P[i].setZero();
  P[0](NCOL-1,NCOL-1)=v;
}

void::conf_t::generate_static_traceless(double v)
{
  //first is zero but for L
  X[0].setZero();
  
  //fill with L
  auto J=generate_L(NCOL);
  for(int i=0;i<nX;i++) X[i]=J[i];
  
  //metti a 0 le P tranne p[0]
  for(int i=0;i<glb_N;i++) P[i].setZero();
  P[0](NCOL-1,NCOL-1)=v;
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
  vector<double> ang_A(NCOL);
  for(int i=0;i<NCOL;i++) if(!(in>>ang_A[i])) CRASH("unable to read A[%d]",i);
  
  //read Z
  vector<double> ang_Z(NCOL);
  ang_Z[NCOL-1]=0;
  for(int i=0;i<NCOL-1;i++)
    {
      if(!(in>>ang_Z[i])) CRASH("unable to read Z[%d]",i);
      ang_Z[NCOL-1]-=ang_Z[i];
    }
  
  //fill X[0], X+ and P+
  matr_t XPlus,PPlus;
  for(int ic=0;ic<NCOL;ic++)
    {
      X[0](ic,ic)=ang_Z[ic];
      XPlus(ic,(ic+1)%NCOL)=ang_A[ic]; //there is a mismatch of a factor 2 because of the way X[1] and X[2] are defined in terms of XPlus
      PPlus(ic,(ic+1)%NCOL)=ang_A[ic]*I*ang_J/(4*NCOL*square(ang_A[ic]));
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
      for(int ir=0;ir<NCOL;ir++)
	for(int ic=ir+1;ic<NCOL;ic++)
	  {
	    complex<double> delta_y=get_comp_rand_gauss(h);
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

void conf_t::hermitianize()
{
  for(int i=0;i<glb_N;i++)
    {
      X[i]=(X[i]+X[i].adjoint()).eval()/2;
      P[i]=(P[i]+P[i].adjoint()).eval()/2;
    }
}

void conf_t::gauge_transf(matr_t fix)
{
  for(int i=0;i<glb_N;i++)
    {
      X[i]=fix.adjoint()*X[i]*fix;
      P[i]=fix.adjoint()*P[i]*fix;
    }
  hermitianize();
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

void conf_t::write(string path)
{
  cout<<"Opening "<<path<<" to write"<<endl;
  
  ofstream out(path);
  if(!out.good()) CRASH("opening %s",path.c_str());
  for(int i=0;i<NCOL;i++)
    for(int j=0;j<NCOL;j++)
      {
	for(auto &Xi : this->X) out.write((char*)&Xi(i,j),sizeof(complex<double>));
	for(auto &Pi : this->P) out.write((char*)&Pi(i,j),sizeof(complex<double>));
      }
  //time
  out<<this->t<<endl;
  //generator
  out<<gen<<endl;
  
  //close
  out.close();
}

void conf_t::read(string path)
{
  cout<<"Opening "<<path<<" to read"<<endl;
  
  ifstream out(path);
  if(!out.good()) CRASH("opening %s",path.c_str());
  for(int i=0;i<NCOL;i++)
    for(int j=0;j<NCOL;j++)
      {
	for(auto &Xi : this->X) if(!(out.read((char*)&Xi(i,j),sizeof(complex<double>)))) CRASH("reading X, i=%d j=%d",i,j);
	for(auto &Pi : this->P) if(!(out.read((char*)&Pi(i,j),sizeof(complex<double>)))) CRASH("reading P, i=%d j=%d",i,j);
      }
  //time
  if(!(out>>this->t)) CRASH("reading time");
  //generator
  if(!(out>>gen)) CRASH("Reading random number generator");
  
  //close
  out.close();
}
