#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_CONF
#include "conf.hpp"

#include "random.hpp"
#include "tools.hpp"

#include <iostream>

void conf_t::hermitianize()
{
  for(int i=0;i<glb_N;i++)
    {
      X[i]=(X[i]+X[i].adjoint()).eval()/2;
      P[i]=(P[i]+P[i].adjoint()).eval()/2;
    }
}

void conf_t::gauge_transf(matr_t trans)
{
  for(int i=0;i<glb_N;i++)
    {
      X[i]=trans.adjoint()*X[i]*trans;
      P[i]=trans.adjoint()*P[i]*trans;
    }
  hermitianize();
}

conf_t conf_t::get_gauge_transformed(matr_t trans)
{
  conf_t out=(*this);
  out.gauge_transf(trans);
  
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
