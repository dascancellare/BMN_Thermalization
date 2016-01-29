#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#define EXTERN_MATR
#include "matr.hpp"

using namespace std;

complex<double> I(0.0,1.0);

void fill_generators()
{
  for(int i=0;i<N;i++)
    for(int j=i+1;j<N;j++)
      {
	matr_t gen;
	
	//real generator
	gen.Zero();
	gen(i,j)=gen(j,i)=1;
	generators.push_back(gen);
	
      	//real generator
	gen.Zero();
	gen(i,j)=+I;
	gen(j,i)=-I;
	generators.push_back(gen);
      }
  
  //diagonal generator
  for(int i=1;i<N;i++)
    {
      matr_t gen;
      
      gen.Zero();
      double norm=1/sqrt(i*(i+1)/2);
      for(int j=0;j<i;j++) gen(j,j)=norm;
      gen(i,i)=-i*norm;
      generators.push_back(gen);
    }
}

vector<Matrix<complex<double>,Dynamic,Dynamic> > generate_L(int dimrep)
{
  double j=(dimrep-1.0)/2;
  
  vector<Matrix<complex<double>,Dynamic,Dynamic> > J(3);
  for(auto &Ji : J) Ji.resize(dimrep,dimrep);
  
  Matrix<complex<double>,Dynamic,Dynamic> JPlus(dimrep,dimrep);
  JPlus.setZero();
  
  for(int indn=0;indn+1<dimrep;indn++) JPlus(indn,indn+1)=sqrt((indn+1)*(2*j-indn));
  
  J[1]=  0.5*(JPlus+JPlus.transpose());
  J[2]=-I*0.5*(JPlus-JPlus.transpose());
  
  J[0].setZero();
  for(int indn=0;indn<dimrep;indn++) J[0](indn,indn)=j-indn;
  
  return J;
}
