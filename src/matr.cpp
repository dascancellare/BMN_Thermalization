#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#define EXTERN_MATR
#include "matr.hpp"

using namespace std;

vector<matr_t> get_generators(int n)
{
  vector<matr_t> out;
  
  //real generators
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      {
	matr_t gen;
	gen.Zero();
	gen(i,j)=gen(j,i)=1;
	out.push_back(gen);
      }
  
  //diagonal generators
  for(int i=1;i<n;i++)
    {
      matr_t gen;
      gen.Zero();
      double norm=1/sqrt(i*(i+1)/2);
      for(int j=0;j<i;j++) gen(j,j)=norm;
      gen(i,i)=-i*norm;
      out.push_back(gen);
    }
  
  //imag generators
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      {
	matr_t gen;
	gen.Zero();
	gen(i,j)=+I;
	gen(j,i)=-I;
	out.push_back(gen);
      }
  
  return out;
}
