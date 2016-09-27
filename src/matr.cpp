#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#define EXTERN_MATR
#include "matr.hpp"

using namespace std;

void fill_generators()
{
  for(int i=0;i<NCOL;i++)
    for(int j=i+1;j<NCOL;j++)
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
  for(int i=1;i<NCOL;i++)
    {
      matr_t gen;
      
      gen.Zero();
      double norm=1/sqrt(i*(i+1)/2);
      for(int j=0;j<i;j++) gen(j,j)=norm;
      gen(i,i)=-i*norm;
      generators.push_back(gen);
    }
}
