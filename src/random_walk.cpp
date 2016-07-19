#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "conf.hpp"
#include "matr.hpp"
#include "observables.hpp"
#include "random.hpp"
#include "tools.hpp"
#include "theory.hpp"
#include "update.hpp"

#include <mpi.h>

int main(int narg,char **arg)
{
  //init MPI
  int nranks,rank;
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  fill_generators();
  gen.seed(2345);
  
  matr_t M;
  M.setZero();
  const size_t ngen=generators.size();
  
  array<vector<double>,NCOL> o;
  const double h=0.01;

  matr_t id;
  id.setIdentity();
  
  //double next_t=-0;
  for(double t=0;t<1000;t+=h)
    {
      //generate new matr
      //if(t+h/2>=next_t)
      //{
	  //next_t+=10;
      matr_t rand;
      rand.setZero();
      for(size_t a=0;a<ngen;a++) rand+=h*generators[a]*get_real_rand_gauss(1);
      rand+=h*id*get_real_rand_gauss(1);
	  //cout<<"new"<<endl;
	  //}
	  //cout<<rand<<endl;
	
	  M+=rand*h;
	//M=(M+rand*h)/(1+h);
      
      //M/=1+h;
      SelfAdjointEigenSolver<matr_t> es;
      auto e=es.compute(M).eigenvalues();
      for(int i=0;i<NCOL;i++) o[i].push_back(e(i));
    }

  ofstream out("/tmp/eig");
  for(int i=0;i<NCOL;i++)
    {
      double t=0;
      for(auto &x : o[i])
	{
	  out<<t<<" "<<x<<endl;
	  t+=h;
	}
      
      out<<"&"<<endl;
    }
  
  return 0;
}
