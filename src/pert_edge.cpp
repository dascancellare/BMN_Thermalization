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

init_setup_pars_t init_pars;
theory_t theory;

////////// parameters //////////

double therm_time;
double meas_time;

int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s input",arg[0]);
  //init MPI
  int nranks,rank;
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nranks);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  glb_N=9;
  double seed;
  double dt;
  
  fill_generators();
  gen.seed(seed);
  
  obs_pars_t obs;
  
  double h;
  int niters;
  double eps;
  ifstream input(arg[1]);
  string base_out;
  if(!input.good()) CRASH("unable to open \"input\"");
  read(seed,input,"Seed");
  read(therm_time,input,"ThermTime");
  read(meas_time,input,"MeasTime");
  read(dt,input,"dt");
  read(theory.mass,input,"Mass");
  read(niters,input,"NIters");
  read(h,input,"h");
  read(eps,input,"Eps");
  read(base_out,input,"BaseOut");
  
  int nper_node=niters/nranks;
  int istart=rank*nper_node;
  int iend=istart+nper_node;
  if(rank==nranks-1) iend=niters;
  if(nranks>niters) CRASH("niters %d smaller than number of ranks %d",niters,nranks);
  
  for(int iiter=0;iiter<niters;iiter++)
    {
      //generate initial conf
      conf_t conf;
      const size_t ngen=generators.size();
      for(auto &X : conf.X) X.setZero();
      for(auto &P : conf.P)
	{
	  P.setZero();
	  for(size_t a=0;a<ngen;a++) P+=generators[a]*get_real_rand_gauss(h/(N*N-1));
	}
      
      if(iiter>=istart &&iiter<iend)
	{
	  ///////////////////   init   //////////////////////
	  cout<<"Rank "<<rank<<", Iter "<<iiter+1<<"/"<<niters<<endl;
	  
	  ////////////////// thermalize /////////////////////
	  
	  //set the evolver and evolve
	  update_t evolver(dt);
	  evolver.integrate(conf,theory,therm_time,obs);
	  
	  //go to the base in which X0 is diagonal
	  SelfAdjointEigenSolver<matr_t> es;
	  conf.gauge_transf(es.compute(conf.X[0]).eigenvectors());
	  
	  //shift the largest eigenvalue
	  //define correction for X
	  matr_t dX;
	  dX.setZero();
	  dX(0,0)=-eps;
	  dX(N-1,N-1)=+eps;
	  //define correction for P
	  matr_t dP;
	  dP.setZero();
	  for(int i=1;i<N-1;i++)
	    {
	      dP(0,i)=-eps*conf.P[0](0,i)/(eps-conf.X[0](0,0)+conf.X[0](i,i));
	      dP(i,N-1)=-eps*conf.P[0](i,N-1)/(eps-conf.X[0](i,i)+conf.X[0](N-1,N-1));
	    }
	  dP(0,N-1)=-2*eps*conf.P[0](0,N-1)/(2*eps-conf.X[0](0,0)+conf.X[0](N-1,N-1));
	  dP=(dP+dP.adjoint()).eval();
	  //make the transformation
	  conf.X[0]+=dX;
	  conf.P[0]+=dP;
	  
	  evolver.integrate(conf,theory,meas_time,obs);
	}
    }
  
  obs.write(base_out);
  
  MPI_Finalize();
  
  return 0;
}
