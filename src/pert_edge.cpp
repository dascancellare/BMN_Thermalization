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
  
  double seed;
  double dt;
  
  fill_generators();
  gen.seed(seed);
  
  obs_pars_t obs;
  
  ifstream input(arg[1]);
  double h;
  int niters;
  int iX_pert;
  int non_null_min_mom;
  int non_null_max_mom;
  double eps;
  string base_out;
  double init_random_shift;
  if(!input.good()) CRASH("unable to open \"input\"");
  read(seed,input,"Seed");
  read(init_random_shift,input,"InitRandomShift");
  read(therm_time,input,"ThermTime");
  read(meas_time,input,"MeasTime");
  read(dt,input,"dt");
  read(theory.mass,input,"Mass");
  read(niters,input,"NIters");
  read(h,input,"h");
  read(eps,input,"Eps");
  read(iX_pert,input,"iXPert");
  read(non_null_min_mom,input,"NonNullMinMom");
  read(non_null_max_mom,input,"NonNullMaxMom");
  read(base_out,input,"BaseOut");
  
  int nper_node=max(1,niters/nranks);
  int istart=min(rank*nper_node,niters);
  int iend=min(istart+nper_node,niters);
  if(rank==nranks-1) iend=niters;
  if(nranks>niters) CRASH("Cannot work with %d ranks and %d iters",nranks,niters);
  
  ofstream out_eig_xpert("X0_eigenvalues_prequench");
  for(int iiter=0;iiter<niters;iiter++)
    {
      //generate initial conf
      conf_t conf;
      double shift_time=int(get_real_rand_gauss(square(init_random_shift))/dt)*dt;
      conf.t=shift_time;
      conf.meas_t=-100000*obs.meas_each;
      cout<<conf.t<<endl;
      const size_t ngen=generators.size();
      for(auto &X : conf.X) X.setZero();
      for(int i=0;i<glb_N;i++)
	{
	  conf.P[i].setZero();
	  if(i>=non_null_min_mom && i<non_null_max_mom)
	    for(size_t a=0;a<ngen;a++) conf.P[i]+=generators[a]*get_real_rand_gauss(h/(N*N-1));
	}
      
      if(iiter>=istart &&iiter<iend)
	{
	  ///////////////////   init   //////////////////////
	  cout<<"Rank "<<rank<<", Iter "<<iiter+1<<"/"<<niters<<endl;
	  
	  ////////////////// thermalize /////////////////////
	  
	  //set the evolver and evolve
	  update_t evolver(dt);
	  string path=combine("conf_%d",iiter);
	  if(!file_exists(path))
	    {
	      evolver.integrate(conf,theory,therm_time-shift_time,obs);
	      conf.write(path);
	    }
	  else conf.read(path);
	  
	  //go to the base in which X0 is diagonal
	  SelfAdjointEigenSolver<matr_t> es;
	  conf.gauge_transf(es.compute(conf.X[iX_pert]).eigenvectors());
	  
	  //store the eigenvalues of Y0
	  auto ei=es.compute(conf.X[nX]).eigenvalues();
	  for(int i=0;i<N;i++) out_eig_xpert<<ei(i)<<endl;
	  
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
	      dP(0,i)=-eps*conf.P[iX_pert](0,i)/(eps-conf.X[iX_pert](0,0)+conf.X[iX_pert](i,i));
	      dP(i,N-1)=-eps*conf.P[iX_pert](i,N-1)/(eps-conf.X[iX_pert](i,i)+conf.X[iX_pert](N-1,N-1));
	    }
	  dP(0,N-1)=-2*eps*conf.P[iX_pert](0,N-1)/(2*eps-conf.X[iX_pert](0,0)+conf.X[iX_pert](N-1,N-1));
	  dP=(dP+dP.adjoint()).eval();
	  //make the transformation
	  conf.X[iX_pert]+=dX;
	  conf.P[iX_pert]+=dP;
	  
	  //mark the trace to subtracty
	  sq_X_trace_ref=conf.sq_X_trace();
	  
	  evolver.integrate(conf,theory,meas_time,obs);
	}
    }
  
  obs.write(base_out);
  
  MPI_Finalize();
  
  return 0;
}
