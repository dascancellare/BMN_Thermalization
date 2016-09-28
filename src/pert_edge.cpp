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

theory_t theory;

////////// parameters //////////

double dt;
int seed;
int niters;
int iX_pert;
double eps;
double walltime;
string base_out;
double h;
double evol_therm_time;
double meas_time;
int non_null_min_mom;
int non_null_max_mom;

obs_pars_t obs;
int nmulti=1;
int nprocs,proc;
int istart,iend;

//timings used to sed nmulti
double init_time;
double init_meas_time;
double cpu_therm_time;
double avail_meas_time;

void read_input(string path)
{
  ifstream input(path);
  if(!input.good()) CRASH("unable to open \"%s\"",path.c_str());
  read(seed,input,"Seed");
  read(evol_therm_time,input,"ThermTime");
  read(meas_time,input,"MeasTime");
  read(dt,input,"dt");
  read(theory.mass,input,"Mass");
  read(niters,input,"NIters");
  read(walltime,input,"Walltime");
  if(walltime<0) read(nmulti,input,"NMulti");
  read(h,input,"h");
  read(eps,input,"Eps");
  read(iX_pert,input,"iXPert");
  read(non_null_min_mom,input,"NonNullMinMom");
  read(non_null_max_mom,input,"NonNullMaxMom");
  read(base_out,input,"BaseOut");
}

void init(int narg,char **arg)
{
  cout<<"Compiled for NCOL="<<NCOL<<endl;
  
  //init MPI
  MPI_Init(&narg,&arg);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&proc);
  init_time=MPI_Wtime();
  
  fill_generators();
  gen.seed(seed);
  //set the range of configurations relative to the proc
  int nper_node=max(1,niters/nprocs);
  istart=min(proc*nper_node,niters);
  iend=min(istart+nper_node,niters);
  if(proc==nprocs-1) iend=niters;
  if(niters%nprocs) CRASH("Cannot work with %d procs and %d iters",nprocs,niters);
}

void thermalize_or_load(vector<conf_t> &conf,update_t &evolver)
{
  for(int iiter=0;iiter<niters;iiter++)
    {
      //generate initial conf
      conf[iiter].meas_t=-100000*obs.meas_each;
      const size_t ngen=generators.size();
      for(auto &X : conf[iiter].X) X.setZero();
      for(int i=0;i<glb_N;i++)
	{
	  conf[iiter].P[i].setZero();
	  if(i>=non_null_min_mom && i<non_null_max_mom)
	    for(size_t a=0;a<ngen;a++) conf[iiter].P[i]+=generators[a]*get_real_rand_gauss(h/(NCOL*NCOL-1));
	}
      
      if(iiter>=istart &&iiter<iend)
	{
	  ///////////////////   init   //////////////////////
	  cout<<"Proc "<<proc<<", Iter "<<iiter+1<<"/"<<niters<<endl;
	  
	  ////////////////// thermalize /////////////////////
	  
	  //set the evolver and evolve
	  string path=combine("conf_%d",iiter);
	  if(!file_exists(path))
	    {
	      cout<<"File "<<path<<"does not exists, creating it"<<endl;
	      evolver.integrate(conf[iiter],theory,evol_therm_time,obs);
	      conf[iiter].write(path);
	    }
	  else
	    {
	      cout<<"File "<<path<<" exists, loading it"<<endl;
	      conf[iiter].read(path);
	    }
	}
    }
  
  //take note of all time
  init_meas_time=MPI_Wtime();
  cpu_therm_time=init_meas_time-init_time;
  avail_meas_time=walltime-cpu_therm_time;
}

conf_t perturb(conf_t que_conf)
{
  //store the eigenvalues of Y0 before the transformation
  // auto ei=es.compute(que_conf.X[nX]).eigenvalues();
  // for(int i=0;i<NCOL;i++) Y0_pre[i+NCOL*(imulti+nmulti*iiter)]=ei(i);
  
  //go to the base in which X0 is diagonal
  SelfAdjointEigenSolver<matr_t> es;
  que_conf.gauge_transf(es.compute(que_conf.X[iX_pert]).eigenvectors());
  
#if 1
  
  //shift the largest eigenvalue
  //define correction for X
  matr_t dX;
  dX.setZero();
  dX(0,0)=-eps;
  dX(NCOL-1,NCOL-1)=+eps;
  //define correction for P
  matr_t dP;
  dP.setZero();
  for(int i=1;i<NCOL-1;i++)
    {
      dP(0,i)=-eps*que_conf.P[iX_pert](0,i)/(eps-que_conf.X[iX_pert](0,0)+que_conf.X[iX_pert](i,i));
      dP(i,NCOL-1)=-eps*que_conf.P[iX_pert](i,NCOL-1)/(eps-que_conf.X[iX_pert](i,i)+que_conf.X[iX_pert](NCOL-1,NCOL-1));
    }
  dP(0,NCOL-1)=-2*eps*que_conf.P[iX_pert](0,NCOL-1)/(2*eps-que_conf.X[iX_pert](0,0)+que_conf.X[iX_pert](NCOL-1,NCOL-1));
  dP=(dP+dP.adjoint()).eval();
  //make the transformation
  que_conf.X[iX_pert]+=dX;
  que_conf.P[iX_pert]+=dP;
  
#else
  
  que_conf.X[iX_pert]*=(1+eps);
  que_conf.P[iX_pert]*=1.0/(1+eps);
  
#endif
  
  //store the eigenvalues of Y0 at the transformation
  // ei=es.compute(que_conf.X[nX]).eigenvalues();
  // for(int i=0;i<NCOL;i++) Y0_at[i+NCOL*(imulti+nmulti*iiter)]=ei(i);
  
  return que_conf;
}

void update_nmulti(int imulti)
{
  if(walltime>0)
    {
      double curr_meas_time=MPI_Wtime();
      double used_meas_time=curr_meas_time-init_meas_time;
      double time_per_multi=used_meas_time/imulti;
      if(proc==0) cout<<"Time per multi: "<<time_per_multi<<", Time available for meas: "<<avail_meas_time<<", Nmulti "<<nmulti;
      nmulti=avail_meas_time/time_per_multi*0.9;
      MPI_Bcast(&nmulti,1,MPI_INT,0,MPI_COMM_WORLD);
      if(proc==0) cout<<"-> "<<nmulti<<endl;
    }
}

int main(int narg,char **arg)
{
  if(narg<2) CRASH("Use %s input",arg[0]);
  
  read_input(arg[1]);
  init(narg,arg);
  
  // double *Y0_pre=new double[NCOL*niters*nmulti];
  // double *Y0_at=new double[NCOL*niters*nmulti];
  // double *Y0_aft=new double[NCOL*niters*nmulti];
  // memset(Y0_pre,0,sizeof(double)*NCOL*niters*nmulti);
  // memset(Y0_at,0,sizeof(double)*NCOL*niters*nmulti);
  // memset(Y0_aft,0,sizeof(double)*NCOL*niters*nmulti);
  
  //perform thermalization
  vector<conf_t> conf(niters);
  update_t evolver(dt);
  thermalize_or_load(conf,evolver);
  
  for(int imulti=1;imulti<=nmulti;imulti++)
    {
      //write info on where we arrived
      if(proc==0) cout<<"Imulti: "<<imulti<<"/"<<nmulti<<endl;
      
      for(int iiter=istart;iiter<iend;iiter++)
	{
	  //quench the configuration
	  conf_t que_conf=perturb(conf[iiter]);
	  
	  //evolve perturbed
	  evolver.integrate(que_conf,theory,meas_time,obs);
	  
	  //store the eigenvalues of Y0 after retermalization
	  // ei=es.compute(que_conf.X[nX]).eigenvalues();
	  // for(int i=0;i<NCOL;i++) Y0_aft[i+NCOL*(imulti+nmulti*iiter)]=ei(i);
	  
	  //evlove unperturbed
	  double stored_time=conf[iiter].t;
	  double stored_meas_time=conf[iiter].meas_t;
	  evolver.integrate(conf[iiter],theory,meas_time);
	  conf[iiter].t=stored_time;
	  conf[iiter].meas_t=stored_meas_time;
	}
      
      //update the number of multi
      update_nmulti(imulti);
    }
  
  // MPI_Allreduce(MPI_IN_PLACE,Y0_pre,NCOL*niters*nmulti,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // MPI_Allreduce(MPI_IN_PLACE,Y0_at,NCOL*niters*nmulti,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // MPI_Allreduce(MPI_IN_PLACE,Y0_aft,NCOL*niters*nmulti,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // if(proc==0)
  //   {
  //     ofstream out_eig_Y0_pre("Y0_eigenvalues_prequench");
  //     for(int i=0;i<niters*NCOL*nmulti;i++) out_eig_Y0_pre<<Y0_pre[i]<<endl;
  //     ofstream out_eig_Y0_at("Y0_eigenvalues_atquench");
  //     for(int i=0;i<niters*NCOL*nmulti;i++) out_eig_Y0_at<<Y0_at[i]<<endl;
  //     ofstream out_eig_Y0_aft("Y0_eigenvalues_aftquenc");
  //     for(int i=0;i<niters*NCOL*nmulti;i++) out_eig_Y0_aft<<Y0_aft[i]<<endl;
  //   }
  
  obs.write(base_out);
  
  MPI_Finalize();
  
  return 0;
}
