#ifndef _OBSERVABLES_HPP
#define _OBSERVABLES_HPP

#ifndef EXTERN_OBSERVABLES
 #define EXTERN_OBSERVABLES extern
 #define INIT_TO(a)
#else
 #define INIT_TO(a) =a
#endif

#include <fstream>
#include <map>

#include "conf.hpp"
#include "theory.hpp"

#include <mpi.h>

using namespace std;

const int ngrav=nX*(nX+1)/2-1;
const int nL=(nX*(nX-1))/2+((glb_N-nX)*((glb_N-nX)-1))/2;

struct obs_t : pair<int,double>
{
  obs_t(){first=second=0;}
  void add(double x){first++;second+=x;}
  double ave(){return second/first;}
};

struct obs_vec_t : map<double,obs_t>
{
  void print(string path)
  {
    
    //take rank
    int rank,nranks;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    
    //open
    fstream outf(path,(rank==0)?fstream::out:(fstream::app|fstream::out));
    outf.precision(16);
    if(!outf.good()) CRASH("check paths");
    
    //loop print
    for(int irank=0;irank<nranks;irank++)
      {
	MPI_Barrier(MPI_COMM_WORLD);
	if(irank==rank)
	  for(auto &obs : *this)
	    outf<<irank<<" "<<obs.first<<" "<<obs.second.first<<" "<<obs.second.second<<endl;
      }
  }
};

//! keep all files for a given measure set
struct obs_pars_t
{
  double meas_each;  //!< interval between measurement
  
  obs_pars_t(double meas_each) : meas_each(meas_each) {}
  
  void write(string path="")
  {
    ener.print(path+"energy");
    kin_ener.print(path+"kin_energy");
    common_pot.print(path+"common_pot");
    sq_Y_trace.print(path+"sq_Y_trace");
    //fo_Y_trace.print(path+"fo_Y_trace");
    sq_Y_trace_ch2.print(path+"sq_Y_trace_ch2");
    //q_PY_trace_ch2.print(path+"sq_PY_trace_ch2");
    //sq_Y_trace_ch_extra.print(path+"sq_Y_trace_ch_extra");
    sq_Y_trace_ch_modulo.print(path+"sq_Y_trace_ch_modulo");
    constraint.print(path+"constraint");
    //for(int iL=0;iL<nL;iL++) L[iL].print(path+"L_"+to_string(iL));
    for(int i=0;i<ngrav;i++)
      {
	grav[i+0*ngrav].print(path+"grav_x_"+to_string(i));
	grav[i+1*ngrav].print(path+"grav_p_"+to_string(i));
	grav[i+2*ngrav].print(path+"grav_x_tilde_"+to_string(i));
      }
  }
  
  //perform all measurement
  void measure_all(double t,theory_t &theory,conf_t &conf);
  
private:
  obs_vec_t kin_ener;
  obs_vec_t common_pot;
  // obs_vec_t mass_pot;
  obs_vec_t ener;
  obs_vec_t constraint;
  // obs_vec_t sq_X_trace;
  // obs_vec_t sq_X_trace_sub;
  obs_vec_t sq_Y_trace;
  // obs_vec_t fo_Y_trace;
  // obs_vec_t sq_Ymom_trace;
  // obs_vec_t sq_Y_trace_ch1;
  obs_vec_t sq_Y_trace_ch2;
  // obs_vec_t sq_PY_trace_ch2;
  // obs_vec_t sq_Y_trace_ch_extra;
  obs_vec_t sq_Y_trace_ch_modulo;
  array<obs_vec_t,3*ngrav> grav;
  // obs_vec_t sq_Ymom_trace_ch_modulo;
  // obs_vec_t trace;
  // map<int,array<obs_t,NCOL> > eig_x0;
  // map<int,array<obs_t,NCOL> > eig_x1;
  // map<int,array<obs_t,NCOL> > eig_x01;
  // map<int,array<obs_t,N> > eig_y0;
  array<obs_vec_t,nL> L;
};

#undef EXTERN_OBSERVABLES
#undef INIT_TO

#endif
