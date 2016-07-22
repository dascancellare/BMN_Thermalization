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

EXTERN_OBSERVABLES double sq_X_trace_ref INIT_TO(0);

const int nL=(nX*(nX-1))/2+((glb_N-nX)*((glb_N-nX)-1))/2;

struct obs_t
{
  void add(double x)
  {
    p.first+=x;
    p.second+=x*x;
    n++;
  }
  
  obs_t()
  {
    p.first=p.second=0;
    n=0;
  }
  
  int ave_err(double &ave,double &err)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    int ntot;
    
    MPI_Allreduce(&n,&ntot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&p.first,&ave,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&p.second,&err,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    ave/=ntot;
    err/=ntot;
    
    err-=ave*ave;
    err=sqrt(err/(ntot-1));
    
    return ntot;
  };
  
  string ave_err_str()
  {
    ostringstream os;
    os.precision(16);
    double ave,err;
    int ntot=ave_err(ave,err);
    os<<ave;
    if(ntot>1) os<<" "<<err;
    return os.str();
  }
  
private:
  pair<double,double> p;
  int n;
};

//! keep all files for a given measure set
struct obs_pars_t
{
  double meas_each;  //!< interval between measurement
  
  obs_pars_t(double meas_each=0.1) : meas_each(meas_each) {}
  
  void write(string path="")
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    //open
    ofstream kin_ener_out(rank==0?(path+"kinetic_energy"):"/dev/null");
    ofstream common_pot_out(rank==0?(path+"common_potential"):"/dev/null");
    ofstream mass_pot_out(rank==0?(path+"mass_potential"):"/dev/null");
    ofstream ener_out(rank==0?(path+"energy"):"/dev/null");
    ofstream constraint_out(rank==0?(path+"constraint"):"/dev/null");
    ofstream trace_out(rank==0?(path+"trace"):"/dev/null");
    ofstream sq_X_trace_out(rank==0?(path+"sq_X_trace"):"/dev/null");
    ofstream sq_X_trace_sub_out(rank==0?(path+"sq_X_trace_sub"):"/dev/null");
    ofstream sq_Y_trace_out(rank==0?(path+"sq_Y_trace"):"/dev/null");
    ofstream sq_Y_trace_ch1_out(rank==0?(path+"sq_Y_trace_ch1"):"/dev/null");
    ofstream sq_Y_trace_ch2_out(rank==0?(path+"sq_Y_trace_ch2"):"/dev/null");
    ofstream sq_Y_trace_ch_extra_out(rank==0?(path+"sq_Y_trace_ch_extra"):"/dev/null");
    ofstream sq_Y_trace_ch_modulo_out(rank==0?(path+"sq_Y_trace_ch_modulo"):"/dev/null");
    ofstream eig_x0_out(rank==0?(path+"eigenvalues_x0"):"/dev/null");
    ofstream eig_x1_out(rank==0?(path+"eigenvalues_x1"):"/dev/null");
    // ofstream eig_x1_out(path+"eigenvalues_x1");
    // ofstream eig_y0_out(path+"eigenvalues_y0");
    ofstream L_out(rank==0?(path+"L"):"/dev/null");
    
    //check
    if(!kin_ener_out.good()) CRASH("check paths");
    
    //set precision
    kin_ener_out.precision(16);
    common_pot_out.precision(16);
    mass_pot_out.precision(16);
    ener_out.precision(16);
    constraint_out.precision(16);
    sq_X_trace_out.precision(16);
    sq_X_trace_sub_out.precision(16);
    sq_Y_trace_out.precision(16);
    sq_Y_trace_ch1_out.precision(16);
    sq_Y_trace_ch2_out.precision(16);
    sq_Y_trace_ch_extra_out.precision(16);
    trace_out.precision(16);
    eig_x0_out.precision(16);
    eig_x1_out.precision(16);
    // eig_y0_out.precision(16);
    L_out.precision(16);
    
    for(auto &x : kin_ener) kin_ener_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : common_pot) common_pot_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : mass_pot) mass_pot_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : ener) ener_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : constraint) constraint_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_X_trace) sq_X_trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_X_trace_sub) sq_X_trace_sub_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace) sq_Y_trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace_ch1) sq_Y_trace_ch1_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace_ch2) sq_Y_trace_ch2_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace_ch_extra) sq_Y_trace_ch_extra_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : trace) trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(int i=0;i<NCOL;i++)
      {
	for(auto &x : eig_x0) eig_x0_out<<x.first*meas_each<<" "<<x.second[i].ave_err_str()<<endl;
	eig_x0_out<<"&"<<endl;
	for(auto &x : eig_x1) eig_x1_out<<x.first*meas_each<<" "<<x.second[i].ave_err_str()<<endl;
	eig_x1_out<<"&"<<endl;
      }
    //for(auto &x : eig_x1) eig_x1_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    //for(auto &x : eig_y0) eig_y0_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(int i=0;i<nL;i++)
      {
	for(auto &x : L) L_out<<x.first*meas_each<<" "<<x.second[i].ave_err_str()<<endl;
	L_out<<"&"<<endl;
      }
  }
  
  //perform all measurement
  void measure_all(double t,theory_t &theory,conf_t &conf);
  
private:
  map<int,obs_t> kin_ener;
  map<int,obs_t> common_pot;
  map<int,obs_t> mass_pot;
  map<int,obs_t> ener;
  map<int,obs_t> constraint;
  map<int,obs_t> sq_X_trace;
  map<int,obs_t> sq_X_trace_sub;
  map<int,obs_t> sq_Y_trace;
  map<int,obs_t> sq_Y_trace_ch1;
  map<int,obs_t> sq_Y_trace_ch2;
  map<int,obs_t> sq_Y_trace_ch_extra;
  map<int,obs_t> sq_Y_trace_ch_modulo;
  map<int,obs_t> trace;
  map<int,array<obs_t,NCOL> > eig_x0;
  map<int,array<obs_t,NCOL> > eig_x1;
  // map<int,array<obs_t,N> > eig_y0;
  map<int,array<obs_t,nL> > L;
};

#undef EXTERN_OBSERVABLES
#undef INIT_TO

#endif
