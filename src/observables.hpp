#ifndef _OBSERVABLES_HPP
#define _OBSERVABLES_HPP

#ifndef EXTERN_OBSERVABLES
 #define EXTERN_OBSERVABLES extern
#endif

#include <fstream>
#include <map>

#include "conf.hpp"
#include "theory.hpp"

using namespace std;

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
  
  void ave_err(double &ave,double &err)
  {
    ave=p.first/n;
    err=p.second/n;
    err-=ave*ave;
    err=sqrt(err/(n-1));
  };
  
  string ave_err_str()
  {
    ostringstream os;
    os.precision(16);
    double ave,err;
    ave_err(ave,err);
    os<<ave;
    if(n>1) os<<" "<<err;
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
    //open
    ofstream kin_ener_out(path+"kinetic_energy");
    ofstream common_pot_out(path+"common_potential");
    ofstream mass_pot_out(path+"mass_potential");
    ofstream ener_out(path+"energy");
    ofstream constraint_out(path+"constraint");
    ofstream trace_out(path+"trace");
    ofstream sq_X_trace_out(path+"sq_X_trace");
    ofstream sq_Y_trace_out(path+"sq_Y_trace");
    ofstream sq_Y_trace_ch1_out(path+"sq_Y_trace_ch1");
    ofstream sq_Y_trace_ch2_out(path+"sq_Y_trace_ch2");
    ofstream eig_x0_out(path+"eigenvalues_x0");
    // ofstream eig_x1_out(path+"eigenvalues_x1");
    // ofstream eig_y0_out(path+"eigenvalues_y0");
    
    //check
    if(!kin_ener_out.good()) CRASH("check paths");
    
    //set precision
    kin_ener_out.precision(16);
    common_pot_out.precision(16);
    mass_pot_out.precision(16);
    ener_out.precision(16);
    constraint_out.precision(16);
    sq_X_trace_out.precision(16);
    sq_Y_trace_out.precision(16);
    sq_Y_trace_ch1_out.precision(16);
    sq_Y_trace_ch2_out.precision(16);
    trace_out.precision(16);
    eig_x0_out.precision(16);
    // eig_x1_out.precision(16);
    // eig_y0_out.precision(16);
    
    for(auto &x : kin_ener) kin_ener_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : common_pot) common_pot_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : mass_pot) mass_pot_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : ener) ener_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : constraint) constraint_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_X_trace) sq_X_trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace) sq_Y_trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace_ch1) sq_Y_trace_ch1_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : sq_Y_trace_ch2) sq_Y_trace_ch2_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(auto &x : trace) trace_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    for(int i=0;i<N;i++)
      {
	for(auto &x : eig_x0) eig_x0_out<<x.first*meas_each<<" "<<x.second[i].ave_err_str()<<endl;
	eig_x0_out<<"&"<<endl;
      }
    //for(auto &x : eig_x1) eig_x1_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    //for(auto &x : eig_y0) eig_y0_out<<x.first*meas_each<<" "<<x.second.ave_err_str()<<endl;
    
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
  map<int,obs_t> sq_Y_trace;
  map<int,obs_t> sq_Y_trace_ch1;
  map<int,obs_t> sq_Y_trace_ch2;
  map<int,obs_t> trace;
  map<int,array<obs_t,N> > eig_x0;
  // map<int,array<obs_t,N> > eig_x1;
  // map<int,array<obs_t,N> > eig_y0;
};

#endif
