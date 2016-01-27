#ifndef _OBSERVABLES_HPP
#define _OBSERVABLES_HPP

#ifndef EXTERN_OBSERVABLES
 #define EXTERN_OBSERVABLES extern
#endif

#include <fstream>

#include "conf.hpp"
#include "theory.hpp"

using namespace std;

//! keep all files for a given measure set
struct obs_pars_t
{
  double meas_each;  //!< interval between measurement
  
  obs_pars_t(string path="",double meas_each=0.1) : meas_each(meas_each)
  {
    //open
    kin_ener.open(path+"kinetic_energy");
    common_pot.open(path+"common_potential");
    mass_pot.open(path+"mass_potential");
    ener.open(path+"energy");
    constraint.open(path+"constraint");
    trace.open(path+"trace");
    eig_x0.open(path+"eigenvalues_x0");
    eig_x1.open(path+"eigenvalues_x1");
    eig_y0.open(path+"eigenvalues_y0");
    
    //check
    if(!kin_ener.good()) CRASH("check paths");
    
    //set precision
    kin_ener.precision(16);
    common_pot.precision(16);
    mass_pot.precision(16);
    ener.precision(16);
    constraint.precision(16);
    trace.precision(16);
    eig_x0.precision(16);
    eig_x1.precision(16);
    eig_y0.precision(16);
  }
  
  //perform all measurement
  void measure_all(double t,theory_t &theory,conf_t &conf);
  
private:
  ofstream kin_ener;
  ofstream common_pot;
  ofstream mass_pot;
  ofstream ener;
  ofstream constraint;
  ofstream trace;
  ofstream eig_x0;
  ofstream eig_x1;
  ofstream eig_y0;
  
  obs_pars_t(){}
};

#endif
