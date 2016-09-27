#ifndef _UPDATE_HPP
#define _UPDATE_HPP

#include "conf.hpp"
#include "gauge_fix.hpp"
#include "observables.hpp"
#include "theory.hpp"
#include "tools.hpp"

//! evolver
struct update_t
{
  double dt; //!< integration step
  vector<matr_t> F; //!< force
  enum update_method_t{LEAPFROG,OMELYAN}; //!< integration schemes
  update_method_t method;  //!< selected integration scheme
  
  //! constructor
  update_t(double dt) : dt(dt),method(OMELYAN)
  {
    check_gln_N_set();
    F.resize(glb_N);
  }
  
  //! update according to the set method
  double update(conf_t &conf,theory_t &theory,double t);
  
  //! integrate and measure
  void integrate(conf_t &conf,theory_t &theory,double DT,obs_pars_t &obs);
private:
  //! update the positions on the base of momenta
  void update_positions(conf_t &conf,double step);
  
  //! update the momenta
  void update_momenta(conf_t &conf,theory_t &theory,double step,double t);
  
  //! leapfrog integrator
  void update_leapfrog(conf_t &conf,theory_t &theory,double t);
  
  //! Omelyan integrator
  void update_Omelyan(conf_t &conf,theory_t &theory,double t);
  
  update_t(){}
};

#endif
