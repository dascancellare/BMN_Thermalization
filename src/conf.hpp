#ifndef _CONF_HPP
#define _CONF_HPP

#ifndef EXTERN_CONF
 #define EXTERN_CONF extern
#endif

#include "matr.hpp"
#include "tools.hpp"

///////////////////////////////////// types & globals //////////////////////////////////////

//! init mode
enum init_setup_kind_t{init_static,init_angular};

//! get the init mode from string
inline init_setup_kind_t init_setup_find_from_string(string what)
{
  if(what=="static") return init_static;
  if(what!="angular") CRASH("use static or angular");
  return init_angular;
}

//! contains all parameters to start a simulation
struct init_setup_pars_t
{
  init_setup_kind_t kind;
  double v;
  string path;
  double hx,hy;
  init_setup_pars_t() : kind(init_static),v(0),hx(0.001),hy(0.001){}
};

//! configuration
struct conf_t
{
  vector<matr_t> X; //!< positions
  vector<matr_t> P; //!< momenta
  double t; //!< simulation time
  double meas_t; //!< last time a measure has been done
  
  //! constructors
  conf_t()
  {
    check_gln_N_set();
    X.resize(glb_N);
    P.resize(glb_N);
    t=meas_t=0;
  }
  
  //! make a gauge transformation
  void gauge_transf(matr_t transf);
  
  //! return gauge transformed
  conf_t get_gauge_transformed(matr_t transf);
  
  //! compute the norm with another conf
  double get_norm_with(conf_t oth);
  
  //! common generation
  void generate(init_setup_pars_t &pars);
  
  //! compute the kinetic energy of the trace of the momenta
  double kinetic_energy_trace();
  
  //! compute the kinetic energy
  double kinetic_energy();
private:
  //! put everything to random and then overwrite with L
  void generate_static(double v);
  
  //! solution with angular momenta
  void generate_angular(string &path);
};

#undef EXTERN_CONF

#endif
