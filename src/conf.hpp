#ifndef _CONF_HPP
#define _CONF_HPP

#ifndef EXTERN_CONF
 #define EXTERN_CONF extern
#endif

#include <functional>

#include "matr.hpp"
#include "tools.hpp"

///////////////////////////////////// types & globals //////////////////////////////////////

//! init mode
enum init_setup_kind_t{init_static,init_static_traceless,init_angular};

//! get the init mode from string
inline init_setup_kind_t init_setup_find_from_string(string what)
{
  if(what=="static") return init_static;
  if(what=="static_traceless") return init_static_traceless;
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
  
  //! generate and fill
  conf_t(init_setup_pars_t &pars) : conf_t() {generate(pars);}
  
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
  
  //! compute the summ of the trace of the squares of all X
  double sq_X_trace();
  
  //! write to a file
  void write(string path);
  
  //! read to a file
  void read(string path);
  
  //! compute the summ of the trace of the squares of all Y
  double sq_Y_trace_weighted(double *coef);
  double sq_Y_trace();
  double sq_Ymom_trace();
  double sq_Y_trace_ch1();
  double sq_Y_trace_ch2();
  double sq_Y_trace_ch_extra();
  double sq_Y_trace_ch_modulo();
  double sq_Ymom_trace_ch_modulo();
  
  //! difference
  conf_t operator-(conf_t oth)
  {
    conf_t out;
    transform(X.begin(),X.end(),oth.X.begin(),out.X.begin(),minus<matr_t>());
    transform(P.begin(),P.end(),oth.P.begin(),out.P.begin(),minus<matr_t>());
    return out;
  }
  
  //! put to hermitian
  void hermitianize();
  
  //! return the norm (sum of the square norm of X and P)
  double squared_norm();
  
  //! return the sqrt of the sum of the square norm of X and P
  double norm(){return sqrt(squared_norm());}
private:
  //! put everything to random and then overwrite with L (of 1 dim less)
  void generate_static(double v);
  
  //! same but with full L
  void generate_static_traceless(double v);
  
  //! solution with angular momenta
  void generate_angular(string &path);
};

#undef EXTERN_CONF

#endif
