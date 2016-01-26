#ifndef _CONF_HPP
#define _CONF_HPP

#ifndef EXTERN_CONF
 #define EXTERN_CONF extern
#endif

#include "matr.hpp"
#include "tools.hpp"

///////////////////////////////////// types & globals //////////////////////////////////////

enum init_setup_kind_t{init_static,init_angular};

inline init_setup_kind_t init_setup_find_from_string(string what)
{
  if(what=="static") return init_static;
  if(what!="angular") CRASH("use static or angular");
  return init_angular;
}

struct init_setup_pars_t
{
  init_setup_kind_t kind;
  double v;
  string path;
  double hx,hy;
};

//configuration
struct conf_t
{
  //members
  vector<matr_t> X,P;
  
  //constructors
  conf_t()
  {
    if(glb_N<0) CRASH("please init glb_N before");
    X.resize(glb_N);
    P.resize(glb_N);
  }
  
  //common generation
  void generate(init_setup_pars_t &pars);
  
  //compute the kinetic energy of the trace of the momenta
  double kinetic_energy_trace();
  
  //compute the kinetic energy
  double kinetic_energy();
private:
  //put everything to random and then overwrite with L
  void generate_static(double v);
  
  //solution with angular momenta
  void generate_angular(string &path);
  
};

#undef EXTERN_CONF

#endif
