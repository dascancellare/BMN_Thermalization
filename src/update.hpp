#ifndef _UPDATE_HPP
#define _UPDATE_HPP

#include "conf.hpp"
#include "theory.hpp"
#include "tools.hpp"

struct update_t
{
  double dt;
  vector<matr_t> F;
  enum update_method_t{LEAPFROG,OMELYAN};
  update_method_t method;
  
  //constructor
  update_t()
  {
    dt=0.1;
    method=LEAPFROG;
    if(glb_N<0) CRASH("please init glb_N before");
    F.resize(glb_N);
  }
  
  //update according to the set method
  double update(conf_t &conf,theory_t &theory,double t);
  
private:
  //update the positions on the base of momenta
  void update_positions(conf_t &conf,double step);
  
  //update the momenta
  void update_momenta(conf_t &conf,theory_t &theory,double step,double t);
  
  //full update of positions and momenta
  void update_leapfrog(conf_t &conf,theory_t &theory,double t);
};

#endif
