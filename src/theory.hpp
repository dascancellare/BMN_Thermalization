#ifndef _THEORY_HPP
#define _THEORY_HPP

#include "conf.hpp"

struct theory_t
{
  double mass; //!< mass parameter
  bool pert; //!< flag for perturbation
  double c1,c2; //!< gaussian coefficient of perturbation
  
  //! compute the force
  void get_force(vector<matr_t> &F,conf_t &conf,double t);
  
  //! compute the constraint (trace of the square)
  double constraint(conf_t &conf);
  
  //! piece of the potential common to x and y
  double common_potential(vector<matr_t> &X);
  
  //! compute the potential due to mass term
  double mass_potential(vector<matr_t> &X,double t);
  
  //! full potential
  double potential(vector<matr_t> &X,double t);
  
  //! compute the total hamiltonian
  double hamiltonian(conf_t &conf,double t);
  
  theory_t() : mass(0),pert(false),c1(1e-8),c2(1e-8) {}
};

#endif
