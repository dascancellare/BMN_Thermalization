#include <complex>
#include <cstdarg>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include "conf.hpp"
#include "matr.hpp"
#include "observables.hpp"
#include "random.hpp"
#include "tools.hpp"
#include "theory.hpp"
#include "update.hpp"

init_setup_pars_t init_pars;
theory_t theory;

struct meas_vect : vector<double>
{
  string print_ave_err()
  {
    double s=0,s2=0;
    for(auto &meas : *this)
      {
	s+=meas;
	s2+=square(meas);
      }
    s/=size();
    s2/=size();
    s2-=s*s;
    
    ostringstream os;
    os<<s<<" "<<sqrt(s2/size());
    
    return os.str();
  }
};

////////// parameters //////////

const double therm_time=100;
const double meas_time=10;
const double pert_time=1;
const double pert_magn=1e-8;
const int niters=2;
const bool gf=true;

int main()
{
  double seed=7859634;
  double dt=0.05;
  theory.mass=1;
  
  fill_generators();
  gen.seed(seed);
  
  init_pars.kind=init_static_traceless;
  init_pars.hx=init_pars.hy=0.01;
  
  map<double,meas_vect> meas_diverge;
  
  //! store the final result
  vector<conf_t> stored_end_perturbed;
  vector<conf_t> stored_end_nonperturbed;
  
  for(int iiter=0;iiter<niters;iiter++)
    {
      ///////////////////   init   //////////////////////
      
      obs_pars_t obs_common;
      obs_pars_t obs_nonperturbed;
      obs_pars_t obs_perturbed;
      
      //generate random matrices+L
      conf_t conf;
      conf.generate(init_pars);
      
      ////////////////// thermalize /////////////////////
      
      //set the evolver and evolve for a common interval
      update_t evolver(dt);
      evolver.integrate(conf,theory,therm_time,obs_common);
      
      ///////////////  perturbate or not ///////////////
      
      //switch on the perturbation
      theory_t pert_theory=theory;
      pert_theory.pert=true;
      pert_theory.c1=get_rand_double(-pert_magn,+pert_magn);
      pert_theory.c2=get_rand_double(-pert_magn,+pert_magn);
      
      //perturbate
      conf_t perturbed_conf=conf;
      gauge_fix_pars_t *perturbed_fixer=NULL;
      if(gf) perturbed_fixer=new gauge_fix_pars_t(perturbed_conf);
      evolver.integrate(perturbed_conf,pert_theory,pert_time,obs_perturbed,perturbed_fixer);
      
      //not perturbated
      gauge_fix_pars_t *fixer=NULL;
      if(gf) fixer=new gauge_fix_pars_t(conf);
      evolver.integrate(conf,theory,pert_time,obs_nonperturbed,fixer);
      
      //////////////////  evolve both  ///////////////////
      
      for(double tdiverg=0;tdiverg<meas_time;tdiverg++)
	{
	  cout<<iiter<<" "<<conf.t<<endl;
	  
	  //get diff
	  double diff=(perturbed_conf-conf).norm();
	  meas_diverge[tdiverg].push_back(diff);
	  
	  evolver.integrate(perturbed_conf,theory,1,obs_perturbed,perturbed_fixer);
	  evolver.integrate(conf,theory,1,obs_nonperturbed,fixer);
	  
	}
      
      //store last conf
      stored_end_nonperturbed.push_back(conf);
      stored_end_perturbed.push_back(perturbed_conf);
      
      obs_common.write("/tmp/common_evolution_");
      obs_nonperturbed.write("/tmp/nonperturbed_");
      obs_perturbed.write("/tmp/perturbed_");
    }
  
  //print the divergence and its error
  ofstream out("/tmp/div");
  for(auto &it : meas_diverge) out<<it.first<<" "<<it.second.print_ave_err()<<endl;
  
  //estimate the average distance of confs at last time
  meas_vect typical_diff;
  for(size_t i=0;i<stored_end_perturbed.size();i++)
    for(size_t j=i+1;i<stored_end_perturbed.size();i++)
      typical_diff.push_back((stored_end_perturbed[i]-stored_end_perturbed[j]).norm());
  
  cout<<typical_diff.print_ave_err()<<endl;
  
  //the final difference cannot be gaugefixed
  
  return 0;
}
