#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdarg>

#include "conf.hpp"
#include "matr.hpp"
#include "observables.hpp"
#include "random.hpp"
#include "tools.hpp"
#include "theory.hpp"
#include "update.hpp"

init_setup_pars_t init_pars;
theory_t theory;

int main()
{
  glb_N=9;
  double seed=7859634;
  double dt=0.005;
  theory.mass=0;
  
  fill_generators();
  gen.seed(seed);
  
  //////////////////////////////////////////////////////
  
  //open the two stream observable files
  obs_pars_t obs_nonfixed("/tmp/nonfixed_");
  obs_pars_t obs_fixed("/tmp/fixed_");
  
  //generate random matrices+L
  conf_t conf;
  conf.generate(init_pars);
  
  //set the evolver and evolve for a common interval
  update_t evolver(dt);
  evolver.integrate(conf,theory,10,obs_nonfixed);
  
  //set fixing to original configuration and evolve more
  gauge_fix_pars_t fixer(conf);
  evolver.integrate(conf,theory,10,obs_nonfixed);
  
  evolver.dt/=1;
  
  //evolve without gaugefixed
  conf_t conf_nonfixed=conf;
  evolver.integrate(conf_nonfixed,theory,10,obs_nonfixed);
  
  //evolve with gaugefixed
  conf_t conf_fixed=conf;
  evolver.integrate(conf_fixed,theory,10,obs_fixed,&fixer);
  
  return 0;
}
