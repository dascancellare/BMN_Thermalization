#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "update.hpp"

void update_t::update_positions(conf_t &conf,double step)
{for(int i=0;i<glb_N;i++) conf.X[i]+=conf.P[i]*step;}

void update_t::update_momenta(conf_t &conf,theory_t &theory,double t,double step)
{
  theory.get_force(F,conf,t);
  for(int i=0;i<glb_N;i++) conf.P[i]+=F[i]*step;
}

void update_t::update_leapfrog(conf_t &conf,theory_t &theory,double t)
{
  update_momenta(conf,theory,t,dt/2);
  update_positions(conf,dt);
  update_momenta(conf,theory,t,dt/2);
}

double update_t::update(conf_t &conf, theory_t &theory,double t)
{
  switch(method)
    {
    case LEAPFROG:update_leapfrog(conf,theory,t);break;
    case OMELYAN:CRASH("not implemented yet");
    }
  
  return dt;
}
