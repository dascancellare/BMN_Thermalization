#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

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

void update_t::update_Omelyan(conf_t &conf,theory_t &theory,double t)
{
  const double lambda=0.1931833;
  update_momenta(conf,theory,t,lambda*dt);
  update_positions(conf,dt/2);
  update_momenta(conf,theory,t,(1-2*lambda)*dt);
  update_positions(conf,dt/2);
  update_momenta(conf,theory,t,lambda*dt);
}

double update_t::update(conf_t &conf,theory_t &theory,double t)
{
  switch(method)
    {
    case LEAPFROG:update_leapfrog(conf,theory,t);break;
    case OMELYAN:update_Omelyan(conf,theory,t);break;
    }
  
  return dt;
}

void update_t::integrate(conf_t &conf,theory_t &theory,double DT,obs_pars_t &obs,gauge_fix_pars_t *gauge_fixer)
{
  //compute number of steps needed to integrate
  int nt=DT/dt;
  //cout<<"Number of integration steps: "<<nt<<" to integrate "<<DT<<" in steps of "<<dt<<endl;
  
  for(int it=0;it<nt;it++)
    {
      //meas
      if(conf.t+obs.meas_each/10>=conf.meas_t)
	{
	  obs.measure_all(conf.t,theory,conf);
	  conf.meas_t=conf.t+obs.meas_each;
	}
      
      //update
      double ret=update(conf,theory,conf.t);
      conf.t+=ret;
      
      //subtract the trace
      if(0)
      for(int i=0;i<glb_N;i++)
	{
	  complex<double> tr;
	  tr=conf.X[i].trace()*(1.0/N);
	  for(int ic=0;ic<N;ic++) conf.X[i](ic,ic)-=tr;
	  tr=conf.P[i].trace()*(1.0/N);
	  for(int ic=0;ic<N;ic++) conf.P[i](ic,ic)-=tr;
	}
      
      //fix the gauge if needed
      if(gauge_fixer)
	{
	  double diff=gauge_fixer->fix(conf);
	  gauge_fixer->ref_conf=conf;
	  cerr<<conf.t<<" "<<diff<<" "<<gauge_fixer->get_pars()[0]<<endl;
	}
    }
  
  //last meas
  if(conf.t>=conf.meas_t)
    {
      obs.measure_all(conf.t,theory,conf);
      conf.meas_t=conf.t+obs.meas_each;
    }
}
