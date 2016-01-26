#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdarg>

#include "matr.hpp"
#include "observables.hpp"
#include "random.hpp"
#include "tools.hpp"
#include "theory.hpp"
#include "update.hpp"

//length of trajectory
double T,dt;

init_setup_pars_t init_pars;
obs_file obs;
theory_t theory;

//integrate for a given time
void integrate(conf_t &conf,theory_t &theory,double DT,obs_file &obs,double meas_each)
{
  static double t=0;
  static double meas_t=0;
  
  update_t updater;
  updater.dt=dt;
  
  //compute number of steps needed to integrate
  int nt=DT/updater.dt;
  cout<<"Number of integration steps: "<<nt<<" to integrate "<<T<<" in steps of "<<updater.dt<<endl;  
  
  for(int it=0;it<nt;it++)
    {
      cout<<"Integration step "<<it+1<<"/"<<nt<<endl;
      
      //meas
      if(t>=meas_t)
	{
	  obs.measure_all(t,theory,conf);
	  meas_t=t+meas_each;
	}
      
      //update
      t+=updater.update(conf,theory,t);
    }
}

int main()
{
  //open the input
  ifstream input("input");
  if(!input.good()) CRASH("unable to open \"input\"");
  
  read(glb_N,input,"glb_N");
  read(seed,input,"seed");
  read(T,input,"T");
  read(dt,input,"dt");
  read(theory.mass,input,"mass");
  read(init_pars.hx,input,"hx");
  read(init_pars.hy,input,"hy");
  read(init_pars.v,input,"v");
  string init_format_str;
  read(init_format_str,input,"init_format");
  init_pars.kind=init_setup_find_from_string(init_format_str);
  if(init_pars.kind==init_angular) read(init_pars.path,input,"path");
  
  conf_t conf;
  
  //set things that depends on glb_N
  if(glb_N<=3) CRASH("glb_N must be >3");
  fill_generators();
  
  //init random generator device
  gen.seed(seed);
  
  //generate random matrices+L
  conf.generate(init_pars);
  
  //putting to zero the mass
  theory.mass=0;
  
  integrate(conf,theory,T,obs,0.1);
  theory.mass=1e-5;
  integrate(conf,theory,1,obs,0.1);
  theory.mass=0;
  integrate(conf,theory,T,obs,0.1);
  
  return 0;
}
