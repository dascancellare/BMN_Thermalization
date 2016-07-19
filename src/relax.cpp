#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

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
const double meas_time=100;
const double pert_time=1;
const double pert_magn=1e-8;
const int niters=2;
const bool gf=true;

int main()
{
  double seed=time(0);//7859634;
  double dt=0.05;
  theory.mass=1;
  
  fill_generators();
  gen.seed(seed);
  
  ///////////////////   init   //////////////////////
  
  obs_pars_t obs;
  
  conf_t conf;
  for(auto &Xi : conf.X) Xi.setZero();
  for(auto &Pi : conf.P) Pi.setZero();
  
  for(int i=0;i<nX+1;i++)
    for(int ir=0;ir<NCOL-1;ir++)
      for(int ic=ir+1;ic<NCOL-1;ic++)
	{
	  complex<double> delta_y=get_comp_rand_gauss(0.1);
	  conf.X[i](ir,ic)=delta_y;
	  conf.X[i](ic,ir)=conj(delta_y);
	}
  
  ////////////////// thermalize /////////////////////
  
  //set the evolver and evolve for a common interval
  update_t evolver(dt);
  evolver.integrate(conf,theory,therm_time,obs);
  
  for(int i=0;i<NCOL;i++)
    {
      conf.X[nX](i,NCOL-1)=get_comp_rand_gauss(1e-8);
      conf.X[nX](NCOL-1,i)=conj(conf.X[nX](i,NCOL-1));
    }
  conf.X[nX](NCOL-1,NCOL-1)=0.001;
  
  evolver.integrate(conf,theory,meas_time,obs);
  
  obs.write("/tmp/obs_");
  
  return 0;
}
