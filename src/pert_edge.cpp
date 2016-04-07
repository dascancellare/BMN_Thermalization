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

////////// parameters //////////

const double therm_time=100;
const double meas_time=250;

int main()
{
  glb_N=9;
  double seed=7859634;
  double dt=0.05;
  theory.mass=1;
  
  fill_generators();
  gen.seed(seed);
  
  obs_pars_t obs;
  
  const int niters=3;
  for(int iiter=0;iiter<niters;iiter++)
    {
      ///////////////////   init   //////////////////////
      cout<<"Iter "<<iiter+1<<"/"<<niters<<endl;
      
      //generate initial conf
      double h=0.05;
      conf_t conf;
      const size_t ngen=generators.size();
      for(auto &X : conf.X) X.setZero();
      for(auto &P : conf.P)
	{
	  P.setZero();
	  for(size_t a=0;a<ngen;a++) P+=generators[a]*get_real_rand_gauss(h/sqrt(N));
	}
      
      ////////////////// thermalize /////////////////////
      
      //set the evolver and evolve
      update_t evolver(dt);
      evolver.integrate(conf,theory,therm_time,obs);
      
      //go to the base in which X0 is diagonal
      SelfAdjointEigenSolver<matr_t> es;
      conf.gauge_transf(es.compute(conf.X[0]).eigenvectors());
      
      //shift the largest eigenvalue
      double eps=1;
      //define correction for X
      matr_t dX;
      dX.setZero();
      dX(0,0)=-eps;
      dX(N-1,N-1)=+eps;
      //define correction for P
      matr_t dP;
      dP.setZero();
      for(int i=1;i<N-1;i++)
	{
	  dP(0,i)=-eps*conf.P[0](0,i)/(eps-conf.X[0](0,0)+conf.X[0](i,i));
	  dP(i,N-1)=-eps*conf.P[0](i,N-1)/(eps-conf.X[0](i,i)+conf.X[0](N-1,N-1));
	}
      dP(0,N-1)=-2*eps*conf.P[0](0,N-1)/(2*eps-conf.X[0](0,0)+conf.X[0](N-1,N-1));
      dP=(dP+dP.adjoint()).eval();
      //make the transformation
      conf.X[0]+=dX;
      conf.P[0]+=dP;
      
      evolver.integrate(conf,theory,meas_time,obs);
    }
  
  obs.write("/tmp/obs_");
  
  return 0;
}
