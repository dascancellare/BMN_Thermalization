#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#include "gauge_fix.hpp"

using namespace std;

namespace
{
  conf_t *transforming;
  conf_t *ref_fix;
  
  void dist(int &npar,double *fuf,double &out,double *p,int flag)
  {
    vector<double> w(npar);
    w.assign(p,p+npar);
    out=-transforming->get_gauge_transformed(generate_sun(w)).get_norm_with(*ref_fix);
  }
}

gauge_fix_pars_t::gauge_fix_pars_t(conf_t ref_conf) : ref_conf(ref_conf)
{
  //init minuit
  check_gln_N_set();
  minu=new TMinuit(generators.size());
  minu->SetFCN(dist);
  minu->SetPrintLevel(-1);
  //do not touch
  for(size_t i=0;i<generators.size();i++) minu->DefineParameter(i,"par",0,0.001,0,0);
  
  //minize after setting tolerance
  double tol=1e-16;
  int iflag;
  minu->SetMaxIterations(10000000);
  minu->mnexcm("SET ERR",&tol,1,iflag);
}

vector<double> gauge_fix_pars_t::get_pars()
{
  vector<double> pars(generators.size());
  for(size_t i=0;i<generators.size();i++)
    {
      double dum;
      minu->GetParameter(i,pars[i],dum);
    }
  return pars;
}

double gauge_fix_pars_t::find_gaugefixing(conf_t &conf)
{
  //set the reference
  ref_fix=&ref_conf;
  transforming=&conf;
  
  //get the norm
  double ch2_bef;
  minu->Eval(generators.size(),NULL,ch2_bef,get_pars().data(),0);
  
  //increase the error by 100
  for(size_t i=0;i<generators.size();i++)
    {
      double par,err;
      minu->GetParameter(i,par,err);
      minu->DefineParameter(i,"par",par,err*100,0,0);
    }
  
  //minimize
  minu->Migrad();
  
  //get the maximized norm
  double ch2_aft;
  minu->Eval(generators.size(),NULL,ch2_aft,get_pars().data(),0);
  return ch2_aft/ch2_bef-1;
}

double gauge_fix_pars_t::fix(conf_t &conf)
{
  double diff=find_gaugefixing(conf);
  conf.gauge_transf(generate_sun(get_pars()));
  return diff;
}
