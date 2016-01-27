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
    
    //cerr<<"==============================================================="<<endl<<"  Tranforming matrix: "<<endl<<generate_sun(w)<<endl<<"  orig matrix: "<<endl<<(*transforming).X[0]<<endl<<"  transformed: "<<endl<<transforming->get_gauge_transformed(generate_sun(w)).X[0]<<endl<<"  functional: "<<out<<endl;
  }
}

gauge_fix_pars_t::gauge_fix_pars_t(conf_t ref_conf) : ref_conf(ref_conf)
{
  //initialize vector of pars
  check_gln_N_set();
  this->resize(generators.size());
  for(auto &it : *this) it=0;
  
  //init minuit
  minu=new TMinuit(generators.size());
  minu->SetFCN(dist);
  minu->SetPrintLevel(-1);
  for(size_t i=0;i<generators.size();i++) minu->DefineParameter(i,"par",(*this)[i],0.01,0,0);
  
  //minize after setting tolerance
  double tol=1e-16;
  int iflag;
  minu->SetMaxIterations(10000000);
  minu->mnexcm("SET ERR",&tol,1,iflag);
}

void gauge_fix_pars_t::find_gaugefixing(conf_t &conf)
{
  //set the reference
  ref_fix=&ref_conf;
  transforming=&conf;
  
  //cerr<<"Fixing: "<<endl;
  //cerr<<transforming->X[0]<<endl;
  //cerr<<"Against: "<<endl;
  //cerr<<ref_fix->X[0]<<endl;
  
  //evaluate ch2
  //double pars[generators.size()];
  //for(size_t i=0;i<generators.size();i++) pars[i]=0;
  //double ch2_bef;
  //minu->Eval(generators.size(),NULL,ch2_bef,pars,0);
  
  //minimize
  minu->Migrad();
  
  //get pars
  double dum;
  for(size_t i=0;i<generators.size();i++)
    //{
      minu->GetParameter(i,(*this)[i],dum);
      //pars[i]=(*this)[i];
  //}
  
  //double ch2_aft;
  //minu->Eval(generators.size(),NULL,ch2_aft,pars,0);
  //cerr<<"Difference in ch2: "<<ch2_aft-ch2_bef<<endl;
  //cerr<<" "<<(*this)[0]<<endl;
  //cerr<<generate_sun((*this))<<endl;
}

void gauge_fix_pars_t::fix(conf_t &conf)
{
  find_gaugefixing(conf);
  //conf.gauge_transf(generate_sun((*this)));
}
