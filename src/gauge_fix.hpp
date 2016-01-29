#ifndef _GAUGE_FIX_HPP
#define _GAUGE_FIX_HPP

#include <vector>

#include <TMinuit.h>

#include "conf.hpp"

using namespace std;

//! hold pars to fix a gauge
struct gauge_fix_pars_t
{
  conf_t ref_conf; //!< reference conf
  TMinuit *minu; //<! minimizer
  
  //! init fixing a reference
  gauge_fix_pars_t(conf_t ref_conf);
  
  //! get pars
  vector<double> get_pars();
  
  //! find the transformed configuration that maximize the trace w.r.t ref
  double find_gaugefixing(conf_t &ref);
  
  //! perform gauge fixing w.r.t store conf
  double fix(conf_t &conf);
  
  ~gauge_fix_pars_t(){delete minu;}
private:
  gauge_fix_pars_t(){}
};

#endif