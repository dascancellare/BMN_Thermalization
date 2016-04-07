#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>

#define EXTERN_RANDOM
#include "random.hpp"

#include "matr.hpp"

double get_real_rand_gauss(double h)
{
  normal_distribution<double> gauss(0,sqrt(h));
  return gauss(gen);
}

complex<double> get_comp_rand_gauss(double h)
{return get_real_rand_gauss(h)+I*get_real_rand_gauss(h);}

double get_rand_double(double min,double max)
{
  uniform_real_distribution<double> unif(min,max);
  return unif(gen);
}
