#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>

#define EXTERN_RANDOM
#include "random.hpp"

#include "matr.hpp"

complex<double> get_rand_gauss(double h)
{
  normal_distribution<double> gauss(0,sqrt(h/n));
  return gauss(gen)+I*gauss(gen);
}

double get_rand_double(double min,double max)
{
  uniform_real_distribution<double> unif(min,max);
  return unif(gen);
}
