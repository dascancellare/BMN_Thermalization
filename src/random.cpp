#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>

#define EXTERN_RANDOM
#include "random.hpp"

#include "matr.hpp"

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss(double h)
{
  normal_distribution<double> gauss(0,sqrt(h/n));
  return gauss(gen)+I*gauss(gen);
}

