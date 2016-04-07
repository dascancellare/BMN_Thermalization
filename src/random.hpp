#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
#endif

#include <complex>
#include <random>

using namespace std;

EXTERN_RANDOM mt19937_64 gen;

//! return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_comp_rand_gauss(double h);
double get_real_rand_gauss(double h);

//! return a double in the range [min,max)
double get_rand_double(double min,double max);

#undef EXTERN_RANDOM

#endif
