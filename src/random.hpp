#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
#endif

#include <complex>
#include <random>

using namespace std;

EXTERN_RANDOM int seed;
EXTERN_RANDOM mt19937_64 gen;

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss(double h);

#undef EXTERN_RANDOM

#endif
