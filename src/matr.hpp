#ifndef _MATR_HPP
#define _MATR_HPP

#ifndef EXTERN_MATR
 #define EXTERN_MATR extern
 #define ONLY_INSTANTIATION
#endif

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include "tools.hpp"

using namespace Eigen;
using namespace std;

/////////////////////////////////////// types /////////////////////////////////

//! type of matrix
const int N=6,n=N-1;
const int nX=3; //do not touch
EXTERN_MATR int glb_N
#ifndef ONLY_INSTANTIATION
=-1
#endif
;
typedef Matrix<complex<double>,N,N> matr_t;

/////////////////////////////////////// globals ///////////////////////////////

//! imaginary unit
extern complex<double> I;

//! generators
EXTERN_MATR vector<matr_t> generators;

///////////////////////////////////// prototypes //////////////////////////////

//! fill the generators of SU(N)
void fill_generators();

//! define the matrices of su(2) representation of dimension j
vector<Matrix<complex<double>,Dynamic,Dynamic> > generate_L(int dimrep);

///////////////////////////////////// templates ////////////////////////////////

//! compute the eignevalues of X[0]
EXTERN_MATR SelfAdjointEigenSolver<matr_t> es;
template <typename D> auto eigenvalues(const MatrixBase<D> &x) ->decltype(es.eigenvalues().transpose())
{
  es.compute(x);
  return es.eigenvalues().transpose();
}

//! commutation between two matrices
template <typename Da,typename Db> auto comm(const MatrixBase<Da> &a,const MatrixBase<Db> &b) -> decltype(a*b-b*a)
{return a*b-b*a;}

//! define the square
template <class T> T square(T a)
{return a*a;}

//! return the trace of the square
template <typename D> double trace_square(const MatrixBase<D> &M)
{return (M*M).trace().real();}

//! check that glb_N has been fixed
inline void check_gln_N_set()
{if(glb_N<0) CRASH("please init glb_N before");}

//! generate a matrix on the base of the arguments
inline auto generate_sun(vector<double> &w) -> decltype(generators[0].exp())
{
  matr_t res;
  res.setZero();
  
  if(generators.size()!=w.size()) CRASH("generators and coefficients have different size, %u vs %u",generators.size(),w.size());
  
  for(size_t i=0;i<generators.size();i++) res+=I*generators[i]*w[i];
  
  return res.exp();
}

#undef EXTERN_MATR

#endif
