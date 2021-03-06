#ifndef _MATR_HPP
#define _MATR_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_MATR
 #define EXTERN_MATR extern
 #define ONLY_INSTANTIATION
#endif

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#include <vector>
#include "tools.hpp"

using namespace Eigen;
using namespace std;

/////////////////////////////////////// types /////////////////////////////////

//! type of matrix
const int nX=3; //do not touch
const int glb_N=9;

typedef Matrix<complex<double>,NCOL,NCOL> matr_t;

/////////////////////////////////////// globals ///////////////////////////////

//! imaginary unit
#define I complex<double>(0.0,1.0)

//! generators
EXTERN_MATR vector<matr_t> generators;

///////////////////////////////////// prototypes //////////////////////////////

//! fill the generators of SU(n)
template <int n> vector<Matrix<complex<double>,n,n>> get_generators()
{
  vector<Matrix<complex<double>,n,n>> out;
  
  //real generators
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      {
	Matrix<complex<double>,n,n> gen;
	gen.Zero();
	gen(i,j)=gen(j,i)=1;
	out.push_back(gen);
      }
  
  //diagonal generators
  for(int i=1;i<n;i++)
    {
      Matrix<complex<double>,n,n> gen;
      gen.Zero();
      double norm=1/sqrt(i*(i+1)/2);
      for(int j=0;j<i;j++) gen(j,j)=norm;
      gen(i,i)=-i*norm;
      out.push_back(gen);
    }
  
  //imag generators
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      {
	Matrix<complex<double>,n,n> gen;
	gen.Zero();
	gen(i,j)=+I;
	gen(j,i)=-I;
	out.push_back(gen);
      }
  
  return out;
}

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

//! anticommutation between two matrices
template <typename Da,typename Db> auto anti_comm(const MatrixBase<Da> &a,const MatrixBase<Db> &b) -> decltype(a*b+b*a)
{return a*b+b*a;}

//! define the square
template <class T> T square(T a)
{return a*a;}

//! return the trace of the square
template <typename D> double trace_square(const MatrixBase<D> &M)
{return (M*M).trace().real();}

//! return the trace of the fourth power
template <typename D> double trace_fourth(const MatrixBase<D> &M)
{return (M*M*M*M).trace().real();}

//! return the norm2 of the matrix
template <typename D> double trace_norm2(const MatrixBase<D> &M)
{return (M*M.adjoint()).trace().real();}

//! check that glb_N has been fixed
inline void check_gln_N_set()
{if(glb_N<0) CRASH("please init glb_N before");}

//! generate a matrix on the base of the arguments
inline auto generate_sun(vector<double> w) -> decltype(generators[0].exp())
{
  matr_t res;
  res.setZero();
  
  if(generators.size()!=w.size()) CRASH("generators and coefficients have different size, %u vs %u",generators.size(),w.size());
  
  for(size_t i=0;i<generators.size();i++) res+=I*generators[i]*w[i];
  
  return res.exp();
}

#undef EXTERN_MATR

#endif
