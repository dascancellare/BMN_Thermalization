#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace Eigen;
using namespace std;

//define the type
const int N=8,n=N-1;
const int glb_N=9,nX=3;
typedef Matrix<complex<double>,N,N> matr_t;

//length of trajectory
double T=100;
double t=0;
double dt=0.001;

//parameters of the initial conditions
const double h=0.001;
const double v=20;
//const double h=0.00;
//const double v=0;

//mass squared(depending on time)
inline double sqm(double t)
{return 2;}

//random stuff
int seed=51768324;
mt19937_64 gen(seed);
normal_distribution<double> gauss(0,sqrt(h/n));

//degrees of freedom
vector<matr_t> X(glb_N),P(glb_N);
// need to define Jplus

//imaginary unit
complex<double> I(0.0,1.0);

//comutation between two matrices
template <typename Da,typename Db> auto comm(const MatrixBase<Da> &a,const MatrixBase<Db> &b) -> decltype(a*b-b*a)
{return a*b-b*a;}

//define the square
template <class T> T square(T a)
{return a*a;}

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss()
{return gauss(gen)+I*gauss(gen);}

ofstream energy_file("energy");
ofstream constraint_file("constraint");
ofstream trace_file("trace");
ofstream eigenvalues_x0_file("eigenvalues_x0");

//define the matrices of su(2) representation of dimension j
vector<Matrix<complex<double>,Dynamic,Dynamic> > generate_L(int dimrep)
{
  double j=(dimrep-1.0)/2;
  
  vector<Matrix<complex<double>,Dynamic,Dynamic> > J(3);
  for(auto &Ji : J) Ji.resize(dimrep,dimrep);
  
  Matrix<complex<double>,Dynamic,Dynamic> Jplus(dimrep,dimrep);
  Jplus.setZero();
  
  for(int indn=0;indn+1<dimrep;indn++) Jplus(indn,indn+1)=sqrt((indn+1)*(2*j-indn));
  
  J[1]=  0.5*(Jplus+Jplus.transpose());
  J[2]=-I*0.5*(Jplus-Jplus.transpose());
  
  J[0].setZero();
  for(int indn=0;indn<dimrep;indn++) J[0](indn,indn)=j-indn;
  
  return J;
}

//put everything to random and then overwrite with L
void generate_matrices()
{
  for(auto &Xi : X) Xi.setZero();
  
  //fluttuazione
  for(int i=1;i<glb_N;i++) //first X is 0 apart from L
    for(int ir=0;ir<N;ir++)
      for(int ic=ir+1;ic<N;ic++)
	{
	  complex<double> delta_y=get_gauss();
	  X[i](ir,ic)=delta_y;
	  X[i](ic,ir)=conj(delta_y);
	}
  
  //riempire con L
  auto J=generate_L(n);
  for(int i=0;i<nX;i++) X[i].block(0,0,n,n)=J[i];
  
  //metti a 0 le P tranne p[0]
  for(int i=1;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

//return the trace of the square
template <typename D> double trace_square(const MatrixBase<D> &M)
{return (M*M).trace().real();}

double potential()
{
  double V=0;
  
  //potential X
  for(int i=0;i<nX;i++) V+=sqm(t)*(trace_square(X[i])+
				   2.0*(I*X[i]*comm(X[(i+1)%nX],X[(i+2)%nX])).trace().real());
  
  //potential Y
  for(int a=nX;a<glb_N;a++) V+=sqm(t)*trace_square(X[a])/4;
  
  //common part
  for(int a=0;a<glb_N;a++)
    for(int b=0;b<glb_N;b++)
      V+=-trace_square(comm(X[a],X[b]))/2;
  
  //add final 1/2
  V/=2;
  
  return V;
}

//compute hamiltonian
double hamiltonian()
{
  //momenta
  double K=0;
  for(int i=0;i<glb_N;i++) K+=trace_square(P[i]);
  K/=2;
  
  double V=potential();
  
  return K+V;
}

//update the positions on the base of momenta
void update_positions(double dt)
{for(int i=0;i<glb_N;i++) X[i]+=P[i]*dt;}

void update_momenta(double dt)
{
  // calcola la forza
  vector<matr_t> F(glb_N);
  for(auto &Fi : F) Fi.setZero();
  
  //primo pezzo diverso per X ed Y
  for(int i=0;i<nX;i++) F[i]=sqm(t)*(-X[i]-3.0*I*comm(X[(i+1)%nX],X[(i+2)%nX]));
  for(int a=nX;a<glb_N;a++) F[a]=-0.25*sqm(t)*X[a];
  
  //secondo pezzo uguale per tutti
  for(int i=0;i<glb_N;i++) for(int a=0;a<glb_N;a++) F[i]+=comm(comm(X[a],X[i]),X[a]);
  
  //#define DEBUG
  
#ifdef DEBUG
  int icheck=1;
  double e_bef=hamiltonian();

  double eps=1e-6;
  X[icheck](1,0)+=eps/2;
  X[icheck](0,1)+=eps/2;
  
  double e_aft=hamiltonian();
  double f_num=-(e_aft-e_bef)/eps;
  double f_exa=(F[icheck](1,0)+F[icheck](0,1)).real()/2.0;
  cout<<"Exa: "<<f_exa<<", num: "<<f_num<<endl;
  
  X[icheck](1,0)-=eps;
  X[icheck](0,1)-=eps;
#endif
  
  ////
  
  //aggiorna i momenti
  for(int i=0;i<glb_N;i++) P[i]+=F[i]*dt;
}

//main leapfrog integration step
void integration_step()
{
  update_positions(dt);
  update_momenta(dt);
}

//compute the constraint (trace of the square)
double constraint()
{
  matr_t C;
  C.setZero();
  
  for(int i=0;i<glb_N;i++) C+=comm(X[i],P[i]);
  
  return trace_square(C);
}

//compute the eignevalues of X[0]
SelfAdjointEigenSolver<matr_t> es;
template <typename D> auto eigenvalues(const MatrixBase<D> &x) ->decltype(es.eigenvalues().transpose())
{
  es.compute(x);
  return es.eigenvalues().transpose();
}

//perform all measurement
void measure_observables()
{
  energy_file<<t<<" "<<hamiltonian()<<endl;
  constraint_file<<t<<" "<<constraint()<<endl;
  trace_file<<t<<" "<<X[0].trace().real()<<endl;
  eigenvalues_x0_file<<t<<" "<<eigenvalues(X[0])<<endl;
}

int main()
{
  generate_L(1);
  
  energy_file.precision(16);
  constraint_file.precision(16);
  trace_file.precision(16);
  eigenvalues_x0_file.precision(16);
  
  //generate random matrices+L
  generate_matrices();
  
  //compute number of steps needed to integrate
  int nt=T/dt;
  cout<<"Number of integration steps: "<<nt<<" to integrate "<<T<<" in steps of "<<dt<<endl;
  
  //first half step: p(dt/2)= p(0)+ F(0)*dt/2
  update_momenta(dt/2);
  
  //integrate equation of motion using leapfrog
  for(int it=0;it<nt;it++)
    {
      cout<<"Integration step "<<it+1<<"/"<<nt<<endl;
      
      update_positions(dt);
      if(it%10==0)
	{
	  update_momenta(dt/2);
	  measure_observables();
	  update_momenta(dt/2);
	}
      else update_momenta(dt);
      
      t+=dt;
    }
  
  return 0;
}
