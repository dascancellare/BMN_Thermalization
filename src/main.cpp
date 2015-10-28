#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

using namespace Eigen;
using namespace std;

//define the type
const int N=4,n=N-1;
const int glb_N=9,nX=3;
typedef Matrix<complex<double>,N,N> matr_t;

//length of trajectory
double T=10;
double t=0;
double dt=0.01;

//parameters of the initial conditions
const double h=0.001;
const double v=10;

//random stuff
int seed=57683245;
mt19937_64 gen(seed);
normal_distribution<double> gauss(0,sqrt(h/(2*n)));

//degrees of freedom
vector<matr_t> X(glb_N),P(glb_N);
// need to define Jplus

//imaginary uniti
complex<double> I(0.0,1.0);

//comutation between two matrices
matr_t comm(matr_t a,matr_t b)
{return a*b-b*a;}

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss()
{return gauss(gen)+I*gauss(gen);}

ofstream energy_file("energy");
//define the matrices of su(2) representation of dimension j
/*
void generate_L (int j)
{
  int dimrep=2j+1;
  vector<matr_t> Jplus(dimrep),J1(dimrep),J2(dimrep),J3(dimrep);
  Jplus.setZero();
  J1.setZero();
  J2.setZero();
  J3.setZero();
  
  for(int indn=0,indn+1<dimrep,indn++)
    Jplus(indn,indn+1)=sqrt((indn+1)(2j+indn));
  
  J1=0.5*( Jplus + transpose(Jplus));
  J2=-I*0.5*( Jplus - transpose(Jplus));

  for(int indn=0,ind<dimrep,ind++)
    J3(ind,ind)=j-ind;
}
*/
//put everything to random and then overwrite with L
void generate_matrices()
{
  //fluttuazione
  for(int i=1;i<glb_N;i++) //first X is 0 apart from L
    for(int ir=0;ir<N;ir++)
      {
	X[i](ir,ir)=0;
	for(int ic=ir+1;ic<N;ic++)
	  {
	    complex<double> delta_y=get_gauss();
	    X[i](ir,ic)=delta_y;
	    X[i](ic,ir)=conj(delta_y);
	  }
      }
  
  //riempire con L
  
  
  //metti a 0 le P tranne p[0]
  for(int i=1;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

double trace_square(matr_t M)
{return (M*M).trace().real();}

//compute hamiltonian
double hamiltonian()
{
  //momenta
  double K=0;
  for(int i=0;i<glb_N;i++) K+=trace_square(P[i]);
  K/=2;
  
  //potential X
  double V=0;
  for(int i=0;i<nX;i++) V+=trace_square(X[i]+I*comm(X[(i+1)%nX],X[(i+2)%nX]));
  
  //potential Y
  for(int a=nX;a<glb_N;a++) V+=trace_square(X[a])/4;
  
  //potential X-Y
  for(int i=0;i<nX;i++)
    for(int a=nX;a<glb_N;a++)
      V+=-trace_square(comm(X[i],X[a]));
  
  //potential Y-Y
  for(int a=nX;a<glb_N;a++)
    for(int b=nX;b<glb_N;b++)
      V+=-trace_square(comm(X[a],X[b]))/2;
  
  //add final 1/2
  V/=2;
  
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
  for(int i=0;i<nX;i++) F[i]=-X[i]-3.0*I*comm(X[(i+1)%nX],X[(i+2)%nX]);
  for(int a=nX;a<glb_N;a++) F[a]=-0.25*X[a];
  
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

matr_t constraint()
{
  matr_t C;
  C.setZero();
  
  for(int i=0;i<glb_N;i++) C+=comm(X[i],P[i]);
  
  return C;
}

//
void measure_observables()
{
  energy_file<<t<<" "<<hamiltonian()<<endl;
  //cout<<"Constraint: "<<endl<<constraint()<<endl<<endl;
}

int main()
{
  energy_file.precision(16);
  
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
      update_momenta(dt/2);
      
      measure_observables();
      
      update_momenta(dt/2);
      
      t+=dt;
    }
  
  return 0;
}
