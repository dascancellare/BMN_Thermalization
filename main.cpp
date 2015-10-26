#include <Eigen/Dense>
#include <complex>
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
double dt=0.1;

//parameters of the initial conditions
const double h=0.001;
const double v=10;

//random stuff
int seed=57683245;
mt19937_64 gen(seed);
normal_distribution<double> gauss(0,sqrt(h/(2*n)));

//degrees of freedom
vector<matr_t> X(glb_N),P(glb_N);

//imaginary uniti
complex<double> I(0.0,1.0);

//comutation between two matrices
matr_t comm(matr_t a,matr_t b)
{return a*b-b*a;}

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss()
{return gauss(gen)+I*gauss(gen);}

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

//update the positions on the base of momenta
void update_positions(double dt)
{for(int i=0;i<glb_N;i++) X[i]+=P[i]*dt;}

void update_momenta(double dt)
{
  // calcola la forza
  vector<matr_t> F(glb_N);
  
  //primo pezzo diverso per X ed Y
  for(int i=0;i<nX;i++) F[i]=-X[i]-6.0*I*comm(X[(i+1)%nX],X[(i+2)%nX]);
  for(int i=nX;i<glb_N;i++) F[i]=-0.25*X[i];
  
  //secondo pezzo uguale per tutti
  for(int i=0;i<glb_N;i++) for(int a=0;a<glb_N;a++) F[i]+=comm(comm(X[a],X[i]),X[i]);
  
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

double trace_square(matr_t M)
{return (M*M).trace().real();}

//compute hamiltonian
double hamiltonian()
{
  double H=0;
  
  //momenta
  for(int i=0;i<glb_N;i++) H+=(P[0]*P[0]).trace().real();
  
  //potential X
  for(int i=0;i<nX;i++) H+=trace_square(X[i]+I*comm(X[(i+1)%nX],X[(i+2)%nX]));
  
  //potential Y
  for(int a=nX;a<glb_N;a++) H+=trace_square(X[a]*X[a])/4;
  
  //potential Y-Y
  for(int i=0;i<nX;i++)
    for(int a=nX;a<glb_N;a++)
      H+=-trace_square(comm(X[i],X[a]));
  
  //potential Y-Y
  for(int a=nX;a<glb_N;a++)
    for(int b=nX;b<glb_N;b++)
      H+=-trace_square(comm(X[a],X[b]))/2;
  
  //add final 1/2
  H/=2;
  
  return H;
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
  cout<<X[0]<<endl;
  cout<<"Hamiltonian: "<<endl<<hamiltonian()<<endl<<endl;
  cout<<"Constraint: "<<endl<<constraint()<<endl<<endl;
}

int main()
{
  //generate random matrices+L
  generate_matrices();
  
  //compute number of steps needed to integrate
  int nt=T/dt;
  cout<<"Number of integration steps: "<<nt<<" to integrate "<<T<<" in steps of "<<dt<<endl;
  
  //first half step: p(dt/2)= p(0)+ F(0)*dt/2
  //update_momenta(dt/2);
  
  //integrate equation of motion using leapfrog
  for(int it=0;it<nt;it++)
    {
      cout<<"Integration step "<<it+1<<"/"<<nt<<endl;
      
      measure_observables();
      
      integration_step();
      
    }
  
  return 0;
}
