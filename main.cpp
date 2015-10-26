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

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss()
{return gauss(gen)+I*gauss(gen);}

//
void generate_matrices()
{
  //riempire con L
  
  //aggiungere fluttuazione sulle x
  for(int i=1;i<nX;i++)
    {
      complex<double> delta_x=get_gauss();
      X[i](0,N-1)=delta_x;
      X[i](N-1,0)=conj(delta_x);
    }
  
  //aggiungere fluttuazione sulle x
  for(int i=nX;i<glb_N;i++)
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
  
  //metti a 0 le P tranne p[0]
  for(int i=1;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

//update the positions on the base of momenta
void update_positions(double dt)
{
}

void update_momenta(double dt)
{
// calcola la forza
  vector<matr_t> F(glb_N);
  
//aggiorna i momenti
  
}

void integration_step()
{
  update_positions(dt);
  update_momenta(dt);
  
}

void measure_observables()
{
}

int main()
{
  //generate random variables
  generate_random_matrices();
  
  //compute number of steps needed to integrate
  int nt=T/dt;
  
  //first half step: p(dt/2)= p(0)+ F(0)*dt/2
  update_momenta(dt/2);
  
  
  //integrate equation of motion
  for(int it=0;it<nt;it++)
    {
      integration_step();
      
      measure_observables();
    }
  
  return 0;
}
