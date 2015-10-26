#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

//length of trajectory
double T=10;
double dt=0.1;
//parameters of the initial conditions
double h;
double v;
int seed=57683245;

//define the type
const int N=4;
const int glb_N=9;
typedef Matrix<complex<double>,N,N> matr_t;

vector<matr_t> X(glb_N),P(glb_N);

//
void generate_random_matrices()
{
}



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
  cout<<X[0]<<endl;
  
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
