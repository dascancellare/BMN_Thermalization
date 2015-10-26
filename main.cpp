#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

double T=10;
double dt=0.1;

//define the type
const int N=4;
const int glb_N=9;
typedef Matrix<complex<double>,N,N> matr_t;

vector<matr_t> X(glb_N),P(glb_N);

//
void generate_random_matrices()
{
}

void integration_step()
{
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
  
  //first half step
  
  //integrate equation of motion
  for(int it=0;it<nt;it++)
    {
      integration_step();
      
      measure_observables();
    }
  
  return 0;
}
