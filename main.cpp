#include <Eigen/Dense>

double T=10;
double dt=0.1;

int N=4;

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
  //generate random variables
  generate_random_matrices();
  
  //compute number of steps needed to integrate
  int nt=T/dt;
  
  //integrate equation of motion
  for(int it=0;it<nt;it++)
    {
      integration_step();
      
      measure_observables();
    }
  
  return 0;
}
