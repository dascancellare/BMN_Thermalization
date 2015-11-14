#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <cstdarg>

using namespace Eigen;
using namespace std;

//define the type
const int N=8,n=N-1;
const int nX=3;
int glb_N;
typedef Matrix<complex<double>,N,N> matr_t;

//length of trajectory
double T;
double t=0;
double dt;

//parameters of the initial conditions
double mass;
double hx,hy;
double v;

//initial format
int init_format;
const int init_static=0,init_angular=1;
double ang_J;
vector<double> ang_Z,ang_A;

//mass squared(depending on time)
inline double sqm(double t)
{return mass*mass;}

//random stuff
int seed;
mt19937_64 gen;

//degrees of freedom
vector<matr_t> X,P;

//imaginary unit
complex<double> I(0.0,1.0);

////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define crash(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

void internal_crash(int line,const char *file,const char *temp,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);
  
  cerr<<"ERROR at line "<<line<<" of file "<<file<<": "<<buffer<<endl;
  exit(1);
}

//read an element from input file
template <class T> void read(T &out,ifstream &in,string is)
{
  string s;
  if(!(in>>s)) crash("impossible to read expect string \"%s\"",is.c_str());
  if(s!=is) crash("obtained %s while reading %s",s.c_str(),is.c_str());
  if(!(in>>out)) crash("reading data");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//comutation between two matrices
template <typename Da,typename Db> auto comm(const MatrixBase<Da> &a,const MatrixBase<Db> &b) -> decltype(a*b-b*a)
{return a*b-b*a;}

//define the square
template <class T> T square(T a)
{return a*a;}

//return a complex gaussian with standard deviation sqrt(h/(2n))
complex<double> get_gauss(double h)
{
  normal_distribution<double> gauss(0,sqrt(h/n));
  return gauss(gen)+I*gauss(gen);
}

ofstream kinetic_energy_file("kinetic_energy");
ofstream common_potential_file("common_potential");
ofstream mass_potential_file("mass_potential");
ofstream energy_file("energy");
ofstream constraint_file("constraint");
ofstream trace_file("trace");
ofstream eigenvalues_x0_file("eigenvalues_x0");
ofstream eigenvalues_x1_file("eigenvalues_x1");
ofstream eigenvalues_y0_file("eigenvalues_y0");

//define the matrices of su(2) representation of dimension j
vector<Matrix<complex<double>,Dynamic,Dynamic> > generate_L(int dimrep)
{
  double j=(dimrep-1.0)/2;
  
  vector<Matrix<complex<double>,Dynamic,Dynamic> > J(3);
  for(auto &Ji : J) Ji.resize(dimrep,dimrep);
  
  Matrix<complex<double>,Dynamic,Dynamic> JPlus(dimrep,dimrep);
  JPlus.setZero();
  
  for(int indn=0;indn+1<dimrep;indn++) JPlus(indn,indn+1)=sqrt((indn+1)*(2*j-indn));
  cout<<JPlus<<endl;
  
  J[1]=  0.5*(JPlus+JPlus.transpose());
  J[2]=-I*0.5*(JPlus-JPlus.transpose());
  
  J[0].setZero();
  for(int indn=0;indn<dimrep;indn++) J[0](indn,indn)=j-indn;
  
  return J;
}

//put everything to random and then overwrite with L
void generate_static()
{
  //First is zero but for L
  X[0].setZero();
  
  //riempire con L
  auto J=generate_L(n);
  for(int i=0;i<nX;i++) X[i].block(0,0,n,n)=J[i];
  
  //original value for bottom-right corner of X[0]
  X[0](N-1,N-1)=X[0](N-2,N-2)-1.0;
  
  //metti a 0 le P tranne p[0]
  for(int i=0;i<glb_N;i++) P[i].setZero();
  P[0](N-1,N-1)=v;
}

//solution with angular momenta
void generate_angular()
{
  //put everything to zero
  for(int i=0;i<nX;i++) X[i].setZero();
  
  //fill X[0], X+ and P+
  matr_t XPlus,PPlus;
  for(int ic=0;ic<N;ic++)
    {
      X[0](ic,ic)=ang_Z[ic];
      XPlus(ic,(ic+1)%N)=ang_A[ic]; //there is a mismatch of a factor 2 because of the way X[1] and X[2] are defined in terms of XPlus
      PPlus(ic,(ic+1)%N)=ang_A[ic]*I*ang_J/(4*N*square(ang_A[ic]));
    }
  
  //define X[1] and X[2]
  X[1]=    XPlus+XPlus.transpose();
  X[2]=-I*(XPlus-XPlus.transpose());
  
  //define P
  P[0].setZero();
  P[1]=    PPlus+PPlus.adjoint();
  P[2]=-I*(PPlus-PPlus.adjoint());
}

//wrapper
void generate_matrices()
{
  for(auto &Xi : X) Xi.setZero();
  
  //fluttuazione
  for(int i=0;i<glb_N;i++) //first X is 0 apart from L
    {
      double h;
      if(i<nX) h=hx;
      else     h=hy;
      for(int ir=0;ir<N;ir++)
	for(int ic=ir+1;ic<N;ic++)
	  {
	    complex<double> delta_y=get_gauss(h);
	    X[i](ir,ic)=delta_y;
	    X[i](ic,ir)=conj(delta_y);
	  }
    }
  
  //discriminate
  switch(init_format)
    {
    case init_static:generate_static();break;
    case init_angular:generate_angular();break;
    default: crash("unknown init format %d!",init_format);break;
    }
}

//return the trace of the square
template <typename D> double trace_square(const MatrixBase<D> &M)
{return (M*M).trace().real();}

//compute the potential due to mass term
double mass_potential()
{
  double V=0;
  
  //potential X
  for(int i=0;i<nX;i++) V+=sqm(t)*(trace_square(X[i])+
				   2.0*(I*X[i]*comm(X[(i+1)%nX],X[(i+2)%nX])).trace().real());
  
  //potential Y
  for(int a=nX;a<glb_N;a++) V+=sqm(t)*trace_square(X[a])/4;
  
  return V/2;
}

//piece of the potential common to x and y
double common_potential()
{
  double V=0;
  
  //common part
  for(int a=0;a<glb_N;a++)
    for(int b=0;b<glb_N;b++)
      V+=-trace_square(comm(X[a],X[b]))/2;
  
  return V/2;
}

double potential()
{
  double V=0;
  
  //piece coming from mass
  V+=mass_potential();
  
  //add the common part (commutator)
  V+=common_potential();
  
  return V;
}

//compute the kinetic energy
double kinetic_energy()
{
  //momenta
  double K=0;
  for(int i=0;i<glb_N;i++) K+=trace_square(P[i]);
  K/=2;
  
  return K;
}

//compute hamiltonian
double hamiltonian()
{
  double K=kinetic_energy();
  double V=potential();
  
  return K+V;
}

//compute the kinetic energy of the trace of the momenta
double kinetic_energy_trace()
{
  //momenta
  double K=0;
  for(int i=0;i<glb_N;i++) K+=square(P[i].trace().real())/N;
  K/=2;
  
  return K;
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
  kinetic_energy_file<<t<<" "<<kinetic_energy()-kinetic_energy_trace()<<endl;
  common_potential_file<<t<<" "<<common_potential()<<endl;
  mass_potential_file<<t<<" "<<mass_potential()<<endl;
  energy_file<<t<<" "<<hamiltonian()<<endl;
  constraint_file<<t<<" "<<constraint()<<endl;
  trace_file<<t<<" "<<X[0].trace().real()<<endl;
  eigenvalues_x0_file<<t<<" "<<eigenvalues(X[0])<<endl;
  eigenvalues_x1_file<<t<<" "<<eigenvalues(X[1])<<endl;
  eigenvalues_y0_file<<t<<" "<<eigenvalues(X[3])<<endl;
}

//read the parameters for path
void read_angular_pars(string path)
{
  ifstream in(path);
  if(!in.good()) crash("unable to open path: %s",path.c_str());
  
  //read J
  if(!(in>>ang_J)) crash("unable to read J");
  
  //resize vectors
  ang_A.resize(N);
  ang_Z.resize(N);
  
  //read A
  for(int i=0;i<N;i++) if(!(in>>ang_A[i])) crash("unable to read A[%d]",i);
  
  //read Z
  ang_Z[N-1]=0;
  for(int i=0;i<N-1;i++)
    {
      if(!(in>>ang_Z[i])) crash("unable to read Z[%d]",i);
      ang_Z[N-1]-=ang_Z[i];
    }
}

int main()
{
  //open the input
  ifstream input("input");
  if(!input.good()) crash("unable to open \"input\"");
  
  read(glb_N,input,"glb_N");
  read(seed,input,"seed");
  read(T,input,"T");
  read(dt,input,"dt");
  read(mass,input,"mass");
  read(hx,input,"hx");
  read(hy,input,"hy");
  read(v,input,"v");
  read(init_format,input,"init_format");
  if(init_format==init_angular)
    {
      string path;
      read(path,input,"path");
      read_angular_pars(path);
    }
  
  //set things that depends on glb_N
  if(glb_N<=3) crash("glb_N must be >3");
  X.resize(glb_N);
  P.resize(glb_N);
  
  //init random generator device
  gen.seed(seed);
  
  kinetic_energy_file.precision(16);
  common_potential_file.precision(16);
  mass_potential_file.precision(16);
  energy_file.precision(16);
  constraint_file.precision(16);
  trace_file.precision(16);
  eigenvalues_x0_file.precision(16);
  eigenvalues_x1_file.precision(16);
  eigenvalues_y0_file.precision(16);
  
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
