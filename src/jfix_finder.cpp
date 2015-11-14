#include <TMinuit.h>
#include <fstream>
#include <iostream>

using namespace std;

double pref;
int N;

template <class T> T sqr(T in)
{return in*in;}

void energy(int &npar,double *fuf,double &out,double *p,int flag)
{
  double *a2=p;
  double z[N];
  z[N-1]=0;
  for(int i=0;i<N-1;i++)
    {
      z[i]=p[N+i];
      z[N-1]-=z[i];
    }
  
  //kinetic part: J^2/8N^2*sum(1/ai^2)
  double kin=0;
  
  for(int i=0;i<N;i++) kin+=1/a2[i];
  kin*=pref;
  
  //potential
  double pot=0;
  for(int i=0;i<N;i++)
    {
      int ip=(i+1)%N;
      int im=(i+N-1)%N;
      
      pot+=sqr(z[i]+2*a2[im]-2*a2[i])/2+
	2*sqr(1+z[ip]-z[i])*a2[i];
    }
  
  out=kin+pot;
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Use: "<<arg[0]<<" N"<<endl;
      exit(0);
    }
  N=atoi(arg[1]);
  
  TMinuit minu(2*N+1);
  minu.SetFCN(energy);
  minu.SetPrintLevel(-1);
  
  //define A
  double pars[2*N-1];
  for(int i=0;i<N;i++)
    {
      char a2tag[10];
      sprintf(a2tag,"a2%d",i);
      
      //fixed to jplus
      double j=(N-1.0)/2;
      double val=(i+1)*(2*j-i)/4;
      minu.DefineParameter(i,a2tag,val+1e-3,0.0005,0,1e3);
    }
  
  //define Z
  for(int i=0;i<N-1;i++)
    {
      char ztag[10];
      sprintf(ztag,"z%d",i);
      double val=(N-1.0)/2-i;
      minu.DefineParameter(N+i,ztag,val,0.0005,-1e6,1e6);
    }
  
  //storage
  vector<double> min_val,J;
  vector<vector<double> > A(N);
  vector<vector<double> > Z(N-1);
  
  //set tolerance
  double tol=1e-16;
  int iflag;
  minu.SetMaxIterations(10000000);
  minu.mnexcm("SET ERR",&tol,1,iflag);
  
  //increas J
  double step_J=0.01/2/sqr(N);
  double max_J=10;
  for(double glb_J=0;glb_J<max_J;glb_J+=step_J)
    {
      //store J and set prefactor
      J.push_back(glb_J);
      pref=sqr(glb_J)/(8*N*N);
      
      //minimize and store minimum
      if(glb_J>1e-6) minu.Migrad();
      
      //get back A and Z
      for(int i=0;i<N;i++)
	{
	  double temp,dum;
	  //A
	  minu.GetParameter(i,temp,dum);
	  pars[i]=temp;
	  A[i].push_back(sqrt(temp));
	  //Z
	  if(i!=N-1)
	    {
	      minu.GetParameter(i+N,temp,dum);
	      pars[i+N]=temp;
	      Z[i].push_back(temp);
	    }
	}
      
      //store minimum
      double val;
      minu.Eval(2*N-1,NULL,val,pars,0);
      min_val.push_back(val);
    }
  
  ofstream fout("sol");
  fout.precision(16);
  for(size_t ij=0;ij<J.size();ij++)
    {
      fout<<J[ij]<<"\t"<<min_val[ij]<<"\t";
      for(int k=0;k<N;k++) fout<<A[k][ij]<<"\t";
      fout<<"\t";
      for(int k=0;k<N-1;k++) fout<<Z[k][ij]<<"\t";
      fout<<endl;
    }
  
  return 0;
}
