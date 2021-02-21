#include<math.h>
#include<stdio.h>
#include"funcs_comov.h"


double omega_m,omega_l,omega_k,vhh,cbyho,LL;
int curv;

// Omega_m0, Omega_k0, H0/100 are needed as input
void initialize(double omega_m1,double omega_k1,double hh)
{
  vhh=hh;
  omega_m=omega_m1;
  omega_k=omega_k1;
  omega_l=1.-omega_m-omega_k;
  cbyho=3.e3/hh;    /*coH0 =c/H_0 in terms of Mpc*/

  curv=0; // - 1 if open 
  if(fabs(omega_k)>1.e-8)
    {
      LL=cbyho/sqrt(fabs(omega_k)); // radius of curvature 

      if(omega_k>0.) curv=-1; // open
      else curv=1; // closed
    }
}

/*----------------Returns E= H/H_0------------*/ 

double func_E(double x)
{
  double tt;
  tt=omega_m*powf(x,-3)+ omega_k*powf(x,-2.) + omega_l;
  return(sqrt(tt));
}

// used for ceta

double func( double x)
{
  return(powf(func_E(x)*x*x,-1.));
}

//calculate ceta for given scale factor

double ceta(double x)
{
  double func(double);
  double y;
  y=simp(func,10000,x,1.);
  y=y*cbyho;       
  return(y);
}

// which is conformal distance r for given scale factor  
// r=d_A *(1+z) where d_A is angular diamter distance

double rnu(double x)
{
  double ceta(double x);
    
    if(curv==0)
      {
	return(ceta(x));
      }
    else
      {
	if(curv==-1)	return(LL*sinh(ceta(x)/LL));
	else	return(LL*sin(ceta(x)/LL));

      }
}


// nu derivative of r at given scale factor

double rnup(double x)
{
  double term;

  term=cbyho/(x*x*func_E(x)*1420.);

  if(curv==-1)
    {
      term=term*cosh(ceta(x)/LL);
    }

  if(curv==1)
    {
      term=term*cos(ceta(x)/LL);
    }
  return(term);
}



