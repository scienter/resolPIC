#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include <math.h>
#include <mpi.h>

void RungeKutta(double *W,double *prob,int iter,double dt,int start,int end,int flag)
{
  int n,Z;
  double ddt,N,w,oldN,oldW,oldWN,k1,k2,k3,k4;

  ddt=dt/((double)iter);
  //initialization of prob[Z]
  prob[start]=1.0; for(Z=start+1; Z<end; Z++) prob[Z]=0.0;

  n=0;
  while(n<iter)
  {
    //at start
    N=prob[start], w=W[start];
    k1=-w*N;
    k2=-w*(N+0.5*ddt*k1);
    k3=-w*(N+0.5*ddt*k2);
    k4=-w*(N+ddt*k3);   
    prob[start]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
    //rest of ..
    Z=start+1;
    while(Z<end) {
      oldWN=N*w; N=prob[Z], w=W[Z];
      k1=oldWN-w*N;
      k2=oldWN-w*(N+0.5*ddt*k1);
      k3=oldWN-w*(N+0.5*ddt*k2);
      k4=oldWN-w*(N+ddt*k3);
      prob[Z]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
      Z++;
    }
//    //last
//    k1=w*N;
//    k2=w*(N+0.5*ddt*k1);
//    k3=w*(N+0.5*ddt*k2);
//    k4=w*(N+ddt*k3);   
//    prob[end]=N+ddt/6.0*(k1+2.0*k2+2.0*k3+k4); 
    n++;
  }

  Z=start+1;
  while(Z<end) {
   prob[Z]+=prob[Z-1];
   Z++;
  }
//  if(flag==1) {
//    for(Z=start; Z<=end; Z++) 
//    printf("prob[%d]=%g, ",Z,prob[s][Z]);
//    printf("\n");
//  }
}

double randomValue(double beta) 
{ 
   double r; 
   int intRand, randRange=10000000, rangeDev; 
 
   rangeDev=(int)(randRange*(1.0-beta)); 
   intRand = rand() % (randRange-rangeDev); 
   r = ((double)intRand)/randRange+(1.0-beta); 
 
   return r; 
} 

