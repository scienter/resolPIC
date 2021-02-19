#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


void loadPlasmaLens2D(Domain *D,PlasmaLens *PL,int iteration);
void loadPlasmaLens3D(Domain *D,PlasmaLens *PL,int iteration);

void plasmaLens(Domain *D,PlasmaLens *PL,int iteration)
{
  switch(D->dimension)  {
  //1D
  case 1:
    ;
    break;
  //2D
  case 2:
    loadPlasmaLens2D(D,PL,iteration);
    break;
  //3D
  case 3:
    loadPlasmaLens3D(D,PL,iteration);
    break;
  default:
    ;
  }
}

void loadPlasmaLens2D(Domain *D,PlasmaLens *PL,int iteration)
{
   int s,n,i,j,k,istart,iend,jstart,jend,kstart,kend,cnt,l,ii,index;
   int minYSub;
   double ne,y,I0,xi,unitXi,Bphi,normalB,m_I,coef,posX,sign;
   double c[6],x[6];
   int ss[D->nPlasmaLens];
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;
   LoadList *LL;

   ptclList *p;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   k=0;
   minYSub=D->minYSub;

   c[0]=0.3155/2.0;
   c[1]=-0.0233967/3.0;
   c[2]=0.0786747/4.0;
   c[3]=-0.878334/5.0;
   c[4]=1.32213/6.0;
   c[5]=-0.751897/7.0;
   unitXi=D->dy*D->lambda/PL->radius;
   normalB=1.0/(eMass*D->omega/(-eCharge));
   m_I=0.106493;
   coef=mu0*PL->current*0.5/M_PI/PL->radius/m_I;

   LL=D->loadList;
   index=0;
   s=0;
   while(LL->next)  {
     if(LL->type==Beam) {
	    ss[index]=s;
	    index++;
     } else ;
     s++;
     LL=LL->next;
   }

   //position define      
   for(i=istart; i<iend; i++) 
   {
     posX=(double)(i+D->minXSub-istart);
     for(l=0; l<PL->xnodes-1; l++) {
       if(posX>=PL->xpoint[l] && posX<PL->xpoint[l+1]) 
         ne=((PL->xn[l+1]-PL->xn[l])/(PL->xpoint[l+1]-PL->xpoint[l])*(posX-PL->xpoint[l])+PL->xn[l]);
       else
	      ne=1.0;
     }
     for(j=jstart; j<jend; j++)
     {
       for(ii=0; ii<D->nPlasmaLens; ii++)  {
//       for(s=0; s<D->nPlasmaLens; s++)  {
//         p=particle[i][j][k].head[s]->pt;
         p=particle[i][j][k].head[ss[ii]]->pt;
         while(p)  {
           y=(p->y+j-jstart+minYSub); // z=p->z;
	        xi=fabs(y*unitXi);
	        x[0]=xi;
           x[1]=x[0]*xi;
           x[2]=x[1]*xi;
           x[3]=x[2]*xi;
           x[4]=x[3]*xi;
           x[5]=x[4]*xi;

	        Bphi=0.0;
	        for(n=0; n<6; n++) Bphi+=c[n]*x[n];
	        Bphi*=normalB*coef;
//	   printf("xi=%g, Bphi=%g\n",xi,Bphi);
//	        if(y>=0.0) sign=1.0;
//           else       sign=-1.0;	   
           sign=y/fabs(y);
           p->B3+=sign*Bphi*ne;

           p=p->next;
	      }
       }
     }
   }
}

void loadPlasmaLens3D(Domain *D,PlasmaLens *PL,int iteration)
{
   int s,n,i,j,k,istart,iend,jstart,jend,kstart,kend,cnt,l,ii,index;
   int minYSub,minZSub;
   double ne,y,z,r,I0,xi,unitXi,Bphi,normalB,m_I,coef,posX;
   double c[6],x[6],cosTh,sinTh;
   int ss[D->nPlasmaLens];
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;
   LoadList *LL;

   ptclList *p;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minYSub=D->minYSub;
   minZSub=D->minZSub;

   c[0]=0.3155/2.0;
   c[1]=-0.0233967/3.0;
   c[2]=0.0786747/4.0;
   c[3]=-0.878334/5.0;
   c[4]=1.32213/6.0;
   c[5]=-0.751897/7.0;
   unitXi=D->dy*D->lambda/PL->radius;
   normalB=1.0/(eMass*D->omega/(-eCharge));
   m_I=0.106493;
   coef=mu0*PL->current*0.5/M_PI/PL->radius/m_I;

   LL=D->loadList;
   index=0;
   s=0;
   while(LL->next)  {
     if(LL->type==Beam) {
	    ss[index]=s;
	    index++;
     } else ;
     s++;
     LL=LL->next;
   }

   //position define      
   for(i=istart; i<iend; i++) 
   {
     posX=(double)(i+D->minXSub-istart);
     for(l=0; l<PL->xnodes-1; l++) {
       if(posX>=PL->xpoint[l] && posX<PL->xpoint[l+1]) 
         ne=((PL->xn[l+1]-PL->xn[l])/(PL->xpoint[l+1]-PL->xpoint[l])*(posX-PL->xpoint[l])+PL->xn[l]);
       else
	      ne=1.0;
     }
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++) 
         for(ii=0; ii<D->nPlasmaLens; ii++)  {
//       for(s=0; s<D->nPlasmaLens; s++)  {
//         p=particle[i][j][k].head[s]->pt;
           p=particle[i][j][k].head[ss[ii]]->pt;
           while(p)  {
             y=(p->y+j-jstart+minYSub); 
             z=(p->z+k-kstart+minZSub); 
//				 absY=fabs(y);
//				 absZ=fabs(z);
             r=sqrt(y*y+z*z);
	          xi=r*unitXi;
	          x[0]=xi;
             x[1]=x[0]*xi;
             x[2]=x[1]*xi;
             x[3]=x[2]*xi;
             x[4]=x[3]*xi;
             x[5]=x[4]*xi;

             Bphi=0.0;
	          for(n=0; n<6; n++) Bphi+=c[n]*x[n];
	          Bphi*=normalB*coef;
				 cosTh=y/r;
				 sinTh=z/r;
             p->B2+=cosTh*Bphi*ne;
             p->B3+=sinTh*Bphi*ne;

             p=p->next;
	        }
         }
   }
}

