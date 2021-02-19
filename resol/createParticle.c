#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void storeParticle(Domain *D,int i,int j,int k,int s,double x,double y,double z,double weight);

void createParticles1D(Domain *D,int s)
{
   int i,j,k,l,n,gridNum=2;
   int istart,iend,jstart,jend,kstart,kend;
   double x,gapX,px,py,pz,maxW;
   double x0,px0,py0,pz0,weight,W;
   double Q[2],Px[2],Py[2],Pz[2];
   double w2[2],px2[2],py2[2],pz2[2],gapX2[2],x2[2];

   char fileName[100],dataName[100];
   Particle ***particle;
   particle=D->particle;
   queList *qHead,*q,*prev,*new,*tmp;
   ptclList *p,*head;

   qHead=(queList *)malloc(sizeof(queList ));
   qHead->next=NULL;

   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   kstart=D->kstart; kend=D->kend;
  
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++) 
       {
         Q[0]=D->q1[i][j][k];    Q[1]=D->q2[i][j][k];

         //first step
         W=0.0;
         for(n=0; n<gridNum; n++)  W+=Q[n];
    
         if(W>0)
         {
           x0=Q[1]/W;
           //second step
           w2[0]=Q[0];      w2[1]=Q[1];
           x2[0]=x0*0.5;    x2[1]=(1.0+x0)*0.5;
           gapX2[0]=x0*0.5;  gapX2[1]=(1.0-x0)*0.5;
           for(n=0; n<2; n++)  {
             new = (queList *)malloc(sizeof(queList ));
             new->next = qHead->next;
             qHead->next = new;
             new->x=x2[n];   new->weight=w2[n]; new->gapX=gapX2[n];
           }

           maxW=0.0;
           for(n=0; n<2; n++) { if(maxW<w2[n])  maxW=w2[n]; else; }
           n=(int)(maxW/D->targetW);
  
           for(l=0; l<n; l++)    {
             q=qHead->next;
             while(q)  {
               x=q->x;   weight=q->weight; gapX=q->gapX;
               if(weight>D->targetW) {
                 q->gapX=gapX*0.5; q->x=x-gapX*0.5; q->weight=weight*0.5;
                 tmp=q->next;
                 new = (queList *)malloc(sizeof(queList ));
                 new->x=x+gapX*0.5; new->weight=weight*0.5; new->gapX=gapX*0.5;
                 new->next=q->next;
                 q->next=new;
                 q=tmp;
               }  else {
                 q=q->next;
               }
             }         //End of while(q)
           }           //End of qHead->next

           q=qHead->next;
           head=D->particle[i][j][k].head[s]->pt;
           while(q)  {
             storeParticle(D,i,j,k,s,q->x,0.0,0.0,q->weight);
             qHead->next=q->next;
             q->next=NULL;
             free(q);
             q=qHead->next;
           }

         }  else; 	//End of W>0
       }	//End of i,j,k     

   free(qHead);
}

void storeParticle(Domain *D,int i,int j,int k,int s,double x,double y,double z,double weight)
{
  ptclList *new;
 
  new = (ptclList *)malloc(sizeof(ptclList ));
  new->next = D->particle[i][j][k].head[s]->pt;
  D->particle[i][j][k].head[s]->pt=new;

  new->x=x; new->y=y; new->z=z; new->weight=weight;
}

