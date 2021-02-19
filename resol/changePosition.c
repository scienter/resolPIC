#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void changeParticlePosition(Domain *D)
{
   int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
   double resolX,resolY,resolZ;
   ptclList *p;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->reIend;
   jstart=D->jstart;    jend=D->reJend;
   kstart=D->kstart;    kend=D->reKend;

   if(D->mode==LOW) {
     resolX=1.0/((double)D->resolX);
     resolY=1.0/((double)D->resolY);
     resolZ=1.0/((double)D->resolZ);
   } else if(D->mode==HIGH) {
     resolX=(double)D->resolX;
     resolY=(double)D->resolY;
     resolZ=(double)D->resolZ;
   } else ;


   for(i=istart; i<iend; i++) 
   {
     for(j=jstart; j<jend; j++) 
       for(k=kstart; k<kend; k++)
         for(s=0; s<D->nSpecies; s++)  {
           p=D->particle[i][j][k].head[s]->pt;
           while(p)  {
             p->x+=i-istart+D->minXSub*resolX;
             p->y+=j-jstart+D->minYSub*resolY;
             p->z=0.0;
             p->weight*=resolX*resolY*resolZ;
             p=p->next;
           }
         }
     if(i%10==0 && myrank==0) 
       printf("%d/%d is done\n",i,iend-istart);
     else;
   }

}

