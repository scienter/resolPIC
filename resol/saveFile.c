#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void saveParticle(Domain *D,int s,Particle ***particle)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
   double x,y,z,px,py,pz,weight,charge;
   char fileName[100];
   FILE *out;
   ptclList *p;

   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart; iend=D->reIend;
   jstart=D->jstart; jend=D->reJend;
   kstart=D->kstart; kend=D->reKend;
 
   sprintf(fileName,"reParticle%d_%d",D->step,myrank); 
   out=fopen(fileName,"w");
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++) 
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)  {
           x=i+p->x-istart+D->minXSub;  y=p->y;  z=p->z;
           px=p->px;  py=p->py;  pz=p->pz;
           weight=p->weight;
           charge=p->charge;
           fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,z,px,py,pz,weight,charge);
           p=p->next;
         }
       }	//End of i,j,k     
   fclose(out);
}

