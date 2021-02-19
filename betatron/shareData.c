#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "parameter.h"
#include <mpi.h>

void shareSpectrumData(Parameter *D)
{
   int i,j,start,numberData;
   int myrank,nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   numberData=6;
   for(i=1; i<nTasks; i++)
   {
     if(myrank==i)
     {
       for(i=0; i<3; i++)
         D->shareBC[i]=D->B[i];
       for(i=0; i<3; i++)
         D->shareBC[3+i]=D->C[i];

       MPI_Send(D->shareBC,numberData,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
     else if(myrank==0)
     {
       MPI_Recv(D->shareBC,numberData,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
       for(i=0; i<3; i++)
         D->B[i]+=D->shareBC[i];
       for(i=0; i<3; i++)
         D->C[i]+=D->shareBC[3+i];
     }
     else	; 
   }
}

void shareData(Parameter *D)
{
   int i,j,start,numberData;
   int myrank,nTasks;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   numberData=D->numThX*D->numThY;
   for(i=1; i<nTasks; i++)
   {
     if(myrank==i)
     {
       start=0;
       for(i=0; i<D->numThX; i++)
       {
         for(j=0; j<D->numThY; j++)
           D->share[start+j]=D->data[i][j];
         start+=D->numThY;
       }
       MPI_Send(D->share,numberData,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     }
     else if(myrank==0)
     {
       MPI_Recv(D->share,numberData,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<D->numThX; i++)
       {
         for(j=0; j<D->numThY; j++)
           D->data[i][j]+=D->share[start+j];
         start+=D->numThY;
       }
     }
     else	; 
   }
}

