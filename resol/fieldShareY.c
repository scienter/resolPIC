#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"


void MPI_Transfer_Yminus(Domain *D,double ***f1,int nx,int nz,int iend,int jend,int kend,int share)
{
    int i,j,k,numberData,start,end,num;
    int istart,jstart,kstart;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=jstart=kstart=2;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))%D->M;

    num=nx*nz*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)	
        for(k=0; k<nz; k++) data[start+k]=f1[i][j+jstart][k]; start+=nz;

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
          for(k=0; k<nz; k++) f1[i][jend+j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)	
      for(j=0; j<share; j++)
        for(k=0; k<nz; k++) data[start+k]=f1[i][jstart+j][k]; start+=nz;
        
    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)	
          for(k=0; k<nz; k++) f1[i][jend+j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


void MPI_Transfer_Yplus(Domain *D,double ***f1,int nx,int nz,int iend,int jend,int kend,int share)
{
    int i,j,k,num;
    int istart,jstart,kstart;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
    istart=jstart=kstart=2;

    rank=(myrank%(D->M*D->N))%D->M;
    num=nx*nz*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
        for(k=0; k<nz; k++) data[start+k]=f1[i][jend-j][k]; start+=nz;
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
          for(k=0; k<nz; k++) f1[i][jstart-j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)	
        for(k=0; k<nz; k++) data[start+k]=f1[i][jend-j][k]; start+=nz;
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)	
          for(k=0; k<nz; k++) f1[i][jstart-j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==1 && rank!=D->M-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

