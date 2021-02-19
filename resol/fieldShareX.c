#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_Transfer_Xminus(Domain *D,double ***f1,int ny,int nz,int iend,int jend,int kend,int share)
{
    int i,j,k,num,start,end;
    int istart,jstart,kstart;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=jstart=kstart=2;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=ny*nz*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=0; i<share; i++)
      for(j=0; j<ny; j++)	 
        for(k=0; k<nz; k++) data[start+k]=f1[i+istart][j][k]; start+=nz;
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++)	
          for(k=0; k<nz; k++) f1[iend+i][j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


void MPI_Transfer_Xplus(Domain *D,double ***f1,int ny,int nz,int iend,int jend,int kend,int share)
{
    int i,j,k,num;
    int istart,jstart,kstart;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=jstart=kstart=2;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/(D->M*D->N);

    num=ny*nz*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)	
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=1; i<share; i++)
      for(j=0; j<ny; j++)      
        for(k=0; k<nz; k++) data[start+k]=f1[iend-i][j][k]; start+=nz;
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)	
          for(k=0; k<nz; k++) f1[istart-i][j][k]=data[start+k]; start+=nz;
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

