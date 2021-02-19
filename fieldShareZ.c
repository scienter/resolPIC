#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void MPI_Transfer6F_Zminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;

    num=6*nx*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f4[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f5[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f6[i][j][k+kstart]; start+=share;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++) f1[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f4[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f5[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f6[i][j][k+kend]=data[start+k]; start+=share;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f4[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f5[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f6[i][j][k+kstart]; start+=share;
      }

    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++) f1[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f4[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f5[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f6[i][j][k+kend]=data[start+k]; start+=share;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Zplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int nx,int ny,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    rank=(myrank%(D->M*D->N))/D->M;
    num=6*nx*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f4[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f5[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f6[i][j][kend-k]; start+=2;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++) f1[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f4[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f5[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f6[i][j][kstart-k]=data[start+k-1]; start+=2;
        }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~even ~ cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f4[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f5[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f6[i][j][kend-k]; start+=2;
      }
      
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++) f1[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f4[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f5[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f6[i][j][kstart-k]=data[start+k-1]; start+=2;
        }
    }
    else if(rank%2==1 && rank!=D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer3F_Zminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,numberData,start,end,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;

    num=3*nx*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][k+kstart]; start+=share;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++) f1[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][k+kend]=data[start+k]; start+=share;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][k+kstart]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][k+kstart]; start+=share;
      }

    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++) f1[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][k+kend]=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][k+kend]=data[start+k]; start+=share;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer3F_Zplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            
   
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    rank=(myrank%(D->M*D->N))/D->M;
    num=3*nx*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kend-k]; start+=2;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++) f1[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kstart-k]=data[start+k-1]; start+=2;
        }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~even ~ cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kend-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kend-k]; start+=2;
      }
      
    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++) f1[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kstart-k]=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kstart-k]=data[start+k-1]; start+=2;
        }
    }
    else if(rank%2==1 && rank!=D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_TransferJ_Zminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,start,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;
    num=3*nx*ny*2;
    data = (double *)malloc(num*sizeof(double ));


    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kstart-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kstart-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kstart-k]; start+=2;
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	 {
          for(k=1; k<share; k++) f1[i][j][kend-k]+=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kend-k]+=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kend-k]+=data[start+k-1]; start+=2;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=1; k<share; k++) data[start+k-1]=f1[i][j][kstart-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f2[i][j][kstart-k]; start+=2;
        for(k=1; k<share; k++) data[start+k-1]=f3[i][j][kstart-k]; start+=2;
      }
        
    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	 {
          for(k=1; k<share; k++) f1[i][j][kend-k]+=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f2[i][j][kend-k]+=data[start+k-1]; start+=2;
          for(k=1; k<share; k++) f3[i][j][kend-k]+=data[start+k-1]; start+=2;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_TransferJ_Zplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int ny,int share)
{
    int i,j,k,start,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;
    num=3*nx*ny*share;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++) {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][kend+k]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][kend+k]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][kend+k]; start+=share;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	{
          for(k=0; k<share; k++) f1[i][j][kstart+k]+=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][kstart+k]+=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][kstart+k]+=data[start+k]; start+=share;
        }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++) {
        for(k=0; k<share; k++) data[start+k]=f1[i][j][kend+k]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f2[i][j][kend+k]; start+=share;
        for(k=0; k<share; k++) data[start+k]=f3[i][j][kend+k]; start+=share;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevZrank,D->prevZrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	{
          for(k=0; k<share; k++) f1[i][j][kstart+k]+=data[start+k]; start+=share;
          for(k=0; k<share; k++) f2[i][j][kstart+k]+=data[start+k]; start+=share;
          for(k=0; k<share; k++) f3[i][j][kstart+k]+=data[start+k]; start+=share;
        }
    }
    else if(rank%2==1 && rank!=D->N-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);

}


void MPI_TransferRho_Zminus(Domain *D,double ***f1,int nx,int ny,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank;
    double *Zminus; 

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;
    numShare=nx*ny*(share-1);
    Zminus = (double *)malloc(numShare*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=1; k<share; k++)
          Zminus[start+k-1]=f1[i][j][kstart-k];
        start+=(share-1);
      }

    if(rank%2==0 && rank!=D->N-1)
    {
      MPI_Recv(Zminus,numShare,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++)
            f1[i][j][kend-k]+=Zminus[start+k-1];
          start+=(share-1);
        }  
    }
    else if(rank%2==1)
       MPI_Send(Zminus,numShare,MPI_DOUBLE,D->prevZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=1; k<share; k++)
          Zminus[start+k-1]=f1[i][j][kstart-k];
        start+=(share-1);
      }
        
    if(rank%2==1 && rank!=D->N-1)
    {
      MPI_Recv(Zminus,numShare,MPI_DOUBLE,D->nextZrank,D->nextZrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=1; k<share; k++)
            f1[i][j][kend-k]+=Zminus[start+k-1];
          start+=(share-1);
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(Zminus,numShare,MPI_DOUBLE,D->prevZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    free(Zminus);
}

void MPI_TransferRho_Zplus(Domain *D,double ***f1,int nx,int ny,int share)
{
    int i,j,k,start,numShare;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank; 
    double *Zplus;

    MPI_Status status;         
   
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=(myrank%(D->M*D->N))/D->M;
    numShare=nx*ny*share;
    Zplus = (double *)malloc(numShare*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++)
          Zplus[start+k]=f1[i][j][kend+k];
        start+=share;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(Zplus,numShare,MPI_DOUBLE,D->prevZrank,D->prevZrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++)
            f1[i][j][kstart+k]+=Zplus[start+k];
          start+=share;
        }
    }
    else if(rank%2==0 && rank!=D->N-1)
       MPI_Send(Zplus,numShare,MPI_DOUBLE,D->nextZrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)	
      {
        for(k=0; k<share; k++)
          Zplus[start+k]=f1[i][j][kend+k];
        start+=share;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(Zplus,numShare,MPI_DOUBLE,D->prevZrank,D->prevZrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)	
        {
          for(k=0; k<share; k++)
            f1[i][j][kstart+k]+=Zplus[start+k];
          start+=share;
        }
    }
    else if(rank%2==1 && rank!=D->N-1)
      MPI_Send(Zplus,numShare,MPI_DOUBLE,D->nextZrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    free(Zplus);
}

