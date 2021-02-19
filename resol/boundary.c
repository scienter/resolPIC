#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string.h>


double ***memoryAsign(int nx, int ny, int nz);
void calParameter(int nx,int *istart,int *iend,int *saveNxSub,int rankX,int *biasX,int reIend,int reNxSub,int L);
void HighBoundary(Domain *D);
void LowBoundary(Domain *D);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void boundary(Domain *D)
{
  if(D->mode==HIGH) HighBoundary(D);
  else if(D->mode==LOW) LowBoundary(D);
  else ;
}


void HighBoundary(Domain *D)
{
   int i,j,k,s,nx,ny,nz;
   int subX,subY,subZ,remainX,remainY,remainZ,nxSub1D,nySub2D,nzSub3D;
   int minX,minY,minZ,maxX,maxY,maxZ,rank,rankX,rankY,rankZ,tmpX,tmpY,tmpZ;
   char fileName[100],dataName[100];
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //restore meta data
   sprintf(fileName,"resolE%d.h5",D->step);
//   sprintf(fileName,"dumpResolParticle%d.h5",D->step+1);
   if(myrank==0)  {
     restoreIntMeta(fileName,"/nx",&D->nx,1);
     restoreIntMeta(fileName,"/ny",&D->ny,1);
     restoreIntMeta(fileName,"/nz",&D->nz,1);
     restoreIntMeta(fileName,"/minXDomain",&D->minXDomain,1);
     restoreIntMeta(fileName,"/minYDomain",&D->minYDomain,1);
     restoreIntMeta(fileName,"/minZDomain",&D->minZDomain,1);
   }    else   ;
   MPI_Bcast(&D->nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minYDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minZDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   sprintf(fileName,"dumpResolParticle%d.h5",D->step);
   if(myrank==0)  {
     restoreIntMeta(fileName,"/nSpecies",&D->nSpecies,1);
   }    else   ;
   MPI_Bcast(&D->nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
//   if(D->minXDomain>0)  D->minXDomain-=1;  else ;

   if(D->nx%D->resolX!=0 || D->ny%D->resolY!=0 || D->nz%D->resolZ!=0)  {
       printf("check nx, ny, nz. The remain of resol value is not zero.\n");
       exit(0);
     }  else ;
   nx=D->nx;  ny=D->ny;  nz=D->nz;

   D->nxSub=nx/D->L; subX=D->nxSub; remainX=nx%D->L; minX=maxX=0;
   D->nySub=ny/D->M; subY=D->nySub; remainY=ny%D->M; minY=maxY=0;
   D->nzSub=nz/D->N; subZ=D->nzSub; remainZ=nz%D->N; minZ=maxZ=0;

   minX=maxX=0;
   for(rankX=0; rankX<D->L; rankX++)   {
     if(rankX<remainX)   tmpX=subX+1;
     else                tmpX=subX;
     minX=maxX;
     maxX=minX+tmpX;

     minZ=maxZ=0;
     for(rankZ=0; rankZ<D->N; rankZ++)     {
       if(rankZ<remainZ)   tmpZ=subZ+1;
       else                tmpZ=subZ;
       minZ=maxZ;
       maxZ=minZ+tmpZ;

       minY=maxY=0;
       for(rankY=0; rankY<D->M; rankY++)       {
         if(rankY<remainY)   tmpY=subY+1;
         else                tmpY=subY;
         minY=maxY;
         maxY=minY+tmpY;

         rank=rankY+rankZ*D->M+rankX*(D->M*D->N);
         if(myrank==rank)         {
           D->nxSub=tmpX;
           D->minXSub=minX;
           D->maxXSub=maxX;
           D->nySub=tmpY;
           D->minYSub=minY;
           D->maxYSub=maxY;
           D->nzSub=tmpZ;
           D->minZSub=minZ;
           D->maxZSub=maxZ;
         }
       }
     }
   }
//   D->minXSub+=D->minXDomain;
   D->minYSub+=D->minYDomain;
   D->minZSub+=D->minZDomain;
//   D->maxXSub+=D->minXDomain;
   D->maxYSub+=D->minYDomain;
   D->maxZSub+=D->minZDomain;

   D->istart=2; D->iend=D->nxSub+2;
   D->jstart=0; D->jend=1;
   D->kstart=0; D->kend=1;
   D->reIend=D->nxSub*D->resolX+2;
   D->reJend=1;
   D->reKend=1;
   nxSub1D=D->nxSub*D->resolX+5; nySub2D=1; nzSub3D=1;
   if(D->dimension>1) {
     nySub2D=D->nySub*D->resolY+5; D->jstart=2; D->jend=D->nySub+2;
     D->reJend=D->nySub*D->resolY+2;
   }  else;
   if(D->dimension>2) {
     nzSub3D=D->nzSub*D->resolZ+5; D->kstart=2; D->kend=D->nzSub+2;
     D->reKend=D->nzSub*D->resolZ+2;
   }  else;
   D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
   for(i=0; i<nxSub1D; i++) {
     D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
     for(j=0; j<nySub2D; j++)  {
       D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
     }
   }

   switch (D->dimension) {
   case 1:
     j=k=0;
     for(i=D->istart; i<D->reIend; i++)     {
       D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
       for(s=0; s<D->nSpecies; s++)      {
         D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
         D->particle[i][j][k].head[s]->pt = NULL;
       }
     }
     break;
   case 2:
     k=0;
     for(i=D->istart; i<D->reIend; i++)
       for(j=D->jstart; j<D->reJend; j++)     {
         D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)      {
           D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j][k].head[s]->pt = NULL;
         }
       }
     break;
   }

   //Field restoring
   rankX=myrank/(D->M*D->N);
   rankZ=(myrank%(D->M*D->N))/D->M;
   rankY=(myrank%(D->M*D->N))%D->M;

   D->nextXrank=rankY+rankZ*D->M+(rankX+1)*(D->M*D->N);
   D->prevXrank=rankY+rankZ*D->M+(rankX-1)*(D->M*D->N);
   if(rankX==D->L-1)
     D->nextXrank=rankY+rankZ*D->M;
   else if(rankX==0)
     D->prevXrank=rankY+rankZ*D->M+(D->L-1)*(D->M*D->N);

   D->nextYrank=(rankY+1)+rankZ*D->M+rankX*(D->M*D->N);
   D->prevYrank=(rankY-1)+rankZ*D->M+rankX*(D->M*D->N);
   if(rankY==D->M-1)
     D->nextYrank=rankZ*D->M+rankX*(D->M*D->N);
   else if(rankY==0)
     D->prevYrank=(D->M-1)+rankZ*D->M+rankX*(D->M*D->N);

   D->nextZrank=rankY+(rankZ+1)*D->M+rankX*(D->M*D->N);
   D->prevZrank=rankY+(rankZ-1)*D->M+rankX*(D->M*D->N);
   if(rankZ==D->N-1)
     D->nextZrank=rankY+rankX*(D->M*D->N);
   else if(rankZ==0)
     D->prevZrank=rankY+(D->N-1)*D->M+rankX*(D->M*D->N);




   //setting resolution change
   D->reNx=D->nx*D->resolX;
   D->reNy=D->ny*D->resolY;
   D->reNz=D->nz*D->resolZ;

   D->reNxSub=D->nxSub*D->resolX;
   D->reNySub=D->nySub*D->resolY;
   D->reNzSub=D->nzSub*D->resolZ;

   calParameter(D->reNx+5,&D->saveIstart,&D->saveIend,&D->saveNxSub,rankX,&D->biasX,D->reIend,D->reNxSub,D->L);
   D->saveJstart=0; D->saveJend=1;
   D->saveKstart=0; D->saveKend=1;
   D->saveNySub=1;  D->saveNzSub=1;
   calParameter(D->reNx+5,&D->saveIstart,&D->saveIend,&D->saveNxSub,rankX,&D->biasX,D->reIend,D->reNxSub,D->L);
   D->saveJstart=0; D->saveJend=1;
   D->saveKstart=0; D->saveKend=1;
   D->saveNySub=1;  D->saveNzSub=1;
   if(D->dimension>1)
     calParameter(D->reNy+5,&D->saveJstart,&D->saveJend,&D->saveNySub,D->rankY,&D->biasY,D->reJend,D->reNySub,D->M);
   else ;
   if(D->dimension>2)
     calParameter(D->reNz+5,&D->saveKstart,&D->saveKend,&D->saveNzSub,D->rankZ,&D->biasZ,D->reKend,D->reNzSub,D->N);
   else ;

   nxSub1D=D->nxSub+5; nySub2D=1; nzSub3D=1;
   if(D->dimension>1) nySub2D=D->nySub+5; else;
   if(D->dimension>2) nzSub3D=D->nzSub+5; else;
   D->fieldE=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldNext=memoryAsign(nxSub1D,nySub2D,nzSub3D);

}


void LowBoundary(Domain *D)
{
   int i,j,k,s,nx,ny,nz;
   int subX,subY,subZ,remainX,remainY,remainZ,nxSub1D,nySub2D,nzSub3D;
   int minX,minY,minZ,maxX,maxY,maxZ,rank,rankX,rankY,rankZ,tmpX,tmpY,tmpZ;
   char fileName[100],dataName[100];
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //restore meta data
   sprintf(fileName,"resolE%d.h5",D->step);
//   sprintf(fileName,"dumpResolParticle%d.h5",D->step+1);
   if(myrank==0)  {
     restoreIntMeta(fileName,"/nx",&D->nx,1);
     restoreIntMeta(fileName,"/ny",&D->ny,1);
     restoreIntMeta(fileName,"/nz",&D->nz,1);
     restoreIntMeta(fileName,"/minXDomain",&D->minXDomain,1);
     restoreIntMeta(fileName,"/minYDomain",&D->minYDomain,1);
     restoreIntMeta(fileName,"/minZDomain",&D->minZDomain,1);
   }    else   ;
   MPI_Bcast(&D->nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minYDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&D->minZDomain,1,MPI_INT,0,MPI_COMM_WORLD);
   sprintf(fileName,"dumpResolParticle%d.h5",D->step+1);
   if(myrank==0)  {
     restoreIntMeta(fileName,"/nSpecies",&D->nSpecies,1);
   }    else   ;
   MPI_Bcast(&D->nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
//   if(D->minXDomain>0)  D->minXDomain-=1;  else ;

   if(D->nx%D->resolX!=0 || D->ny%D->resolY!=0 || D->nz%D->resolZ!=0)  {
       printf("check nx, ny, nz. The remain of resol value is not zero.\n");
       exit(0);
     }  else ;
   nx=D->nx/D->resolX;  ny=D->ny/D->resolY;  nz=D->nz/D->resolZ;

   D->nxSub=nx/D->L; subX=D->nxSub; remainX=nx%D->L; minX=maxX=0;
   D->nySub=ny/D->M; subY=D->nySub; remainY=ny%D->M; minY=maxY=0;
   D->nzSub=nz/D->N; subZ=D->nzSub; remainZ=nz%D->N; minZ=maxZ=0;

   minX=maxX=0;
   for(rankX=0; rankX<D->L; rankX++)   {
     if(rankX<remainX)   tmpX=subX+1;
     else                tmpX=subX;
     minX=maxX;
     maxX=minX+tmpX;

     minZ=maxZ=0;
     for(rankZ=0; rankZ<D->N; rankZ++)     {
       if(rankZ<remainZ)   tmpZ=subZ+1;
       else                tmpZ=subZ;
       minZ=maxZ;
       maxZ=minZ+tmpZ;

       minY=maxY=0;
       for(rankY=0; rankY<D->M; rankY++)       {
         if(rankY<remainY)   tmpY=subY+1;
         else                tmpY=subY;
         minY=maxY;
         maxY=minY+tmpY;

         rank=rankY+rankZ*D->M+rankX*(D->M*D->N);
         if(myrank==rank)         {
           D->nxSub=tmpX*D->resolX;
           D->minXSub=minX*D->resolX;
           D->maxXSub=maxX*D->resolX;
           D->nySub=tmpY*D->resolY;
           D->minYSub=minY*D->resolY;
           D->maxYSub=maxY*D->resolY;
           D->nzSub=tmpZ*D->resolZ;
           D->minZSub=minZ*D->resolZ;
           D->maxZSub=maxZ*D->resolZ;
         }
       }
     }
   }
//   D->minXSub+=D->minXDomain;
   D->minYSub+=D->minYDomain;
   D->minZSub+=D->minZDomain;
//   D->maxXSub+=D->minXDomain;
   D->maxYSub+=D->minYDomain;
   D->maxZSub+=D->minZDomain;

   D->istart=2; D->iend=D->nxSub+2;
   D->jstart=0; D->jend=1;
   D->kstart=0; D->kend=1;
   D->reIend=D->nxSub/D->resolX+2;
   D->reJend=1;
   D->reKend=1;
   nxSub1D=D->nxSub/D->resolX+5; nySub2D=1; nzSub3D=1;
   if(D->dimension>1) {
     nySub2D=D->nySub/D->resolY+5; D->jstart=2; D->jend=D->nySub+2;
     D->reJend=D->nySub/D->resolY+2;
   }  else;
   if(D->dimension>2) {
     nzSub3D=D->nzSub/D->resolZ+5; D->kstart=2; D->kend=D->nzSub+2;
     D->reKend=D->nzSub/D->resolZ+2;
   }  else;
   D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
   for(i=0; i<nxSub1D; i++) {
     D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
     for(j=0; j<nySub2D; j++)  {
       D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
     }
   }

   switch (D->dimension) {
   case 1:
     j=k=0;
     for(i=D->istart; i<D->reIend; i++)     {
       D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
       for(s=0; s<D->nSpecies; s++)      {
         D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
         D->particle[i][j][k].head[s]->pt = NULL;
       }
     }
     break;
   case 2:
     k=0;
     for(i=D->istart; i<D->reIend; i++)
       for(j=D->jstart; j<D->reJend; j++)     {
         D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)      {
           D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j][k].head[s]->pt = NULL;
         }
       }
     break;
   }

   //Field restoring
   D->rankX=myrank/(D->M*D->N);
   D->rankZ=(myrank%(D->M*D->N))/D->M;
   D->rankY=(myrank%(D->M*D->N))%D->M;

   D->nextXrank=rankY+rankZ*D->M+(rankX+1)*(D->M*D->N);
   D->prevXrank=rankY+rankZ*D->M+(rankX-1)*(D->M*D->N);
   if(rankX==D->L-1)
     D->nextXrank=rankY+rankZ*D->M;
   else if(rankX==0)
     D->prevXrank=rankY+rankZ*D->M+(D->L-1)*(D->M*D->N);

   D->nextYrank=(rankY+1)+rankZ*D->M+rankX*(D->M*D->N);
   D->prevYrank=(rankY-1)+rankZ*D->M+rankX*(D->M*D->N);
   if(rankY==D->M-1)
     D->nextYrank=rankZ*D->M+rankX*(D->M*D->N);
   else if(rankY==0)
     D->prevYrank=(D->M-1)+rankZ*D->M+rankX*(D->M*D->N);

   D->nextZrank=rankY+(rankZ+1)*D->M+rankX*(D->M*D->N);
   D->prevZrank=rankY+(rankZ-1)*D->M+rankX*(D->M*D->N);
   if(rankZ==D->N-1)
     D->nextZrank=rankY+rankX*(D->M*D->N);
   else if(rankZ==0)
     D->prevZrank=rankY+(D->N-1)*D->M+rankX*(D->M*D->N);

   //setting resolution change
   D->reNx=nx;
   D->reNy=ny;
   D->reNz=nz;

   D->reNxSub=D->nxSub/D->resolX;
   D->reNySub=D->nySub/D->resolY;
   D->reNzSub=D->nzSub/D->resolZ;

   calParameter(D->reNx+5,&D->saveIstart,&D->saveIend,&D->saveNxSub,D->rankX,&D->biasX,D->reIend,D->reNxSub,D->L);
   D->saveJstart=0; D->saveJend=1;
   D->saveKstart=0; D->saveKend=1;
   D->saveNySub=1;  D->saveNzSub=1;
   if(D->dimension>1)
     calParameter(D->reNy+5,&D->saveJstart,&D->saveJend,&D->saveNySub,D->rankY,&D->biasY,D->reJend,D->reNySub,D->M);
   else ;
   if(D->dimension>2)
     calParameter(D->reNz+5,&D->saveKstart,&D->saveKend,&D->saveNzSub,D->rankZ,&D->biasZ,D->reKend,D->reNzSub,D->N);
   else ;

   nxSub1D=D->nxSub+5; nySub2D=1; nzSub3D=1;
   if(D->dimension>1) nySub2D=D->nySub+5; else;
   if(D->dimension>2) nzSub3D=D->nzSub+5; else;
   D->fieldE=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldOld=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldNow=memoryAsign(nxSub1D,nySub2D,nzSub3D);
   D->fieldNext=memoryAsign(nxSub1D,nySub2D,nzSub3D);

}

void calParameter(int nx,int *istart,int *iend,int *saveNxSub,int rankX,int *biasX,int reIend,int reNxSub,int L)
{
  if(L==1) {
    *istart=0;
    *iend=reIend+3;
    *biasX=0;
    *saveNxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0;
      *iend=reIend;
      *saveNxSub=reNxSub+2;
      *biasX=0;
    }  else if(rankX==L-1)  {
      *istart=2;
      *iend=reIend+3;
      *saveNxSub=reNxSub+3;
      *biasX=2;
    } else  {
      *istart=2;
      *iend=reIend;
      *saveNxSub=reNxSub;
      *biasX=2;
    }
  }
}


double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)
   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0.0;
       }
   return field;
}
  
int ***memoryAsignInt(int nx, int ny, int nz)
{
   int i,j,k;
   int ***field;

   field = (int ***)malloc((nx)*sizeof(int **));
   for(i=0; i<nx; i++)
   {
     field[i] = (int **)malloc((ny)*sizeof(int *));
     for(j=0; j<ny; j++)
       field[i][j] = (int *)malloc((nz)*sizeof(int ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0;
       }
   return field;
}

int whatMode(char *str)
{
   if(strstr(str,"high"))           return HIGH;
   else if(strstr(str,"low"))       return LOW;
   else return 0;
}


