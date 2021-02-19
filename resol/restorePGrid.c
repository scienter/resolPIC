#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);


void restoreQGrid(Domain *D,int s)
{
   int i,nx,ny,nz,nxSub,nySub,nzSub,biasX,biasY,biasZ,offset[3];
   int istart,iend,jstart,jend,kstart,kend,rankX,rankY,rankZ;
   char fileName[100],dataName[100];
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   //restore meta data
   sprintf(fileName,"PGrid%d.h5",D->step);

   rankX=myrank/(D->M*D->N);
   rankZ=(myrank%(D->M*D->N))/D->M;
   rankY=(myrank%(D->M*D->N))%D->M;

   switch (D->dimension) {
   case 1 :
     nx=D->nx+5; nxSub=D->nxSub+5; istart=0; iend=D->iend+3;
     ny=1;       nySub=1;          jstart=0; jend=1;
     nz=1;       nzSub=1;          kstart=0; kend=1;
     offset[0]=D->minXSub; offset[1]=0; offset[2]=0;

     sprintf(dataName,"%d/q1",s);
     restoreFieldComp(D->q1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/q2",s);
     restoreFieldComp(D->q2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/px1",s);
     restoreFieldComp(D->px1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/px2",s);
     restoreFieldComp(D->px2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/py1",s);
     restoreFieldComp(D->py1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/py2",s);
     restoreFieldComp(D->py2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/pz1",s);
     restoreFieldComp(D->pz1,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/pz2",s);
     restoreFieldComp(D->pz2,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/sigX",s);
     restoreFieldComp(D->sigX,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/sigY",s);
     restoreFieldComp(D->sigY,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/sigZ",s);
     restoreFieldComp(D->sigZ,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     sprintf(dataName,"%d/energy",s);
     restoreFieldComp(D->energy,fileName,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);
     break;  
   }
}

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L)
{
  if(L==1) {
    *istart=0; *iend+=3; *biasX=0; *nxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0; *nxSub+=2; *biasX=0;
    }  else if(rankX==L-1)  {
      *iend+=3;  *nxSub+=3; *biasX=2;
    } else       *biasX=2;
  }
}
