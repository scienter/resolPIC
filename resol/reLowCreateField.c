#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void resolLowCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int edge);
void resolLowCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int edge);
void calLowBField(double ***fieldOld2,double ***fieldOld1,double ***fieldNow,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int *minXDomain);
void calLowJField(double ***fieldOld2,double ***fieldOld1,double ***fieldNow,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int *minXDomain);
void reLowCreateField1D(Domain *D);
void reLowCreateField2D(Domain *D);
double ***memoryAsign(int nx, int ny, int nz);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void deleteField(double ***field,int nx,int ny,int nz);
void MPI_Transfer_Xminus(Domain *D,double ***f1,int ny,int nz,int iend,int jend,int kend,int share);
void MPI_Transfer_Xplus(Domain *D,double ***f1,int ny,int nz,int iend,int jend,int kend,int share);
void MPI_Transfer_Yminus(Domain *D,double ***f1,int nx,int nz,int iend,int jend,int kend,int share);
void MPI_Transfer_Yplus(Domain *D,double ***f1,int nx,int nz,int iend,int jend,int kend,int share);


void reLowCreateField(Domain D)
{
   int istart,iend,jstart,jend,kstart,kend;

   istart=D.istart;    iend=D.iend;
   jstart=D.jstart;    jend=D.jend;
   kstart=D.kstart;    kend=D.kend;

   switch (D.dimension)  {
   case 1 :
     reLowCreateField1D(&D);
     break;
   case 2 :
     reLowCreateField2D(&D);
     break;
   }

}

void reLowCreateField2D(Domain *D)
{
   int nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend;
   int reIend,reJend,reKend;
   int offset[3],offSet[3],tmp,nxSub1D,nySub2D,nzSub3D,minXDomain[3];
   double ***reField1,***reField2,resolX;
   char fileName[100],outFile[100];
   hid_t file_id;
   
   int myrank, nTasks,shift=0;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   nx=D->nx; ny=D->ny; nz=D->nz;
   nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
   istart=D->istart; iend=D->iend; reIend=D->reIend;
   jstart=D->jstart; jend=D->jend; reJend=D->reJend;
   kstart=D->kstart; kend=D->kend; reKend=D->reKend;

   offset[0]=D->minXSub;
   offset[1]=D->minYSub-D->minYDomain;
   offset[2]=0;

   //for saving field file
   offSet[0]=(D->minXSub)/D->resolX+D->biasX;
   offSet[1]=(D->minYSub-D->minYDomain)/D->resolY+D->biasY;
   offSet[2]=0;

   resolX=1.0/(double)D->resolX;

   //create out file
   sprintf(outFile,"dumpField%g.h5",D->step*resolX);
   if(myrank==0)
   {
     file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
     H5Fclose(file_id);
   }
   else ;

   reField1=memoryAsign(D->reNxSub+5,D->nySub+5,1);
   reField2=memoryAsign(D->reNxSub+5,D->reNySub+5,1);

   //save E field
   sprintf(fileName,"resolE%d.h5",D->step);
   restoreFieldComp(D->fieldE,fileName,"/Ex",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   resolLowCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Ex",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ey",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   resolLowCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Ey",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ez",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   resolLowCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Ez",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);
   
   //save B field
   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   if(myrank==0)
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[0],1);
   else;
   MPI_Bcast(&minXDomain[0],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldOld,fileName,"/Bx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   if(myrank==0)
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[1],1);
   else;
   MPI_Bcast(&minXDomain[1],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldNow,fileName,"/Bx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   if(myrank==0)
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[2],1);
   else;
   MPI_Bcast(&minXDomain[2],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldNext,fileName,"/Bx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Bx",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/By",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Bz",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

    //calculation current
   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Jx",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolLowCalY(reField1,reField2,1,reIend+1,1,jend+1,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Jy",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
//   if(D->M>1)  {
//     MPI_Transfer_Yminus(D,D->fieldNow,nxSub+5,1,iend,jend,1,3);
//     MPI_Transfer_Yplus(D,D->fieldNow,nxSub+5,1,iend,jend,1,3);
//   } else      ;
   resolLowCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolLowCalY(reField1,reField2,2,reIend,2,jend,0,1,D->resolY,1);
//   if(D->L>1)  {
//     MPI_Transfer_Xminus(D,reField1,D->nySub+5,1,reIend,jend,1,3);
//     MPI_Transfer_Xplus(D,reField1,D->nySub+5,1,reIend,jend,1,3);
//   } else      ;
   saveFieldComp(reField2,outFile,"/Jz",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   nxSub1D=D->nxSub+5; nySub2D=1; nzSub3D=1;
   if(D->dimension>1) nySub2D=D->nySub+5;  else;
   if(D->dimension>2) nzSub3D=D->nzSub+5;  else;
   deleteField(reField1,D->reNxSub+5,D->nySub+5,1);
   deleteField(reField2,D->reNxSub+5,D->reNySub+5,1);
   deleteField(D->fieldE,nxSub1D,nySub2D,nzSub3D);
   deleteField(D->fieldOld,nxSub1D,nySub2D,nzSub3D);
   deleteField(D->fieldNow,nxSub1D,nySub2D,nzSub3D);
   deleteField(D->fieldNext,nxSub1D,nySub2D,nzSub3D);

   D->nx=D->nx/D->resolX;
   D->ny=D->ny/D->resolY; 
   D->nz=D->nz/D->resolZ; 
   D->minXDomain/=D->resolX;
   D->minYDomain/=D->resolY;
   D->minZDomain/=D->resolZ;

   if(myrank==0) {
     saveIntMeta(outFile,"minXDomain",&D->minXDomain,1);
     saveIntMeta(outFile,"minYDomain",&D->minYDomain,1);
     saveIntMeta(outFile,"minZDomain",&D->minZDomain,1);
     saveIntMeta(outFile,"nSpecies",&D->nSpecies,1);
     saveIntMeta(outFile,"nx",&D->nx,1);
     saveIntMeta(outFile,"ny",&D->ny,1);
     saveIntMeta(outFile,"nz",&D->nz,1);
     printf("%s is made\n",outFile);
   } else;
}


void reLowCreateField1D(Domain *D)
{
   int nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend;
   int offset[3],offSet[3],tmp,minXDomain[3];
   double ***reField;
   char fileName[100],outFile[100];
   hid_t file_id;
   
   int myrank, nTasks,shift=1;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   nx=D->nx; ny=D->ny; nz=D->nz;
   nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   kstart=D->kstart; kend=D->kend;

   offset[0]=D->minXSub;
//   offset[1]=D->minYSub-D->minYDomain;
   offset[1]=0;
   offset[2]=0;

   //for saving field file
   offSet[0]=(D->minXSub)/D->resolX+D->biasX;
//   offSet[1]=(D->minYSub-D->minYDomain)*D->resolY+D->biasY;
   offSet[1]=0;
   offSet[2]=0;

   reField=memoryAsign(D->reNxSub+5,1,1);

   //create out file
   sprintf(outFile,"dumpField%d.h5",D->step/D->resolX);
   if(myrank==0)
   {
     file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
     H5Fclose(file_id);
   }
   else ;

   //save E field
   sprintf(fileName,"resolE%d.h5",D->step);
   restoreFieldComp(D->fieldE,fileName,"/Ex",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolLowCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Ex",D->reNx+5,1,1,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ey",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolLowCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Ey",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ez",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolLowCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Ez",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);
   
   //save B field
   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   if(myrank==0)  
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[0],1); 
   else;
   MPI_Bcast(&minXDomain[0],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldOld,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   if(myrank==0)  
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[1],1); 
   else;
   MPI_Bcast(&minXDomain[1],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldNow,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   if(myrank==0)  
     restoreIntMeta(fileName,"/minXDomain",&minXDomain[2],1); 
   else;
   MPI_Bcast(&minXDomain[2],1,MPI_INT,0,MPI_COMM_WORLD);
   restoreFieldComp(D->fieldNext,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Bx",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/By",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Bz",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   //calculation current
   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Jx",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Jy",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calLowBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,1,0,1,D->resolX,minXDomain);
   resolLowCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Jz",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   deleteField(reField,D->reNxSub+5,1,1); 

   D->minXDomain/=D->resolX;
   if(myrank==0) {
     saveIntMeta(outFile,"minXDomain",&D->minXDomain,1);
     saveIntMeta(outFile,"minYDomain",&D->minYDomain,1);
     saveIntMeta(outFile,"nSpecies",&D->nSpecies,1);
     tmp=D->nx/D->resolX;
     saveIntMeta(outFile,"nx",&tmp,1);
     saveIntMeta(outFile,"ny",&D->ny,1);
     saveIntMeta(outFile,"nz",&D->nz,1);
     printf("%s is made\n",outFile);
   } else;
}

void resolLowCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int edge)
{
  int index,indexJ,i,j,k,n;
  double y1,y2,y3,a,b,c,x;

  if(resolY>1) {
    if(edge==1)  {
      for(i=istart; i<iend; i++)  
        for(j=jstart; j<jend; j+=resolY)  {
          indexJ=(j-jstart)/resolY+jstart;
          for(k=kstart; k<kend; k++)
            reField[i][indexJ][k]=field[i][j][k];
        }
    }  else {         //edge==0
      if(resolY%2==0)  {
        for(i=istart; i<iend; i++) 
          for(j=jstart; j<jend; j+=resolY)  {
            index=j+resolY*0.5;
            indexJ=(j-jstart)/resolY+jstart;
            for(k=kstart; k<kend; k++)   {
              y1=field[i][index-1][k];
              y2=field[i][index][k];
              y3=field[i][index+1][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              x=0.5;
              reField[i][indexJ][k]=a*x*x+b*x+c;
//              reField[i][indexJ][k]=0.5*(y1+y2);
            }
          }
      } else  {
        for(i=istart; i<iend; i++)  
          for(j=jstart; j<jend; j+=resolY) {
            index=j+(resolY-1)*0.5;
            indexJ=(j-jstart)/resolY+jstart;
            for(k=kstart; k<kend; k++)  {
              y1=field[i][index][k];
              reField[i][indexJ][k]=y1;
            }
          }
      }
    }             //End of edge==0
  } else {
    for(i=istart; i<iend; i++)  
      for(j=jstart; j<jend; j++) 
        for(k=kstart; k<kend; k++)
          reField[i][j][k]=field[i][j][k];
  }
}

void resolLowCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int edge)
{
  int index,indexI,i,j,k,n;
  double y1,y2,y3,a,b,c,x;

//lala
  if(resolX>1)  {
    if(edge==1)  {
      for(i=istart; i<iend; i+=resolX)    {
        indexI=(i-istart)/resolX+istart;
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
            reField[indexI][j][k]=field[i][j][k];
      }
    }
    else  {        //edge==0 
      if(resolX%2==0)  {
        for(i=istart; i<iend; i+=resolX)  {
          index=i+resolX*0.5;
          indexI=(i-istart)/resolX+istart;
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)   {
              y1=field[index-1][j][k];
              y2=field[index][j][k];
              y3=field[index+1][j][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              x=0.5;
              reField[indexI][j][k]=a*x*x+b*x+c;
            }
        }
      } else  {
        for(i=istart; i<iend; i+=resolX)  {
          index=i+(resolX-1)*0.5;
          indexI=(i-istart)/resolX+istart;
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)   {
              y1=field[index][j][k];
              reField[indexI][j][k]=y1;
            }
        }
      }
    }             //End of edge==0
  } else {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
          reField[i][j][k]=field[i][j][k];
  } 
}

void calLowBField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX, int *minXDomain)
{
  int i,j,k,shift1,shift2;
  double y1,y2,y3,a,b,c,x,w;

  shift1=minXDomain[1]-minXDomain[0];
  shift2=minXDomain[2]-minXDomain[1];
  if(resolX>1) {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)  {
          y1=fieldOld[i+shift1][j][k];
          y2=fieldNow[i][j][k];
          y3=fieldNext[i-shift2][j][k];
          a=0.5*(y1+y3)-y2;
          b=2.0*y2-1.5*y1-0.5*y3;
          c=y1;
          x=2.0-0.5*(resolX-1);
          fieldNow[i][j][k]=(a*x*x+b*x+c);
        }
  } else ;
}

void calLowJField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX, int *minXDomain)
{
  int i,j,k,shift1,shift2;
  double y1,y2,y3,a,b,c,x,w;

  shift1=minXDomain[1]-minXDomain[0];
  shift2=minXDomain[2]-minXDomain[1];
  if(resolX>1) {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)  {
          y1=fieldOld[i+shift1][j][k];
          y2=fieldNow[i][j][k];
          y3=fieldNext[i-shift2][j][k];
          a=0.5*(y1+y3)-y2;
          b=2.0*y2-1.5*y1-0.5*y3;
          c=y1;
          x=2.0-0.5*(resolX-1);
          fieldNow[i][j][k]=(a*x*x+b*x+c)*resolX;
        }
  } else ;
}
