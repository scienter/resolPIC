#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void resolHighCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int edge);
void resolHighCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int edge);
void calHighBField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int *minXDomain);
void reHighCreateField1D(Domain *D);
void reHighCreateField2D(Domain *D);
double ***memoryAsign(int nx, int ny, int nz);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void deleteField(double ***field,int nx,int ny,int nz);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void calOneData(Domain *D,char *outFile,char *fileName,char *dataName,int step,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int Iend,int Jend,int Kend,int *offset,double ***reField,int edge,int reNx,int reNy,int reNz,int *offSet);



void reHighCreateField(Domain D)
{
   int istart,iend,jstart,jend,kstart,kend;

   istart=D.istart;    iend=D.iend;
   jstart=D.jstart;    jend=D.jend;
   kstart=D.kstart;    kend=D.kend;

   switch (D.dimension)  {
   case 1 :
//     reHighCreateField1D(&D);
     break;
   case 2 :
     reHighCreateField2D(&D);
     break;
   }

}


void reHighCreateField2D(Domain *D)
{
   int nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend;
   int offset[3],offSet[3],minXDomain[3],Iend,Jend,Kend;
   int nxSub1D,nySub2D,nzSub3D;
   double resolX;
   double ***reField1,***reField2;
   char fileName[100],outFile[100];
   hid_t file_id;
   
   int myrank, nTasks,shift=0;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   nx=D->nx; ny=D->ny; nz=D->nz;
   nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   kstart=D->kstart; kend=D->kend;

   Iend=iend+3; Jend=1; Kend=1;
   if(D->dimension>1) { Jend=jend+3; } else;
   if(D->dimension>2) { Kend=kend+3; } else;

   offset[0]=D->minXSub;
   offset[1]=D->minYSub-D->minYDomain;
   offset[2]=0;

   //for saving field file
   offSet[0]=(D->minXSub)*D->resolX+D->biasX;
   offSet[1]=(D->minYSub-D->minYDomain)*D->resolY+D->biasY;
   offSet[2]=0;

   resolX=(double)D->resolX;


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
   resolHighCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Ex",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ey",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   resolHighCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Ey",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ez",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   resolHighCalX(D->fieldE,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
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
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Bx",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/By",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/By",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Bz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNow,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Bz",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

    //calculation current
   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Jx",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jy",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,0);
   saveFieldComp(reField2,outFile,"/Jy",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jz",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,1);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
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

   D->nx=D->nx*D->resolX;
   D->ny=D->ny*D->resolY;
   D->nz=D->nz*D->resolZ;
   D->minXDomain*=D->resolX;
   D->minYDomain*=D->resolY;
   D->minZDomain*=D->resolZ;

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

/*
void reHighCreateField1D(Domain *D)
{
   int nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend;
   int offset[3],offSet[3];
   double ***reField;
   char fileName[100],outFile[100];
   hid_t file_id;
   
ub1D=D->nxSub+5; nySub2D=1; nzSub3D=1;
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

   restoreFieldComp(D->fieldOld,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jx",nx+5,ny+5,nz,nxSub+5,nySub+5,1,0,iend+3,0,jend+3,0,1,offset);
   calHighBField(D->fieldOld,D->fieldNow,D->fieldNext,1,iend+2,0,jend+3,0,1,D->resolX,minXDomain);
   resolHighCalX(D->fieldNext,reField1,2,iend,2,jend,0,1,D->resolX,0);
   resolHighCalY(reField1,reField2,2,D->reIend,2,jend,0,1,D->resolY,1);
   saveFieldComp(reField2,outFile,"/Jx",D->reNx+5,D->reNy+5,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);
t=0;
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
   offSet[0]=(D->minXSub)*D->resolX+D->biasX;
//   offSet[1]=(D->minYSub-D->minYDomain)*D->resolY+D->biasY;
   offSet[1]=0;
   offSet[2]=0;

   reField=memoryAsign(D->reNxSub+5,1,1);

   //create out file
   sprintf(outFile,"dumpField%d.h5",D->step*D->resolX);
   if(myrank==0)
   {
     file_id=H5Fcreate(outFile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
     H5Fclose(file_id);
   }
   else ;



   //save E field
   sprintf(fileName,"resolE%d.h5",D->step);
   restoreFieldComp(D->fieldE,fileName,"/Ex",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Ex",D->reNx+5,1,1,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ey",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Ey",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   restoreFieldComp(D->fieldE,fileName,"/Ez",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   resolCalX(D->fieldE,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Ez",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);
   
   //save B field
   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Bx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Bx",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/By",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/By",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolB%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolB%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Bz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Bz",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   //calculation current
   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jx",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,0);
   saveFieldComp(reField,outFile,"/Jx",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jy",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Jy",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   sprintf(fileName,"resolJ%d.h5",D->step-1+shift);
   restoreFieldComp(D->fieldOld,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+shift);
   restoreFieldComp(D->fieldNow,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   sprintf(fileName,"resolJ%d.h5",D->step+1+shift);
   restoreFieldComp(D->fieldNext,fileName,"/Jz",nx+5,ny,nz,nxSub+5,1,1,0,iend+3,0,1,0,1,offset);
   calBField(D->fieldOld,D->fieldNow,D->fieldNext,0,iend+3,0,1,0,1,D->resolX);
   resolCalX(D->fieldNow,reField,2,iend,0,1,0,1,D->resolX,1);
   saveFieldComp(reField,outFile,"/Jz",D->reNx+5,D->reNy,D->reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);

   deleteField(reField,D->reNxSub+5,1,1); 

   if(myrank==0) {
     D->minXDomain*=D->resolX;
     D->minYDomain*=D->resolY;
     D->nx*=D->resolX; D->ny*=D->resolY; D->nz*=D->resolZ;
     saveIntMeta(outFile,"minXDomain",&D->minXDomain,1);
     saveIntMeta(outFile,"minYDomain",&D->minYDomain,1);
     saveIntMeta(outFile,"nSpecies",&D->nSpecies,1);
     saveIntMeta(outFile,"nx",&D->nx,1);
     saveIntMeta(outFile,"ny",&D->ny,1);
     saveIntMeta(outFile,"nz",&D->nz,1);
     printf("%s is made\n",outFile);
   } else;
}

void calOneData(Domain *D,char *outFile,char *fileName,char *dataName,int step,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int Iend,int Jend,int Kend,int *offset,double ***reField,int edge,int reNx,int reNy,int reNz,int *offSet)
{
   char name[100],outName[100];

   sprintf(name,"%s%d.h5",fileName,step);
   sprintf(outName,"%s",outFile);
   restoreFieldComp(D->fieldE,name,dataName,nx,ny,nz,nxSub,nySub,nzSub,0,Jend,0,Jend,0,Kend,offset);
   resolCalX(D->fieldE,reField,D->istart,D->iend,D->jstart,D->jend,0,1,D->resolX,edge);
   saveFieldComp(reField,outName,dataName,reNx,reNy,reNz,D->saveNxSub,D->saveNySub,D->saveNzSub,D->saveIstart,D->saveIend,D->saveJstart,D->saveJend,D->saveKstart,D->saveKend,offSet);
}
*/

void resolHighCalX(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int edge)
{
  int ii,i,j,k,n,indexI;
  double y1,y2,y3,a,b,c,dx,x;

  dx=1.0/((double)resolX);

  if(resolX>1)  {
    if(edge==1)    {
      for(i=istart; i<iend; i++)  
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++) 
            for(n=0; n<resolX; n++)  {
              x=dx*n;
              ii=i+(int)(x+0.5);
              y1=field[ii-1][j][k];
              y2=field[ii][j][k];
              y3=field[ii+1][j][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              x=1.0-(int)(0.5+x)+x;
              indexI=(i-istart)*resolX+n+istart;
              reField[indexI][j][k]=a*x*x+b*x+c;
            }
    }  else   {        //edge==0
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
            for(n=0; n<resolX; n++)     {
              y1=field[i-1][j][k];
              y2=field[i][j][k];
              y3=field[i+1][j][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              x=0.5+0.5*dx+dx*n;
              indexI=(i-istart)*resolX+n+istart;
              reField[indexI][j][k]=a*x*x+b*x+c;
            }
    }             //End of edge==0
  }  else {		//End of if(resolX>1)
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
          reField[i][j][k]=field[i][j][k];
  }
}

void resolHighCalY(double ***field,double ***reField,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolY,int edge)
{
  int jj,i,j,k,n,indexJ;
  double y1,y2,y3,a,b,c,dy,y;

  dy=1.0/((double)resolY);

  if(resolY>1)  {
    if(edge==1)    {
      for(i=istart; i<iend; i++)  
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++) 
            for(n=0; n<resolY; n++)  {
              y=dy*n;
              jj=j+(int)(y+0.5);
              y1=field[i][jj-1][k];
              y2=field[i][jj][k];
              y3=field[i][jj+1][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              y=1.0-(int)(0.5+y)+y;
              indexJ=(j-jstart)*resolY+n+jstart;
              reField[i][indexJ][k]=a*y*y+b*y+c;
            }
    }  else   {        //edge==0
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
            for(n=0; n<resolY; n++)     {
              y1=field[i][j-1][k];
              y2=field[i][j][k];
              y3=field[i][j+1][k];
              a=0.5*(y1+y3)-y2;
              b=2.0*y2-1.5*y1-0.5*y3;
              c=y1;
              y=0.5+0.5*dy+dy*n;
              indexJ=(j-jstart)*resolY+n+jstart;
              reField[i][indexJ][k]=a*y*y+b*y+c;
            }
    }             //End of edge==0
  }  else {		//End of if(resolX>1)
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
          reField[i][j][k]=field[i][j][k];
  }
}

void calHighBField(double ***fieldOld,double ***fieldNow,double ***fieldNext,int istart,int iend,int jstart,int jend,int kstart,int kend,int resolX,int *minXDomain)
{
  int i,j,k,shift1,shift2;
  double y1,y2,y3,a,b,c,dx,x;

  shift1=minXDomain[1]-minXDomain[0];
  shift2=minXDomain[2]-minXDomain[1];
  if(resolX>1) {
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)      {
          y1=fieldOld[i+shift1][j][k];
          y2=fieldNow[i][j][k];
          y3=fieldNext[i-shift2][j][k];
          a=0.5*(y1+y3)-y2;
          b=2.0*y2-1.5*y1-0.5*y3;
          c=y1;
          dx=1.0/((double)resolX);
          x=0.5+dx*0.5;
          fieldNow[i][j][k]=(a*x*x+b*x+c);
        }
  } else ;
}
