#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "math.h"
#include "hdf5.h"
#include "parameter.h"
#include "constants.h"


void solveSpectrum(Parameter *D,int idList,int coreList,int species,double freq)
{
  double t,x,y,z,px,py,pz,gamma,phase,prevVx,prevVy,prevVz,step;
  double sinTh,cosTh,dt,n_dot_beta,n_dot_Dbeta,denomitor;
  double vx,vy,vz,DbetaX,DbetaY,DbetaZ,betaX,betaY,betaZ;
  double RR,RRx,RRy,RRz,Rx,Ry,Rz,nx,ny,nz,tmpX;
  double A[3];
  double *dataX,*dataY,*dataZ,*dataPx,*dataPy,*dataPz,*dataStep;
  int iteration,i,j,id,core,totalCnt,cnt;
  char fileName[100],dataName[100];
  void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt);
  void restoreIntMeta(char *fileName,char *dataName,int *data);
  
  dt=D->dt;
  sprintf(fileName,"%dTrack%d_%d.h5",species,idList,coreList);
  restoreIntMeta(fileName,"totalCnt",&totalCnt);

  dataX=(double *)malloc(totalCnt*sizeof(double ));
  dataY=(double *)malloc(totalCnt*sizeof(double ));
  dataZ=(double *)malloc(totalCnt*sizeof(double ));
  dataPx=(double *)malloc(totalCnt*sizeof(double ));
  dataPy=(double *)malloc(totalCnt*sizeof(double ));
  dataPz=(double *)malloc(totalCnt*sizeof(double ));
  dataStep=(double *)malloc(totalCnt*sizeof(double ));

  restoreFloatArray(fileName,"x",dataX,totalCnt);
  restoreFloatArray(fileName,"y",dataY,totalCnt);
  restoreFloatArray(fileName,"z",dataZ,totalCnt);
  restoreFloatArray(fileName,"px",dataPx,totalCnt);
  restoreFloatArray(fileName,"py",dataPy,totalCnt);
  restoreFloatArray(fileName,"pz",dataPz,totalCnt);
  restoreFloatArray(fileName,"step",dataStep,totalCnt);

  tmpX=0.0;
  cnt=0;
  while(tmpX==0.0)
  {
    tmpX=dataX[cnt];
    cnt++;
  }

  Rx=D->det_R;
  Ry=D->det_x;
  Rz=D->det_y;

  px=dataPx[cnt-1];
  py=dataPy[cnt-1];
  pz=dataPz[cnt-1];
  gamma=sqrt(1.0+px*px+py*py+pz*pz);
  prevVx=px/gamma;
  prevVy=py/gamma;
  prevVz=pz/gamma; //uz/gamma;
  
  for(i=cnt; i<totalCnt; i++)
  {
    t=dataStep[i]*dt;
    px=dataPx[i];
    py=dataPy[i];
    pz=dataPz[i];
    x=dataX[i];
    y=dataY[i];
    z=dataZ[i];
    gamma=sqrt(1.0+px*px+py*py+pz*pz);
    vx=px/gamma;
    vy=py/gamma;
    vz=pz/gamma; //uz/gamma;

    DbetaX=(vx-prevVx)/dt;
    DbetaY=(vy-prevVy)/dt;
    DbetaZ=(vz-prevVz)/dt;
    betaX=(vx+prevVx)*0.5;
    betaY=(vy+prevVy)*0.5;
    betaZ=(vz+prevVz)*0.5;
    
    RRx=Rx-x;
    RRy=Ry-y;
    RRz=Rz-z;
    RR=sqrt(RRx*RRx+RRy*RRy+RRz*RRz);
//    RR=sqrt(x*x+y*y+z*z);
    
    nx=RRx/RR;
    ny=RRy/RR;
    nz=RRz/RR;
    n_dot_beta=nx*betaX+ny*betaY+nz*betaZ;
    n_dot_Dbeta=nx*DbetaX+ny*DbetaY+nz*DbetaZ;
    denomitor=(1.0-n_dot_beta)*(1.0-n_dot_beta);
    A[0]=(n_dot_Dbeta*(nx-betaX)-DbetaX*(1.0-n_dot_beta))/denomitor*dt;
    A[1]=(n_dot_Dbeta*(ny-betaY)-DbetaY*(1.0-n_dot_beta))/denomitor*dt;
    A[2]=(n_dot_Dbeta*(nz-betaZ)-DbetaZ*(1.0-n_dot_beta))/denomitor*dt;

    phase=freq*(t+RR/velocityC);

    for(j=0; j<3; j++)
    {
      D->B[j]+=A[j]*cos(phase);
      D->C[j]+=A[j]*sin(phase);
    }

    prevVx=vx;
    prevVy=vy;
    prevVz=vz;
  }
  free(dataX);
  free(dataY);
  free(dataZ);
  free(dataPx);
  free(dataPy);
  free(dataPz);
  free(dataStep);

}

void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

