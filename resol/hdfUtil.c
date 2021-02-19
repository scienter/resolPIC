#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
    int ii,i,j,k,start;
    double *field;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=ny;
    dimsf[1]=nx;
    dimsf[2]=nz;
    filespace=H5Screate_simple(3,dimsf,NULL);

    count[0]=nySub;
    count[1]=nxSub;
    count[2]=nzSub;
    offset[0]=offSet[1];
    offset[1]=offSet[0];
    offset[2]=offSet[2];
    memspace=H5Screate_simple(3,count,NULL);

    field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
    start=0;
      for(j=jstart; j<jend; j++)
        for(i=istart; i<iend; i++)
        {
          for(k=kstart; k<kend; k++)
            field[start+k-kstart]=data[i][j][k];
          start+=nzSub;
        }

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny;
  dimsf[1]=nx;
  dimsf[2]=nz;
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=nySub;
  count[1]=nxSub;
  count[2]=nzSub;
  offset[0]=offSet[1];
  offset[1]=offSet[0];
  offset[2]=offSet[2];
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nxSub*nySub*nzSub*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);

  start=0;
  for(j=jstart; j<jend; j++)
    for(i=istart; i<iend; i++)
    {
      for(k=kstart; k<kend; k++)
        data[i][j][k]=field[start+k-kstart];
      start+=nzSub;
    }

  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}



void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

double randomValue(double beta)
{
   double r;
   int intRand, randRange=10000, rangeDev;

   rangeDev=(int)(randRange*(1.0-beta));
   intRand = rand() % (randRange-rangeDev);
   r = ((double)intRand)/randRange+(1.0-beta);

   return r;
}

