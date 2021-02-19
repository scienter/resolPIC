#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

//void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
//void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
//void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
//void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);

void saveDumpParticleHDF(Domain *D)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,tmp;
    int cnt,totalCnt,index,start,cntList[D->nSpecies],dataCnt=10;
    int minXSub,minYSub,minZSub,nxSub,nySub,nzSub;
    int minXDomain,minYDomain,minZDomain;
    char name[100],dataName[100];
    double *data,resolX,resolY,resolZ,shift;
    int *recv;
    Particle ***particle;
    particle=D->particle;
    ptclList *p;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2],a_dims,metaDim[1];
    herr_t ierr;

    nxSub=D->nxSub;        nySub=D->nySub;        nzSub=D->nzSub;
    minXSub=D->minXSub;    minYSub=D->minYSub;    minZSub=D->minZSub;

    istart=D->istart;    iend=D->reIend;
    jstart=D->jstart;    jend=D->reJend;
    kstart=D->kstart;    kend=D->reKend;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    //create file
    if(D->mode==HIGH) {
      resolX=D->resolX;
      resolY=D->resolY;
      resolZ=D->resolZ;
      shift=0.0;
    } else if(D->mode==LOW) {
      resolX=1.0/((double)D->resolX);
      resolY=1.0/((double)D->resolY);
      resolZ=1.0/((double)D->resolZ);
      shift=0;
    }  else ;
    minXDomain=resolX*D->minXDomain;
    minYDomain=resolY*D->minYDomain;
    minZDomain=resolZ*D->minZDomain;
    sprintf(name,"dumpParticle%g.h5",D->step*resolX);
    
    file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    recv = (int *)malloc(nTasks*sizeof(int ));

    for(s=0; s<D->nSpecies; s++)
    {
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)          {
            p=particle[i][j][k].head[s]->pt;
            while(p)   {
              cnt++;
              p=p->next;
            }
          }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

      start=0;
      for(i=0; i<myrank; i++)        start+=recv[i];
      totalCnt=0;
      for(i=0; i<nTasks; i++)        totalCnt+=recv[i];
      cntList[s]=totalCnt;

      //file space
      dimsf[0]=totalCnt;
      dimsf[1]=10;
      filespace=H5Screate_simple(2,dimsf,NULL);
      sprintf(dataName,"%d",s);
      dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      if(totalCnt>0)
      {
        data = (double *)malloc(cnt*dataCnt*sizeof(double ));

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++)
            for(k=kstart; k<kend; k++)  {
              p=particle[i][j][k].head[s]->pt;
              while(p)    {
                data[index*dataCnt+0]=p->x+minXDomain;
                data[index*dataCnt+1]=p->y;
                data[index*dataCnt+2]=p->z;
                data[index*dataCnt+3]=p->px;
                data[index*dataCnt+4]=p->py;
                data[index*dataCnt+5]=p->pz;
                data[index*dataCnt+6]=p->index;
                data[index*dataCnt+7]=p->core;
                data[index*dataCnt+8]=p->weight;
                data[index*dataCnt+9]=p->charge;
                index++;
                p=p->next;
              }
            }

        //memory space
        dimsf[0]=cnt;
        dimsf[1]=dataCnt;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cnt;
        block[1]=dataCnt;
        offset[0]=start;
        offset[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);

        //hyperslab in memory space
        offset[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);
      
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
        H5Pclose(plist_id);
        H5Sclose(memspace);
      
        free(data);
      }	else	; 	//End of totalCnt>0

      MPI_Barrier(MPI_COMM_WORLD);

      H5Dclose(dset_id);

      H5Sclose(filespace);
    }	//End of nSpecies
    free(recv);

    H5Fclose(file_id);

    if(myrank==0)  {
      saveIntMeta(name,"nSpecies",&D->nSpecies,1);
      for(s=0; s<D->nSpecies; s++)  {
        sprintf(dataName,"%dtotalCnt",s);
        saveIntMeta(name,dataName,&cntList[s],1);
      }
      saveIntMeta(name,"minXDomain",&minXDomain,1);
      saveIntMeta(name,"minYDomain",&minYDomain,1);
      saveIntMeta(name,"minZDomain",&minZDomain,1);
      saveIntMeta(name,"nx",&D->nx,1);
      saveIntMeta(name,"ny",&D->ny,1);
      saveIntMeta(name,"nz",&D->nz,1);
      printf("%s is made\n",name);
    }  else	;

    MPI_Barrier(MPI_COMM_WORLD);
}

