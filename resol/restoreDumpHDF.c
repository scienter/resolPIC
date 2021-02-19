#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void restoreParticleHDF(Domain *D,int iteration,Particle ***particle)
{
    int i,j,k,s,n,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int indexI,indexJ,indexK;
    int nxSub,nySub,nzSub,nSpecies,totalCnt,start,core,dataCnt;
    int L,M,N,flag,shareCnt;
    int startIndexX,startIndexY,startIndexZ,unitX,unitY,unitZ;
    int rankX,rankY,rankZ,tmp,rank,cntSub,remain,sub;
    int *recv,*minXSub,*maxXSub,*minYSub,*maxYSub,*minZSub,*maxZSub;
    int *sharePNum,*recvDataCnt,*dataCore,*coreCnt,*dataIndex,*dataCores;
    double *data,minXDomain,minYDomain,minZDomain,shiftX,shiftY,shiftZ;
    double **sendData,**recvData,*shareData;
    double x,y,z,px,py,pz,resolX,resolY,resolZ;
    int offset[3];
    char name[100],dataName[100];
    void restoreIntMeta();
    void restoreFieldComp();
    ptclList *p;
    FILE *out;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;
    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
    L=D->L;  M=D->M;  N=D->N;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    //particle resotre
    recv=(int *)malloc(nTasks*sizeof(int ));
    minXSub=(int *)malloc(nTasks*sizeof(int ));
    maxXSub=(int *)malloc(nTasks*sizeof(int ));
    minYSub=(int *)malloc(nTasks*sizeof(int ));
    maxYSub=(int *)malloc(nTasks*sizeof(int ));
    minZSub=(int *)malloc(nTasks*sizeof(int ));
    maxZSub=(int *)malloc(nTasks*sizeof(int ));
    sharePNum=(int *)malloc(nTasks*sizeof(int ));
    recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    coreCnt=(int *)malloc(nTasks*sizeof(int ));

    if(D->mode==LOW) {
      resolX=1.0/((double)D->resolX);
      resolY=1.0/((double)D->resolY);
      resolZ=1.0/((double)D->resolZ);
    } else if(D->mode==HIGH) {
      resolX=(double)D->resolX;
      resolY=(double)D->resolY;
      resolZ=(double)D->resolZ;
    } else ;
    //restore minSub, maxSub
    unitX=((int)(D->nx*resolX))/D->L+1;
    unitY=((int)(D->ny*resolY))/D->M+1;
    unitZ=((int)(D->nz*resolZ))/D->N+1;
    MPI_Gather(&D->minXSub,1,MPI_INT,minXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxXSub,1,MPI_INT,maxXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->minYSub,1,MPI_INT,minYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxYSub,1,MPI_INT,maxYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->minZSub,1,MPI_INT,minZSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minZSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxZSub,1,MPI_INT,maxZSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxZSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);

    for(i=0; i<nTasks; i++) { 
      minXSub[i]=(int)(minXSub[i]*resolX);
      maxXSub[i]=(int)(maxXSub[i]*resolX);
      minYSub[i]=(int)(minYSub[i]*resolY);
      maxYSub[i]=(int)(maxYSub[i]*resolY);
      minZSub[i]=(int)(minZSub[i]*resolZ);
      maxZSub[i]=(int)(maxZSub[i]*resolZ);
    }
    minXDomain=D->minXDomain*resolX;
    minYDomain=D->minYDomain*resolY;
    minZDomain=D->minZDomain*resolZ;
    if(D->mode==LOW) {
      if(D->resolY>1) {
        if(D->reNy%2==0) shiftY=0;
        else	     shiftY=D->resolY;
      } else {
        shiftY=0;
      }
      shiftX=0.0;
    } else {
      shiftX=0.0;
      shiftY=0.0;
      shiftZ=0.0;
    }

    //restore particle
    sprintf(name,"dumpResolParticle%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(name,"/nSpecies",&nSpecies,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

    int cntList[nSpecies];

    for(s=0; s<nSpecies; s++)  {
      sprintf(dataName,"%dtotalCnt",s);
      if(myrank==0)  
        restoreIntMeta(name,dataName,&cntList[s],1);
      else	;
      MPI_Bcast(&cntList[s],1,MPI_INT,0,MPI_COMM_WORLD);
    }

    //set HDF5 parameter
    hid_t file_id,dset_id,plist_id,attr_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
    herr_t ierr;


    for(s=0; s<nSpecies; s++)
    {
      totalCnt=cntList[s];
      if(totalCnt>0)
      {
        //open file
        plist_id=H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
        file_id=H5Fopen(name,H5F_ACC_RDWR,plist_id);
        H5Pclose(plist_id);

        //set dataset
        sprintf(dataName,"%d",s);
        dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmp=sub+1;
          else	         tmp=sub;
          if(myrank==rank)
            cntSub=tmp;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvDataCnt[i]=0;
          coreCnt[i]=0;
        }

        data = (double *)malloc(cntSub*10*sizeof(double ));

        //file space
        dimsf[0]=totalCnt;
        dimsf[1]=10;
        filespace=H5Screate_simple(2,dimsf,NULL);

        //memory space
        dimsf[0]=cntSub;
        dimsf[1]=10;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cntSub;
        block[1]=10;
        offSet[0]=start;
        offSet[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);

        //hyperslab in memory space
        offSet[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);

        //read data
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Fclose(file_id);

        //for testing the core of each particle
        dataCore = (int *)malloc(cntSub*sizeof(int ));       
        for(i=0; i<cntSub; i++)
          dataCore[i]=-1;

        for(i=0; i<cntSub; i++)  {          
          x=(data[i*10+0]+shiftX)*resolX-minXDomain;         
          y=(data[i*10+1]+shiftY)*resolY;          
          z=data[i*10+2]*resolZ;          
          startIndexX=((int)x)/unitX;
          startIndexY=((int)(y-minYDomain))/unitY;
          startIndexZ=((int)(z-minZDomain))/unitZ;
          rank=startIndexY+startIndexZ*M+startIndexX*M*N;

          flag=0;
          while(flag==0 && rank<nTasks)  {
            if(minXSub[rank]<=x && x<maxXSub[rank] &&
               minYSub[rank]<=y && y<maxYSub[rank] && 
               minZSub[rank]<=z && z<maxZSub[rank]) { 
              flag=1;           

              if(rank==myrank)  {
                indexI=(int)(x-minXSub[rank])+D->istart;
                indexJ=(int)(y-minYSub[rank])+D->jstart;
                indexK=(int)(z-minZSub[rank])+D->kstart;

                p = (ptclList *)malloc(sizeof(ptclList));
                p->next = particle[indexI][indexJ][indexK].head[s]->pt;
                particle[indexI][indexJ][indexK].head[s]->pt=p;
            
                p->x=(x-minXSub[rank])-(int)(x-minXSub[rank]);
                p->y=(y-minYSub[rank])-(int)(y-minYSub[rank]);
                p->z=(z-minZSub[rank])-(int)(z-minZSub[rank]);
                p->px=data[i*10+3];
                p->py=data[i*10+4];
                p->pz=data[i*10+5];
                p->index=data[i*10+6];
                p->core=data[i*10+7];
                p->weight=data[i*10+8];
                p->charge=data[i*10+9];

              }
              dataCore[i]=rank;
              sharePNum[rank]+=1;
            }
            else if(x<0.0 || x>=D->reNx)  {
              dataCore[i]=-1;
              flag=1;
            } else  {
              rank++	;
            }
          }	//End of while(flag=0)

        }

        //set count '0' at myrank due to no reason for sharing
        for(i=0; i<nTasks; i++)
          if(myrank==i)  sharePNum[i]=0;

        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    
            MPI_Send(&sharePNum[i],1,MPI_INT,i,myrank,MPI_COMM_WORLD);         
        }
        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    {
            MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&mpi_status);
          }  else	;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //set memory for send and recving data
        dataCnt=10;
        for(i=0; i<nTasks; i++)   {          
          sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
          recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
        }        

        for(i=0; i<cntSub; i++)  
        {          
          core=dataCore[i];
          if(myrank!=core && core>=0 && core<nTasks)  {
            n=coreCnt[core];
            sendData[core][n*dataCnt+0]=(data[i*10+0]+shiftX)*resolX-minXDomain;
            sendData[core][n*dataCnt+1]=(data[i*10+1]+shiftY)*resolY;
            sendData[core][n*dataCnt+2]=data[i*10+2]*resolZ;
            sendData[core][n*dataCnt+3]=data[i*10+3];
            sendData[core][n*dataCnt+4]=data[i*10+4];
            sendData[core][n*dataCnt+5]=data[i*10+5];
            sendData[core][n*dataCnt+6]=data[i*10+6];
            sendData[core][n*dataCnt+7]=data[i*10+7];
            sendData[core][n*dataCnt+8]=data[i*10+8];
            sendData[core][n*dataCnt+9]=data[i*10+9];
            coreCnt[core]+=1;
          }	else	;
        }

        for(i=0; i<nTasks; i++)
        {
          if(myrank==i)  {
            for(j=0; j<nTasks; j++)
              if(i!=j)
                MPI_Send(sendData[j],sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
          } 
          else  {
            MPI_Recv(recvData[i],recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&mpi_status);
            for(j=0; j<recvDataCnt[i]; j++)  {
              x=recvData[i][j*dataCnt+0];
              y=recvData[i][j*dataCnt+1];
              z=recvData[i][j*dataCnt+2];
              indexI=(int)(x-minXSub[myrank])+D->istart;
              indexJ=(int)(y-minYSub[myrank])+D->jstart;
              indexK=(int)(z-minZSub[myrank])+D->kstart;
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ][indexK].head[s]->pt;
              particle[indexI][indexJ][indexK].head[s]->pt=p;
            
              p->x=(x-minXSub[myrank])-(int)(x-minXSub[myrank]);
              p->y=(y-minYSub[myrank])-(int)(y-minYSub[myrank]);
              p->z=0.0;
//              p->z=(z-D->minZSub)-(int)(z-D->minZSub);
              p->px=recvData[i][j*dataCnt+3];
              p->py=recvData[i][j*dataCnt+4];
              p->pz=recvData[i][j*dataCnt+5];
              p->index=recvData[i][j*dataCnt+6];
              p->core=recvData[i][j*dataCnt+7];
              p->weight=recvData[i][j*dataCnt+8];
              p->charge=recvData[i][j*dataCnt+9];
            }
          } 
          MPI_Barrier(MPI_COMM_WORLD);
        }

        free(data);
        free(dataCore);
        for(i=0; i<nTasks; i++)   {          
          free(recvData[i]);
          free(sendData[i]);
        }        

      }	else ;	//End of totalCnt>0
    } 		//End of nSpecies 

    free(recv);
    free(minXSub);
    free(maxXSub);
    free(minYSub);
    free(maxYSub);
    free(minZSub);
    free(maxZSub);
    free(sharePNum);
    free(recvDataCnt);
    free(coreCnt);
    free(recvData);
    free(sendData);
//lala

}


