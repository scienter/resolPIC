#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void MPI_TransferRho_Xplus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Xminus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Yminus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_TransferRho_Yplus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_TransferRho_Zminus(Domain *D,double ***f1,int nx,int ny,int share);
void MPI_TransferRho_Zplus(Domain *D,double ***f1,int nx,int ny,int share);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void solveDensity1D(Domain *D,int s,double coef);
void solveDensity2D(Domain *D,int s,double coef);
void solveDensity3D(Domain *D,int s,double coef);
void saveCoordHDF(Domain *D,char *fileName);
void density_xdmf(int dimension,char *fileName,int nx,int ny,int nz,int s);
void setZero(double ***den,int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);
void deleteParticle(Domain *D,int i,int j,int k,int s);


void saveDensityHDF(Domain *D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,offSetY,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    int offset[3];
    double *rho0;
    double charge;
    char name[100],dataName[100],fileName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub; nySub=D->nySub; nzSub=D->nzSub;
    istart=D->istart;    iend=nxSub+2;
    jstart=D->jstart;    jend=nySub+2;
    kstart=D->kstart;    kend=nzSub+2;

    rankX=myrank/(D->M*D->N);
    rankZ=(myrank%(D->M*D->N))/D->M;
    rankY=(myrank%(D->M*D->N))%D->M;

    s=0;
    rho0 = (double *)malloc((D->nSpecies)*sizeof(double ));
    LL=D->loadList;
    while(LL->next)    {
      if(LL->charge<0)	charge=-1.0*LL->charge;
      else 	        charge=LL->charge;
      rho0[s]=1.0*LL->density;
      LL=LL->next;
      s++;
    }

    nx=D->nx;  ny=D->ny;    nz=D->nz;

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=D->minZSub-D->minZDomain;

    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%ddensity%d.h5",s,iteration);

      if(myrank==0)        {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
        saveIntMeta(name,"/nx",&nx,1);
        saveIntMeta(name,"/ny",&ny,1);
        saveIntMeta(name,"/nz",&nz,1);
      } else      ;
      MPI_Barrier(MPI_COMM_WORLD);

      switch(D->dimension) {
      //1D
      case 1:
        solveDensity1D(D,s,rho0[s]);
        if(D->L>1)  {
          MPI_TransferRho_Xplus(D,D->Rho,1,1,3);
          MPI_TransferRho_Xminus(D,D->Rho,1,1,3);
        }  else     ;

        sprintf(dataName,"%d",s);
        saveFieldComp(D->Rho,name,dataName,nx,1,1,nxSub,1,1,istart,iend,0,1,0,1,offset);

        if(myrank==0)   {
          saveCoordHDF(D,name);
          sprintf(fileName,"%ddensity%d",s,iteration);
          density_xdmf(1,fileName,nx,ny,nz,s);
          printf("%s\n",name);
        }  else  ;
        break;

      //2D
      case 2:
        solveDensity2D(D,s,rho0[s]);
        if(D->L>1)  {
          MPI_TransferRho_Xplus(D,D->Rho,nySub+5,1,3);
          MPI_TransferRho_Xminus(D,D->Rho,nySub+5,1,3);
        }  else     ;
        if(D->M>1)  {
          MPI_TransferRho_Yplus(D,D->Rho,nxSub+5,1,3);
          MPI_TransferRho_Yminus(D,D->Rho,nxSub+5,1,3);
        }  else     ;

        sprintf(dataName,"%d",s);
        saveFieldComp(D->Rho,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);

        if(myrank==0)   {
          saveCoordHDF(D,name);
          sprintf(fileName,"%ddensity%d",s,iteration);
          density_xdmf(2,fileName,nx,ny,nz,s);
          printf("%s\n",name);
        }  else  ;

        break;

      //3D
      case 3:
        solveDensity3D(D,s,rho0[s]);
        if(D->L>1)  {
          MPI_TransferRho_Xplus(D,D->Rho,nySub+5,nzSub+5,3);
          MPI_TransferRho_Xminus(D,D->Rho,nySub+5,nzSub+5,3);
        }  else     ;
        if(D->M>1)  {
          MPI_TransferRho_Yplus(D,D->Rho,nxSub+5,nzSub+5,3);
          MPI_TransferRho_Yminus(D,D->Rho,nxSub+5,nzSub+5,3);
        }  else     ;
        if(D->N>1)  {
          MPI_TransferRho_Zplus(D,D->Rho,nxSub+5,nySub+5,3);
          MPI_TransferRho_Zminus(D,D->Rho,nxSub+5,nySub+5,3);
        }  else     ;

        sprintf(dataName,"%d",s);
        saveFieldComp(D->Rho,name,dataName,nx,ny,nz,nxSub,nySub,nzSub,istart,iend,jstart,jend,kstart,kend,offset);

        if(myrank==0)   {
          saveCoordHDF(D,name);
          sprintf(fileName,"%ddensity%d",s,iteration);
          density_xdmf(2,fileName,nx,ny,nz,s);
          printf("%s\n",name);
        }  else  ;

        break;
      }   //End of switch (dimension)

    }		//End of for(s)
    free(rho0);
}

void solveDensity1D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  double Wx[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,weight;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;

  j=k=0;
  for(i=0; i<iend+3; i++)
      D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        weight=p->weight;
        x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        for(ii=0; ii<4; ii++)
            D->Rho[i-1+ii][j][k]+=Wx[ii]*coef*weight;
        p=p->next;
      }
    }
}   

void solveDensity2D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  double Wx[4],Wy[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,weight;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;

  k=0;
  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        weight=p->weight;
        x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        y=p->y; y1=1+y; y2=y; y3=1-y; y4=2-y;
        Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
        Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
        Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
        Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
        for(ii=0; ii<4; ii++)
          for(jj=0; jj<4; jj++)
            D->Rho[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*coef*weight;
        p=p->next;
      }
    }
}   

void solveDensity3D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,kk,istart,iend,jstart,jend,kstart,kend;
  double Wx[4],Wy[4],Wz[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4,weight;;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;  iend=D->iend;
  jstart=D->jstart;  jend=D->jend;
  kstart=D->kstart;  kend=D->kend;

  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      for(k=0; k<kend+3; k++)
        D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {
          weight=p->weight;
          x=p->x; x1=1+x; x2=x; x3=1-x; x4=2-x;
          Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
          Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
          Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
          Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
          y=p->y; y1=1+y; y2=y; y3=1-y; y4=2-y;
          Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
          Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
          Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
          Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
          z=p->z; z1=1+z; z2=z; z3=1-z; z4=2-z;
          Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
          Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
          Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
          Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;
          for(ii=0; ii<4; ii++)
            for(jj=0; jj<4; jj++)
              for(kk=0; kk<4; kk++)
                D->Rho[i-1+ii][j-1+jj][k-1+kk]+=Wx[ii]*Wy[jj]*Wz[kk]*coef*weight;
          p=p->next;
        }
      }
}  

void saveCoordHDF(Domain *D,char *fileName)
{
  int ii,i,nx,ny,nz;
  char name[100];
  double *xtic,*ytic,*ztic;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,tic_id;
  herr_t status;
  hid_t filespace;
  hsize_t dimy[1],dimx[1],dimz[1];
  const char *coorName[] = {"/Y","/X","/Z"};

  nx=D->nx;  ny=D->ny;  nz=D->nz;

  if(myrank==0)
  {
    sprintf(name,"%s",fileName);
    file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);

    switch(D->dimension) {
    //1D
    case 1:
      dimx[0]=nx;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++)
        xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;

      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate2(file_id,coorName[1],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);

      free(xtic);
      break;

    //2D
    case 2:
      dimx[0]=nx;      dimy[0]=ny;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++) xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++) ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      for(ii=0; ii<2; ii++)
      {
        if(ii==0) filespace=H5Screate_simple(1,dimy,NULL);
        else if(ii==1) filespace=H5Screate_simple(1,dimx,NULL);
        dset_id=H5Dcreate2(file_id,coorName[ii],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ii==0 ? ytic : xtic);
        H5Dclose(dset_id);
        H5Sclose(filespace);
      }
      free(xtic);
      free(ytic);
      break;

    //3D
    case 3:
      dimx[0]=nx;      dimy[0]=ny;      dimz[0]=nz;
      xtic=(double *)malloc(nx*sizeof(double));
      for(i=0;i<nx;i++) xtic[i]=(i+D->minXDomain)*D->lambda*D->dx;
      ytic=(double *)malloc(ny*sizeof(double));
      for(i=0;i<ny;i++) ytic[i]=(i+D->minYDomain)*D->lambda*D->dy;
      ztic=(double *)malloc(nz*sizeof(double));
      for(i=0;i<nz;i++) ztic[i]=(i+D->minZDomain)*D->lambda*D->dz;
      filespace=H5Screate_simple(1,dimy,NULL);
      dset_id=H5Dcreate2(file_id,"/Y",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ytic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate2(file_id,"/X",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
      filespace=H5Screate_simple(1,dimz,NULL);
      dset_id=H5Dcreate2(file_id,"/Z",H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ztic);
      H5Dclose(dset_id);
      H5Sclose(filespace);

      free(xtic);
      free(ytic);
      free(ztic);
      break;
    }
    H5Fclose(file_id);
  }
  else ;
}

// The number of cells in the X, Y dimensions
void density_xdmf(int dimension,char *fileName,int nx,int ny,int nz,int s)
{
    FILE *xmf = 0;
    char name[100];
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    switch (dimension)  {
    //1D
    case 1 :
      fprintf(xmf, "     <Topology TopologyType=\"1DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VX\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");

      fprintf(xmf, "     <Attribute Name=\"/%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/%d\n",fileName,s);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;

    //2D
    case 2 :
      fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/X\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s.h5:/Y\n",fileName);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");

      fprintf(xmf, "     <Attribute Name=\"/0%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/%d\n",fileName,s);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
/*
    case 3 :
      fprintf(xmf, "     <Topology TopologyType=\"3DRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",ny,nx,nz);
      fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz);
      fprintf(xmf, "        %s%d.h5:/Z\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s%d.h5:/X\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
      fprintf(xmf, "        %s%d.h5:/Y\n",name,iteration);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Geometry>\n");
      for(s=0; s<nSpecies; s++)
      {
        fprintf(xmf, "     <Attribute Name=\"%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
        fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nz,nx,ny);
        fprintf(xmf, "        %s%d.h5:/%d\n",name,iteration,s);
        fprintf(xmf, "       </DataItem>\n");
        fprintf(xmf, "     </Attribute>\n");
      }
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
      break;
*/
    }
}
                                                          
void setZero(double ***den,int nx, int ny, int nz)
{
   int i,j,k;

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         den[i][j][k]=0.0;
}

void deleteParticle(Domain *D,int i,int j,int k,int s)
{
  ptclList *p,*tmp;

  p=D->particle[i][j][k].head[s]->pt;
  while(p)  {
    tmp=p->next;
    D->particle[i][j][k].head[s]->pt=tmp;
    p->next=NULL;
    free(p);
    p=D->particle[i][j][k].head[s]->pt;
  }
}

