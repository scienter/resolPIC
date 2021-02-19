#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_linalg.h>
#include "mesh.h"

void findB(double *a_data,double *c_data,double *result,int row,int col);
double gaussian(double sigma);

void calMomentum1D(Domain *D,int s)
{
  int i,j,k,l,m,cnt;
  int istart,iend,jstart,jend,kstart,kend;
  double x,invW,tmp,pxSam,pySam,pzSam,sigX,sigY,sigZ;
  double lambda[4],w[2],Px[2],Py[2],Pz[2],bx[2],by[2],bz[2];
  double bbx[2],bby[2],bbz[2],Cx[2],Cy[2],Cz[2];
  ptclList *p;

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  kstart=D->kstart; kend=D->kend;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        for(l=0; l<4; l++)  lambda[l]=0.0;
        Px[0]=D->px1[i][j][k];  Px[1]=D->px2[i][j][k];
        Py[0]=D->py1[i][j][k];  Py[1]=D->py2[i][j][k];
        Pz[0]=D->pz1[i][j][k];  Pz[1]=D->pz2[i][j][k];
        sigX=D->sigX[i][j][k];
        sigY=D->sigY[i][j][k];
        sigZ=D->sigZ[i][j][k];

        p=D->particle[i][j][k].head[s]->pt;
        cnt=0;
        while(p)  {
          x=p->x; w[0]=1.0-x;  w[1]=x;
          for(l=0; l<2; l++)
            for(m=0; m<2; m++)
              lambda[l*2+m]+=w[l]*w[m];

          if(cnt%2==0)  {
            pxSam=gaussian(sigX);
            pySam=gaussian(sigY);
            pzSam=gaussian(sigZ);
          } else   {
            pxSam*=-1.0;
            pySam*=-1.0;
            pzSam*=-1.0;
          }
 
          Cx[0]-=w[0]*pxSam; Cy[0]-=w[0]*pySam; Cz[0]-=w[0]*pzSam;
          Cx[1]-=w[1]*pxSam; Cy[1]-=w[1]*pySam; Cz[1]-=w[1]*pzSam;
          p->pxSam=pxSam; p->pySam=pySam; p->pzSam=pzSam;
          p=p->next; cnt++;
        }

        if(cnt>0)  {
          findB(lambda,Px,bbx,2,2);
          findB(lambda,Py,bby,2,2);
          findB(lambda,Pz,bbz,2,2);
          findB(lambda,Cx,bx,2,2);
          findB(lambda,Cy,by,2,2);
          findB(lambda,Cz,bz,2,2);
        }  else;
  
        p=D->particle[i][j][k].head[s]->pt;
        while(p)  {
          x=p->x; invW=1.0/p->weight;
          w[0]=1.0-x;  w[1]=x;
          tmp=0.0; for(l=0; l<2; l++) tmp+=bbx[l]*w[l];
          p->px=tmp*invW;          
          tmp=0.0; for(l=0; l<2; l++) tmp+=bby[l]*w[l];
          p->py=tmp*invW;          
          tmp=0.0; for(l=0; l<2; l++) tmp+=bbz[l]*w[l];
          p->pz=tmp*invW;          
          tmp=0.0; for(l=0; l<2; l++) tmp+=bx[l]*w[l];
          p->pxSam+=tmp;
          tmp=0.0; for(l=0; l<2; l++) tmp+=by[l]*w[l];
          p->pySam+=tmp;
          tmp=0.0; for(l=0; l<2; l++) tmp+=bz[l]*w[l];
          p->pzSam+=tmp;

          p=p->next;
        }

      }
}

void findB(double *a_data,double *c_data,double *result,int row,int col)
{
  int s,i,j;
  double data[row*col];

  for(i=0; i<row*col; i++)
    data[i]=a_data[i];

  gsl_matrix_view m  =gsl_matrix_view_array(data,row,col);
  gsl_vector_view c  = gsl_vector_view_array(c_data,row);
//  gsl_matrix_view inv=gsl_matrix_view_array(invA,row,col);

  gsl_vector *x = gsl_vector_alloc(col);
  gsl_permutation *p=gsl_permutation_alloc(row);
  gsl_linalg_LU_decomp(&m.matrix,p,&s);
  gsl_linalg_LU_solve(&m.matrix,p,&c.vector,x);

  for(i=0; i<row; i++)
    result[i]=gsl_vector_get(x,i);

  gsl_permutation_free(p);
  gsl_vector_free(x);
}

double gaussian(double sigma)
{
   double r,prob,v,random;
   int intRand,randRange=1000;

   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v=random*sigma;
      prob=exp(-v*v/sigma/sigma);
   }
   return v;
}


