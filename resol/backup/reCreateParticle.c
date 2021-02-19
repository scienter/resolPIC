#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void deleteParticle(Particle ***particle,int i,int j,int k,int s);
void storeParticle(Domain *D,Particle ***particle,int i,int j,int k,int s,double *x2,double *y2,double *w2,double *px2,double *py2,double *pz2,double dx,double dy,int l,int m,double resolX,double resolY,double shift);
void reCreateParticle1D(Domain *D,int i,int j,int k,int s);
void reCreateParticle2D(Domain *D,int i,int j,int k,int s);
void calInterpolation(double ***data,int cnt,double resolX);

void reCreateParticle(Domain D)
{
   int i,j,k,s,istart,iend,jstart,jend,kstart,kend;

   istart=D.istart;    iend=D.iend;
   jstart=D.jstart;    jend=D.jend;
   kstart=D.kstart;    kend=D.kend;

   switch (D.dimension)  {
   case 1 :
     j=k=0;
     for(i=istart; i<iend; i++)
       for(s=0; s<D.nSpecies; s++)  
         reCreateParticle1D(&D,i,j,k,s);
     break;
   case 2 :
     k=0;
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++)
         for(s=0; s<D.nSpecies; s++)  
           reCreateParticle2D(&D,i,j,k,s);
     break;
   }

}

void reCreateParticle2D(Domain *D,int i,int j,int k,int s)
{
   int l,m,n,cnt,totalCnt,index,flag,indexI,indexJ;
   double gapX,W,x,y,px,py,pz,weight,tmp2,Pr,Pl;
   double minX,maxX,minY,maxY,dx,dy,tmpX,tmpY,x0,y0;
   double resolX,resolY,shift;
   double ***Q,***Px,***Py,***Pz;
   double px2[2],py2[2],pz2[2],x2[2],y2[2],w2[2],w[4];
   size_t N,stride;
   Particle ***particle,***befoParticle,***nextParticle;
   particle=D->particle;
   ptclList *p;



   tmp2=sqrt(D->targetW[s]);
   cnt=(int)(1.0/tmp2);
   dx=1.0/((double)cnt);
   dy=dx;

   Q=(double ***)malloc(cnt*sizeof(double **));
   Px=(double ***)malloc(cnt*sizeof(double **));
   Py=(double ***)malloc(cnt*sizeof(double **));
   Pz=(double ***)malloc(cnt*sizeof(double **));
   for(l=0; l<cnt; l++) {
     Q[l]=(double **)malloc(cnt*sizeof(double *));
     Px[l]=(double **)malloc(cnt*sizeof(double *));
     Py[l]=(double **)malloc(cnt*sizeof(double *));
     Pz[l]=(double **)malloc(cnt*sizeof(double *));
     for(m=0; m<cnt; m++) {
       Q[l][m]=(double *)malloc(4*sizeof(double *));
       Px[l][m]=(double *)malloc(4*sizeof(double *));
       Py[l][m]=(double *)malloc(4*sizeof(double *));
       Pz[l][m]=(double *)malloc(4*sizeof(double *));
     }
   }
   for(l=0; l<cnt; l++)
     for(m=0; m<cnt; m++)
       for(n=0; n<4; n++) {
         Q[l][m][n]=0.0; Px[l][m][n]=0.0;
         Py[l][m][n]=0.0; Pz[l][m][n]=0.0;
       }

   p=particle[i][j][k].head[s]->pt;
   while(p)  {
     x=p->x; y=p->y; weight=p->weight;
     indexI=(int)(x/dx);  indexJ=(int)(y/dy);
     minX=indexI*dx;  maxX=(indexI+1)*dx;
     minY=indexJ*dy;  maxY=(indexJ+1)*dy;
     px=p->px;     py=p->py;         pz=p->pz;
     tmp2=weight/dx/dy;
     w[0]=(maxX-x)*(maxY-y);
     w[1]=(x-minX)*(maxY-y);
     w[2]=(maxX-x)*(y-minY);
     w[3]=(x-minX)*(y-minY);
     for(n=0; n<4; n++)  {
       Q[indexI][indexJ][n]+=w[n]*tmp2;    
       Px[indexI][indexJ][n]+=w[n]*tmp2*px;    
       Py[indexI][indexJ][n]+=w[n]*tmp2*py;    
       Pz[indexI][indexJ][n]+=w[n]*tmp2*pz;    
     }
     p=p->next;
   }
   deleteParticle(particle,i,j,k,s);

   if(D->mode==HIGH) {  
     resolX=D->resolX; resolY=D->resolY;
     shift=-1.0/resolX;
   } else if(D->mode==LOW) {
     resolX=1.0/((double)D->resolX);
     resolY=1.0/((double)D->resolY);
     shift=0;
   } else ;
 
   for(l=0; l<cnt; l++)  
     for(m=0; m<cnt; m++)  {
       W=0.0;
       for(n=0; n<4; n++)  W+=Q[l][m][n];

       if(W>0)  {
         x0=(Q[l][m][1]+Q[l][m][3])/W;
         y0=(Q[l][m][2]+Q[l][m][3])/W;
         tmpX=(x0-0.5)*(x0-0.5);
         tmpY=(y0-0.5)*(y0-0.5);

         if(tmpY<tmpX) {
           w2[0]=W*(1.0-y0);      w2[1]=W*y0;
           x2[0]=x0;          x2[1]=x0;
           y2[0]=y0*0.5;     y2[1]=(1.0+y0)*0.5;
           Pr=Px[l][m][2]+Px[l][m][3]; Pl=Px[l][m][0]+Px[l][m][1];
           px2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Py[l][m][2]+Py[l][m][3]; Pl=Py[l][m][0]+Py[l][m][1];
           py2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Pz[l][m][2]+Pz[l][m][3]; Pl=Pz[l][m][0]+Pz[l][m][1];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
         } else {
           w2[0]=W*(1.0-x0);      w2[1]=W*x0;
           x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
           y2[0]=y0;              y2[1]=y0;
           Pr=Px[l][m][1]+Px[l][m][3]; Pl=Px[l][m][0]+Px[l][m][2];
           px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Py[l][m][1]+Py[l][m][3]; Pl=Py[l][m][0]+Py[l][m][2];
           py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Pz[l][m][1]+Pz[l][m][3]; Pl=Pz[l][m][0]+Pz[l][m][2];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         }

         storeParticle(D,particle,i,j,k,s,x2,y2,w2,px2,py2,pz2,dx,dy,l,m,resolX,resolY,shift);
       }	//End of W>0

     }	//Enf for l,m

   for(l=0; l<cnt; l++) {
     for(m=0; m<cnt;m++) {
       free(Q[l][m]);       free(Px[l][m]);
       free(Py[l][m]);      free(Pz[l][m]);
     }
     free(Q[l]);       free(Px[l]);
     free(Py[l]);      free(Pz[l]);
   }
   free(Q);   free(Px);
   free(Py);  free(Pz);
}

void reCreateParticle1D(Domain *D,int i,int j,int k,int s)
{
   int l,m,n,cnt,totalCnt,index,flag,indexI;
   double minX,maxX,gapX,x0,W,x,px,py,pz,tmp2,dx,weight,targetW;
   double ***Q,***Px,***Py,***Pz,shift;
   double px2[2],py2[2],pz2[2],x2[2],y2[2],w2[2],resolX;
   size_t N,stride;
   Particle ***particle,***befoParticle,***nextParticle;
   particle=D->particle;
//   befoParticle=D->befoParticle;
//   nextParticle=D->nextParticle;
   ptclList *p;

   targetW=D->targetW[s];
   cnt=(int)(1.0/targetW);
   dx=1.0/((double)cnt);

   Q=(double ***)malloc(cnt*sizeof(double **));
   Px=(double ***)malloc(cnt*sizeof(double **));
   Py=(double ***)malloc(cnt*sizeof(double **));
   Pz=(double ***)malloc(cnt*sizeof(double **));
   for(l=0; l<cnt; l++) {
     Q[l]=(double **)malloc(2*sizeof(double *));
     Px[l]=(double **)malloc(2*sizeof(double *));
     Py[l]=(double **)malloc(2*sizeof(double *));
     Pz[l]=(double **)malloc(2*sizeof(double *));
     for(m=0; m<2; m++) {
       Q[l][m]=(double *)malloc(3*sizeof(double *));
       Px[l][m]=(double *)malloc(3*sizeof(double *));
       Py[l][m]=(double *)malloc(3*sizeof(double *));
       Pz[l][m]=(double *)malloc(3*sizeof(double *));
     }
   }
   for(l=0; l<cnt; l++)
     for(m=0; m<2; m++)
       for(n=0; n<3; n++) {
         Q[l][m][n]=0.0; Px[l][m][n]=0.0;
         Py[l][m][n]=0.0; Pz[l][m][n]=0.0;
       }
   y2[0]=0.0; y2[1]=0.0;

   p=particle[i][j][k].head[s]->pt;
   while(p)  {
     x=p->x;       weight=p->weight;
     indexI=(int)(x/dx);
     minX=indexI*dx;  maxX=(indexI+1)*dx;
     px=p->px;     py=p->py;         pz=p->pz;
     tmp2=weight/dx;
     Q[indexI][0][1]+=(maxX-x)*tmp2;    
     Q[indexI][1][1]+=(x-minX)*tmp2;
     Px[indexI][0][1]+=(maxX-x)*tmp2*px; 
     Px[indexI][1][1]+=(x-minX)*tmp2*px;
     Py[indexI][0][1]+=(maxX-x)*tmp2*py;
     Py[indexI][1][1]+=(x-minX)*tmp2*py;
     Pz[indexI][0][1]+=(maxX-x)*tmp2*pz;
     Pz[indexI][1][1]+=(x-minX)*tmp2*pz;
     p=p->next;
   }
   deleteParticle(particle,i,j,k,s);

   if(D->mode==HIGH) {
     resolX=D->resolX;
     shift=1.0/resolX;
   }
   else if(D->mode==LOW) {
     resolX=1.0/((double)D->resolX);
     shift=0.0; 
   }
   else ;
 
   for(l=0; l<cnt; l++)  {
     W=Q[l][0][1]+Q[l][1][1];
     if(W>0)  {
       x0=Q[l][1][1]/W;
       w2[0]=Q[l][0][1];      w2[1]=Q[l][1][1];
       x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
       px2[0]=((Px[l][0][1]-Px[l][1][1])+(Px[l][0][1]+Px[l][1][1])*x0)/w2[0];
       py2[0]=((Py[l][0][1]-Py[l][1][1])+(Py[l][0][1]+Py[l][1][1])*x0)/w2[0];
       pz2[0]=((Pz[l][0][1]-Pz[l][1][1])+(Pz[l][0][1]+Pz[l][1][1])*x0)/w2[0];
       px2[1]=(2*Px[l][1][1]-(Px[l][0][1]+Px[l][1][1])*x0)/w2[1];
       py2[1]=(2*Py[l][1][1]-(Py[l][0][1]+Py[l][1][1])*x0)/w2[1];
       pz2[1]=(2*Pz[l][1][1]-(Pz[l][0][1]+Pz[l][1][1])*x0)/w2[1];
       storeParticle(D,particle,i,j,k,s,x2,y2,w2,px2,py2,pz2,dx,dx,l,0,resolX,1.0,shift);
     } else;
   }

   for(l=0; l<cnt; l++) {
     for(m=0; m<2; m++) {
       free(Q[l][m]);       free(Px[l][m]);
       free(Py[l][m]);      free(Pz[l][m]);
     }
     free(Q[l]);       free(Px[l]);
     free(Py[l]);      free(Pz[l]);
   }
   free(Q);
   free(Px);
   free(Py);
   free(Pz);
}

void calInterpolation(double ***data,int cnt,double resolX)
{
   int i,j;
   double y1,y2,y3,a,b,c,dx,x;

   for(i=0; i<cnt; i++)    {
     for(j=0; j<2; j++)       {
       y1=data[i][j][0];
       y2=data[i][j][1];
       y3=data[i][j][2];
       a=0.5*(y1+y3)-y2;
       b=2.0*y2-1.5*y1-0.5*y3;
       c=y1;
       dx=1.0/resolX;
       x=0.5+dx*0.5;
       data[i][j][1]=a*x*x+b*x+c;
     }
   }
}	//lala

void deleteParticle(Particle ***particle,int i,int j,int k,int s)
{
  ptclList *p,*tmp;

  p=particle[i][j][k].head[s]->pt;
  while(p)  {
    tmp=p->next;
    particle[i][j][k].head[s]->pt=tmp;
    p->next=NULL;
    free(p);
    p=particle[i][j][k].head[s]->pt;
  }
/*
  p=D->befoParticle[i][j][k].head[s]->pt;
  while(p)  {
    tmp=p->next;
    D->befoParticle[i][j][k].head[s]->pt=tmp;
    p->next=NULL;
    free(p);
    p=D->befoParticle[i][j][k].head[s]->pt;
  }

  p=D->nextParticle[i][j][k].head[s]->pt;
  while(p)  {
    tmp=p->next;
    D->nextParticle[i][j][k].head[s]->pt=tmp;
    p->next=NULL;
    free(p);
    p=D->nextParticle[i][j][k].head[s]->pt;
  }
*/
}

void storeParticle(Domain *D,Particle ***particle,int i,int j,int k,int s,double *x2,double *y2,double *w2,double *px2,double *py2,double *pz2,double dx,double dy,int l,int m,double resolX,double resolY,double shift)
{
  double x,y,weight;
  ptclList *new;

  if(w2[0]>0)  {
    x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX+shift;
    y=(m*dy+y2[0]*dy)+j-D->jstart+D->minYSub*resolY;
    if(x<D->reNx)  {

      new = (ptclList *)malloc(sizeof(ptclList ));
      new->next = particle[i][j][k].head[s]->pt;
      particle[i][j][k].head[s]->pt=new;

      weight=w2[0];
      new->x=x; new->y=y; new->z=0.0;
      new->weight=weight*resolX*resolY;
      new->px=px2[0]; new->py=py2[0]; new->pz=pz2[0];
    } else;
  }   else ;

  if(w2[1]>0)  {
    x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX+shift;
    y=(m*dy+y2[1]*dy)+j-D->jstart+D->minYSub*resolY;
    if(x<D->reNx)  {

      new = (ptclList *)malloc(sizeof(ptclList ));
      new->next = particle[i][j][k].head[s]->pt;
      particle[i][j][k].head[s]->pt=new;

      weight=w2[1];
      new->x=x; new->y=y; new->z=0.0;
      new->weight=weight*resolX*resolY;
      new->px=px2[1]; new->py=py2[1]; new->pz=pz2[1];
    } else ;
  }   else ;
  
}
 
