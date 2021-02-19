#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

double calBeta(double beta,double db,double refE,double *delE,double *p1,double *p2,double w1,double w2,double *sig,int l,double **min,double **max);
void energyConservation(newP **newPart,int i,int j,double *delE,double **min,double **max,int l);
void deleteParticle(Particle ***particle,int i,int j,int k,int s);
void storeParticle(Particle ***particle,int i,int j,int k,int s,newP **newPart,int l,int m);
void reLowCreateParticle1D(Domain *D,int i,int j,int k,int s);
void reLowCreateParticle2D(Domain *D,int i,int j,int k,int s);
void reHighCreateParticle2D(Domain *D,int i,int j,int k,int s);
void calInterpolation(double ***data,int cnt,double resolX);

void reCreateParticle(Domain D)
{
   int i,j,k,s,istart,iend,jstart,jend,kstart,kend;
   ptclList *p;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D.istart;    iend=D.reIend;
   jstart=D.jstart;    jend=D.reJend;
   kstart=D.kstart;    kend=D.reKend;

   switch (D.dimension)  {
   case 1 :
     j=k=0;
     for(i=istart; i<iend; i++)
       for(s=0; s<D.nSpecies; s++)  
         reLowCreateParticle1D(&D,i,j,k,s);
     break;
   case 2 :
     k=0;
     if(D.mode==LOW) {
       for(i=istart; i<iend; i++) {
         for(j=jstart; j<jend; j++) 
           for(s=0; s<D.nSpecies; s++)  {
             p=D.particle[i][j][k].head[s]->pt;
             if(p) reLowCreateParticle2D(&D,i,j,k,s); else ;
           }
         if(i%10==0 && myrank==0) 
           printf("%d/%d is done\n",i,iend-istart);
         else;
       }
     } else if (D.mode==HIGH)  {
       for(i=istart; i<iend; i++) {
         for(j=jstart; j<jend; j++) 
           for(s=0; s<D.nSpecies; s++)  {
             p=D.particle[i][j][k].head[s]->pt;
             if(p) reHighCreateParticle2D(&D,i,j,k,s); else ;
           }
         if(i%10==0 && myrank==0) 
           printf("%d/%d is done\n",i,iend-istart);
         else;
       }
     }
     break;
   }

}

void reHighCreateParticle2D(Domain *D,int i,int j,int k,int s)
{
   int l,m,n,cntX,cntY,flag,indexI,indexJ;
   double gapX,W,x,y,px,py,pz,weight,tmp1,tmp2,Pr,Pl;
   double minX,maxX,minY,maxY,dx,dy,tmpX,tmpY,x0,y0;
   double resolX,resolY;
   double ***Q1,***Px1,***Py1,***Pz1;
   double px2[2],py2[2],pz2[2],x2[2],y2[2],w2[2],w[4],pp[3];
   size_t N,stride;
   Particle ***particle;
   newP **newPart1;
   particle=D->particle;
   ptclList *p;

   tmp1=sqrt(W/D->targetW[s]);
   tmp2=tmp1-(int)tmp1;
   if(tmp1<1.0) { cntX=1;         cntY=1;      }
   else {
     if(tmp2>=0.5) { cntX=(int)tmp1; cntY=cntX+1; }
     else          { cntX=(int)tmp1; cntY=cntX;   }
   }
   dx=1.0/((double)cntX);
   dy=1.0/((double)cntY);

   Q1=(double ***)malloc(cntX*sizeof(double **));
   Px1=(double ***)malloc(cntX*sizeof(double **));
   Py1=(double ***)malloc(cntX*sizeof(double **));
   Pz1=(double ***)malloc(cntX*sizeof(double **));
   for(l=0; l<cntX; l++) {
     Q1[l]=(double **)malloc(cntY*sizeof(double *));
     Px1[l]=(double **)malloc(cntY*sizeof(double *));
     Py1[l]=(double **)malloc(cntY*sizeof(double *));
     Pz1[l]=(double **)malloc(cntY*sizeof(double *));
     for(m=0; m<cntY; m++) {
       Q1[l][m]=(double *)malloc(4*sizeof(double *));
       Px1[l][m]=(double *)malloc(4*sizeof(double *));
       Py1[l][m]=(double *)malloc(4*sizeof(double *));
       Pz1[l][m]=(double *)malloc(4*sizeof(double *));
     }
   }
   for(l=0; l<cntX; l++)
     for(m=0; m<cntY; m++)
       for(n=0; n<4; n++) {
         Q1[l][m][n]=0.0; Px1[l][m][n]=0.0;
         Py1[l][m][n]=0.0; Pz1[l][m][n]=0.0;
       }

   p=particle[i][j][k].head[s]->pt;
   while(p)  {
     x=p->x; y=p->y; weight=p->weight;
     indexI=(int)(x/dx);  indexJ=(int)(y/dy);
     minX=indexI*dx;  maxX=(indexI+1)*dx;
     minY=indexJ*dy;  maxY=(indexJ+1)*dy;
     pp[0]=px=p->px; pp[1]=py=p->py; pp[2]=pz=p->pz;
     tmp2=weight/dx/dy;
     w[0]=(maxX-x)*(maxY-y);
     w[1]=(x-minX)*(maxY-y);
     w[2]=(maxX-x)*(y-minY);
     w[3]=(x-minX)*(y-minY);

     for(n=0; n<4; n++)  {
       Q1[indexI][indexJ][n]+=w[n]*tmp2;    
       Px1[indexI][indexJ][n]+=w[n]*tmp2*px;    
       Py1[indexI][indexJ][n]+=w[n]*tmp2*py;    
       Pz1[indexI][indexJ][n]+=w[n]*tmp2*pz;    
     }
     p=p->next;
   }
   deleteParticle(particle,i,j,k,s);


   resolX=D->resolX; resolY=D->resolY;

   newPart1=(newP **)malloc(cntX*sizeof(newP *));
   for(l=0; l<cntX; l++)
     newPart1[l]=(newP *)malloc(cntY*sizeof(newP ));

   for(l=0; l<cntX; l++)  
     for(m=0; m<cntY; m++)  {
       W=0.0;
       for(n=0; n<4; n++)  W+=Q1[l][m][n];

       if(W>0)  {
         x0=(Q1[l][m][1]+Q1[l][m][3])/W;
         y0=(Q1[l][m][2]+Q1[l][m][3])/W;
         if(x0==0.0) x0=1e-4; else ;
         if(y0==0.0) y0=1e-4; else ;
         tmpX=(x0-0.5)*(x0-0.5);
         tmpY=(y0-0.5)*(y0-0.5);

         if(tmpY<tmpX) {
           w2[0]=W*(1.0-y0);      w2[1]=W*y0;
           x2[0]=x0;          x2[1]=x0;
           y2[0]=y0*0.5;     y2[1]=(1.0+y0)*0.5;
           Pr=Px1[l][m][2]+Px1[l][m][3]; Pl=Px1[l][m][0]+Px1[l][m][1];
           px2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Py1[l][m][2]+Py1[l][m][3]; Pl=Py1[l][m][0]+Py1[l][m][1];
           py2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Pz1[l][m][2]+Pz1[l][m][3]; Pl=Pz1[l][m][0]+Pz1[l][m][1];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
         } else {
           w2[0]=W*(1.0-x0);      w2[1]=W*x0;
           x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
           y2[0]=y0;              y2[1]=y0;
           Pr=Px1[l][m][1]+Px1[l][m][3]; Pl=Px1[l][m][0]+Px1[l][m][2];
           px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Py1[l][m][1]+Py1[l][m][3]; Pl=Py1[l][m][0]+Py1[l][m][2];
           py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Pz1[l][m][1]+Pz1[l][m][3]; Pl=Pz1[l][m][0]+Pz1[l][m][2];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         }

         x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[0]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart1[l][m].x1=x;
         newPart1[l][m].y1=y;
         x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[1]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart1[l][m].x2=x;
         newPart1[l][m].y2=y;
         newPart1[l][m].px1=px2[0];
         newPart1[l][m].px2=px2[1];
         newPart1[l][m].py1=py2[0];
         newPart1[l][m].py2=py2[1];
         newPart1[l][m].pz1=pz2[0];
         newPart1[l][m].pz2=pz2[1];
         newPart1[l][m].w1=w2[0]*resolX*resolY;
         newPart1[l][m].w2=w2[1]*resolX*resolY;
         newPart1[l][m].flag=1;
       } else
         newPart1[l][m].flag=0;
     }	//Enf for l,m

   for(l=0; l<cntX; l++)
     for(m=0; m<cntY; m++)  {
       if(newPart1[l][m].flag==1) storeParticle(particle,i,j,k,s,newPart1,l,m);
       else ;
     }


   for(l=0; l<cntX; l++) {
     for(m=0; m<cntY;m++) {
       free(Q1[l][m]);       free(Px1[l][m]);
       free(Py1[l][m]);      free(Pz1[l][m]);
     }
     free(Q1[l]); free(Px1[l]); free(Py1[l]); free(Pz1[l]);
     free(newPart1[l]); 
   }
   free(Q1); free(Px1); free(Py1); free(Pz1);
   free(newPart1);
}
void reLowCreateParticle2D(Domain *D,int i,int j,int k,int s)
{
   int l,m,n,cntX,cntY,flag,indexI,indexJ,ind1,ind2,divMode=0,num;
   double gapX,W,x,y,px,py,pz,weight,tmp1,tmp2,Pr,Pl,maxV1,maxV2;
   double minX,maxX,minY,maxY,dx,dy,tmpX,tmpY,x0,y0,dh;
   double resolX,resolY,gamma,err;
   double originE[2],newE[2],delE;
   double minP[3],maxP[3],rngP[3],histo[9],pp0[3],pp[3];
   double ***Q1,***Px1,***Py1,***Pz1,***Q2,***Px2,***Py2,***Pz2;
   double px2[2],py2[2],pz2[2],x2[2],y2[2],w2[2],w[4],**min,**max;
   size_t N,stride;
   Particle ***particle;
   newP **newPart1,**newPart2;
   particle=D->particle;
   ptclList *p;

   W=0.0;
   pp0[0]=pp0[1]=pp0[2]=0.0;
   p=particle[i][j][k].head[s]->pt;
   minP[0]=1e5; maxP[0]=-1e5;
   minP[1]=1e5; maxP[1]=-1e5;
   minP[2]=1e5; maxP[2]=-1e5;
   while(p)  {
     weight=p->weight;
     px=p->px;     py=p->py;         pz=p->pz;
     W+=weight;
     pp0[0]+=weight*px; pp0[1]+=weight*py; pp0[2]+=weight*pz;
     if(px<minP[0]) minP[0]=px; else;
     if(px>maxP[0]) maxP[0]=px; else;
     if(py<minP[1]) minP[1]=py; else;
     if(py>maxP[1]) maxP[1]=py; else;
     if(pz<minP[2]) minP[2]=pz; else;
     if(pz>maxP[2]) maxP[2]=pz; else;
     p=p->next;
   }

   if(W==0.0) { printf("Weight is zero. i=%d,j=%d,k=%d\n",i,j,k); exit(0); }
   else {
     for(l=0; l<3; l++) {
       pp0[l]/=W;
       rngP[l]=maxP[l]-minP[l];
     }
     if(rngP[0]>=rngP[1] && rngP[0]>=rngP[2]) divMode=0;
     else if(rngP[1]>=rngP[0] && rngP[1]>=rngP[2]) divMode=1;
     else if(rngP[2]>=rngP[0] && rngP[2]>=rngP[1]) divMode=2;
     else divMode=1;

     if(rngP[divMode]==0.0)  {
       flag=0;
     } else {
       dh=rngP[divMode]/8.0;
       for(l=0; l<9; l++) histo[l]=0.0;
       p=particle[i][j][k].head[s]->pt;
       while(p)  {
         weight=p->weight;
         pp[0]=p->px; pp[1]=p->py; pp[2]=p->pz;
         l=(int)((pp[divMode]-minP[divMode])/dh);
         histo[l]+=weight;
         p=p->next;
       }

       ind1=ind2=4;
       maxV1=maxV2=histo[ind1];
       for(l=1; l<5; l++) {
         if(maxV1<histo[4-l]) ind1=4-l; else;
         if(maxV2<histo[4+l]) ind2=4+l; else;
       }
       if(ind1<4 && ind2>4) flag=1; else flag=0;
     }
   }

   if(flag==0) pp0[divMode]=minP[divMode];
   else        pp0[divMode]=(ind1+ind2)*0.5*dh+minP[divMode];

   tmp1=sqrt(W/D->targetW[s]*0.5);
   tmp2=tmp1-(int)tmp1;
   if(tmp1<1.0) { cntX=1;         cntY=1;      }
   else {
     if(tmp2>=0.5) { cntX=(int)tmp1; cntY=cntX+1; }
     else          { cntX=(int)tmp1; cntY=cntX;   }
   }
   dx=1.0/((double)cntX);
   dy=1.0/((double)cntY);

   Q1=(double ***)malloc(cntX*sizeof(double **));
   Px1=(double ***)malloc(cntX*sizeof(double **));
   Py1=(double ***)malloc(cntX*sizeof(double **));
   Pz1=(double ***)malloc(cntX*sizeof(double **));
   Q2=(double ***)malloc(cntX*sizeof(double **));
   Px2=(double ***)malloc(cntX*sizeof(double **));
   Py2=(double ***)malloc(cntX*sizeof(double **));
   Pz2=(double ***)malloc(cntX*sizeof(double **));
   for(l=0; l<cntX; l++) {
     Q1[l]=(double **)malloc(cntY*sizeof(double *));
     Px1[l]=(double **)malloc(cntY*sizeof(double *));
     Py1[l]=(double **)malloc(cntY*sizeof(double *));
     Pz1[l]=(double **)malloc(cntY*sizeof(double *));
     Q2[l]=(double **)malloc(cntY*sizeof(double *));
     Px2[l]=(double **)malloc(cntY*sizeof(double *));
     Py2[l]=(double **)malloc(cntY*sizeof(double *));
     Pz2[l]=(double **)malloc(cntY*sizeof(double *));
     for(m=0; m<cntY; m++) {
       Q1[l][m]=(double *)malloc(4*sizeof(double *));
       Px1[l][m]=(double *)malloc(4*sizeof(double *));
       Py1[l][m]=(double *)malloc(4*sizeof(double *));
       Pz1[l][m]=(double *)malloc(4*sizeof(double *));
       Q2[l][m]=(double *)malloc(4*sizeof(double *));
       Px2[l][m]=(double *)malloc(4*sizeof(double *));
       Py2[l][m]=(double *)malloc(4*sizeof(double *));
       Pz2[l][m]=(double *)malloc(4*sizeof(double *));
     }
   }
   for(l=0; l<cntX; l++)
     for(m=0; m<cntY; m++)
       for(n=0; n<4; n++) {
         Q1[l][m][n]=0.0; Px1[l][m][n]=0.0;
         Py1[l][m][n]=0.0; Pz1[l][m][n]=0.0;
         Q2[l][m][n]=0.0; Px2[l][m][n]=0.0;
         Py2[l][m][n]=0.0; Pz2[l][m][n]=0.0;
       }
   min=(double **)malloc(2*sizeof(double *));
   max=(double **)malloc(2*sizeof(double *));
   for(l=0; l<2; l++) {
     min[l]=(double *)malloc(3*sizeof(double ));
     max[l]=(double *)malloc(3*sizeof(double ));
   }
   for(l=0; l<2; l++)
     for(m=0; m<3; m++) {
       min[l][m]=1e5;
       max[l][m]=-1e5;
     }

   originE[0]=0.0;
   p=particle[i][j][k].head[s]->pt;
   while(p)  {
     x=p->x; y=p->y; weight=p->weight;
     indexI=(int)(x/dx);  indexJ=(int)(y/dy);
     minX=indexI*dx;  maxX=(indexI+1)*dx;
     minY=indexJ*dy;  maxY=(indexJ+1)*dy;
     pp[0]=px=p->px; pp[1]=py=p->py; pp[2]=pz=p->pz;
     gamma=sqrt(1.0+px*px+py*py+pz*pz);
     tmp2=weight/dx/dy;
     w[0]=(maxX-x)*(maxY-y);
     w[1]=(x-minX)*(maxY-y);
     w[2]=(maxX-x)*(y-minY);
     w[3]=(x-minX)*(y-minY);
     if(pp[divMode]<pp0[divMode]) {
       for(l=0; l<3; l++) {
         if(pp[l]<min[0][l]) min[0][l]=pp[l]; else;
         if(pp[l]>max[0][l]) max[0][l]=pp[l]; else;
       }
       for(n=0; n<4; n++)  {
         Q1[indexI][indexJ][n]+=w[n]*tmp2;    
         Px1[indexI][indexJ][n]+=w[n]*tmp2*px;    
         Py1[indexI][indexJ][n]+=w[n]*tmp2*py;    
         Pz1[indexI][indexJ][n]+=w[n]*tmp2*pz;    
       }
       originE[0]+=weight*gamma;
     } else {
       for(l=0; l<3; l++) {
         if(pp[l]<min[1][l]) min[1][l]=pp[l]; else;
         if(pp[l]>max[1][l]) max[1][l]=pp[l]; else;
       }
       for(n=0; n<4; n++)  {
         Q2[indexI][indexJ][n]+=w[n]*tmp2;    
         Px2[indexI][indexJ][n]+=w[n]*tmp2*px;    
         Py2[indexI][indexJ][n]+=w[n]*tmp2*py;    
         Pz2[indexI][indexJ][n]+=w[n]*tmp2*pz;    
       }
       originE[1]+=weight*gamma;
     }
     p=p->next;
   }
   deleteParticle(particle,i,j,k,s);


   resolX=1.0/((double)D->resolX);
   resolY=1.0/((double)D->resolY);

   //lower division
   newPart1=(newP **)malloc(cntX*sizeof(newP *));
   for(l=0; l<cntX; l++)
     newPart1[l]=(newP *)malloc(cntY*sizeof(newP ));

   newE[0]=0.0;
   for(l=0; l<cntX; l++)  
     for(m=0; m<cntY; m++)  {
       W=0.0;
       for(n=0; n<4; n++)  W+=Q1[l][m][n];

       if(W>0)  {
         x0=(Q1[l][m][1]+Q1[l][m][3])/W;
         y0=(Q1[l][m][2]+Q1[l][m][3])/W;
         if(x0==0.0) x0=1e-4; else ;
         if(y0==0.0) y0=1e-4; else ;
         tmpX=(x0-0.5)*(x0-0.5);
         tmpY=(y0-0.5)*(y0-0.5);

         if(tmpY<tmpX) {
           w2[0]=W*(1.0-y0);      w2[1]=W*y0;
           x2[0]=x0;          x2[1]=x0;
           y2[0]=y0*0.5;     y2[1]=(1.0+y0)*0.5;
           Pr=Px1[l][m][2]+Px1[l][m][3]; Pl=Px1[l][m][0]+Px1[l][m][1];
           px2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Py1[l][m][2]+Py1[l][m][3]; Pl=Py1[l][m][0]+Py1[l][m][1];
           py2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Pz1[l][m][2]+Pz1[l][m][3]; Pl=Pz1[l][m][0]+Pz1[l][m][1];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
         } else {
           w2[0]=W*(1.0-x0);      w2[1]=W*x0;
           x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
           y2[0]=y0;              y2[1]=y0;
           Pr=Px1[l][m][1]+Px1[l][m][3]; Pl=Px1[l][m][0]+Px1[l][m][2];
           px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Py1[l][m][1]+Py1[l][m][3]; Pl=Py1[l][m][0]+Py1[l][m][2];
           py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Pz1[l][m][1]+Pz1[l][m][3]; Pl=Pz1[l][m][0]+Pz1[l][m][2];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         }

         x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[0]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart1[l][m].x1=x;
         newPart1[l][m].y1=y;
         x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[1]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart1[l][m].x2=x;
         newPart1[l][m].y2=y;
         newPart1[l][m].px1=px2[0];
         newPart1[l][m].px2=px2[1];
         newPart1[l][m].py1=py2[0];
         newPart1[l][m].py2=py2[1];
         newPart1[l][m].pz1=pz2[0];
         newPart1[l][m].pz2=pz2[1];
         newPart1[l][m].w1=w2[0]*resolX*resolY;
         newPart1[l][m].w2=w2[1]*resolX*resolY;
         newPart1[l][m].flag=1;
         for(n=0; n<2; n++)
           newE[0]+=w2[n]*sqrt(1.0+px2[n]*px2[n]+py2[n]*py2[n]+pz2[n]*pz2[n]);
       } else
         newPart1[l][m].flag=0;
     }	//Enf for l,m

   //upper division
   newPart2=(newP **)malloc(cntX*sizeof(newP *));
   for(l=0; l<cntX; l++)
     newPart2[l]=(newP *)malloc(cntY*sizeof(newP ));
 
   newE[1]=0.0;
   for(l=0; l<cntX; l++)  
     for(m=0; m<cntY; m++)  {
       W=0.0;
       for(n=0; n<4; n++)  W+=Q2[l][m][n];

       if(W>0)  {
         x0=(Q2[l][m][1]+Q2[l][m][3])/W;
         y0=(Q2[l][m][2]+Q2[l][m][3])/W;
         if(x0==0.0) x0=1e-4; else ;
         if(y0==0.0) y0=1e-4; else ;
         tmpX=(x0-0.5)*(x0-0.5);
         tmpY=(y0-0.5)*(y0-0.5);

         if(tmpY<tmpX) {
           w2[0]=W*(1.0-y0);      w2[1]=W*y0;
           x2[0]=x0;          x2[1]=x0;
           y2[0]=y0*0.5;     y2[1]=(1.0+y0)*0.5;
           Pr=Px2[l][m][2]+Px2[l][m][3]; Pl=Px2[l][m][0]+Px2[l][m][1];
           px2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Py2[l][m][2]+Py2[l][m][3]; Pl=Py2[l][m][0]+Py2[l][m][1];
           py2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
           Pr=Pz2[l][m][2]+Pz2[l][m][3]; Pl=Pz2[l][m][0]+Pz2[l][m][1];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*y0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*y0)/w2[1];
         } else {
           w2[0]=W*(1.0-x0);      w2[1]=W*x0;
           x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
           y2[0]=y0;              y2[1]=y0;
           Pr=Px2[l][m][1]+Px2[l][m][3]; Pl=Px2[l][m][0]+Px2[l][m][2];
           px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Py2[l][m][1]+Py2[l][m][3]; Pl=Py2[l][m][0]+Py2[l][m][2];
           py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
           Pr=Pz2[l][m][1]+Pz2[l][m][3]; Pl=Pz2[l][m][0]+Pz2[l][m][2];
           pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
           pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         }

         x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[0]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart2[l][m].x1=x;
         newPart2[l][m].y1=y;
         x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX;
         y=(m*dy+y2[1]*dy)+j-D->jstart+D->minYSub*resolY;
         newPart2[l][m].x2=x;
         newPart2[l][m].y2=y;
         newPart2[l][m].px1=px2[0];
         newPart2[l][m].px2=px2[1];
         newPart2[l][m].py1=py2[0];
         newPart2[l][m].py2=py2[1];
         newPart2[l][m].pz1=pz2[0];
         newPart2[l][m].pz2=pz2[1];
         newPart2[l][m].w1=w2[0]*resolX*resolY;
         newPart2[l][m].w2=w2[1]*resolX*resolY;
         newPart2[l][m].flag=1;
         for(n=0; n<2; n++)
           newE[1]+=w2[n]*sqrt(1.0+px2[n]*px2[n]+py2[n]*py2[n]+pz2[n]*pz2[n]);
       } else
         newPart2[l][m].flag=0;
     }	//Enf for l,m
/*
   delE=sqrt(delE*delE);
   err=delE/originE[0];
   num=0;
   while(err>0.01 && num<10) {
     for(l=0; l<cntX; l+=2)
       for(m=0; m<cntY; m+=2)  {
         if(newPart1[l][m].flag==1) 
           energyConservation(newPart1,l,m,&delE,min,max,0);
         else ;
         err=delE/originE[0];         
         if(err<0.01) { l=cntX; m=cntY; } else ;
       }
     for(l=1; l<cntX; l+=2)
       for(m=1; m<cntY; m+=2)  {
         if(newPart1[l][m].flag==1) 
           energyConservation(newPart1,l,m,&delE,min,max,0);
         else ;
         err=delE/originE[0];         
         if(err<0.01) { l=cntX; m=cntY; } else ;
       }
     err=delE/originE[0];
     num++;
   }

   delE=originE[1]-newE[1];
   delE=sqrt(delE*delE);
   err=delE/originE[1];
   num=0;
   while(err>0.01 && num<10) {
     for(l=0; l<cntX; l+=2)
       for(m=0; m<cntY; m+=2)  {
         if(newPart2[l][m].flag==1) 
           energyConservation(newPart2,l,m,&delE,min,max,0);
         else ;
         err=delE/originE[1];         
         if(err<0.01) { l=cntX; m=cntY; } else ;
       }
     for(l=1; l<cntX; l+=2)
       for(m=1; m<cntY; m+=2)  {
         if(newPart2[l][m].flag==1) 
           energyConservation(newPart2,l,m,&delE,min,max,0);
         else ;
         err=delE/originE[1];         
         if(err<0.01) { l=cntX; m=cntY; } else ;
       }
     err=delE/originE[1];
     num++;
   }
*/

   for(l=0; l<cntX; l++)
     for(m=0; m<cntY; m++)  {
       if(newPart1[l][m].flag==1) storeParticle(particle,i,j,k,s,newPart1,l,m);
       else ;
     }
   for(l=0; l<cntX; l++)
     for(m=0; m<cntY; m++)  {
       if(newPart2[l][m].flag==1) storeParticle(particle,i,j,k,s,newPart2,l,m);
       else ;
     }



   for(l=0; l<cntX; l++) {
     for(m=0; m<cntY;m++) {
       free(Q1[l][m]);       free(Px1[l][m]);
       free(Py1[l][m]);      free(Pz1[l][m]);
       free(Q2[l][m]);       free(Px2[l][m]);
       free(Py2[l][m]);      free(Pz2[l][m]);
     }
     free(Q1[l]); free(Px1[l]); free(Py1[l]); free(Pz1[l]);
     free(Q2[l]); free(Px2[l]); free(Py2[l]); free(Pz2[l]);
     free(newPart1[l]); free(newPart2[l]);
   }
   free(Q1); free(Px1); free(Py1); free(Pz1);
   free(Q2); free(Px2); free(Py2); free(Pz2);
   free(newPart1); free(newPart2);
   for(l=0; l<2; l++) { free(min[l]); free(max[l]); }
   free(min); free(max);

}

void reLowCreateParticle1D(Domain *D,int i,int j,int k,int s)
{
   int l,m,n,cnt,flag,indexI,indexJ,ind1,ind2,divMode=0,num;
   double gapX,W,x,y,px,py,pz,weight,tmp2,Pr,Pl,maxV1,maxV2,gamma,err;
   double minX,maxX,minY,maxY,dx,dy,tmpX,tmpY,x0,y0,dh;
   double resolX,resolY,newE[2],originE[2];
   double minP[3],maxP[3],rngP[3],histo[9],pp0[3],pp[3];
   double ***Q1,***Px1,***Py1,***Pz1,***Q2,***Px2,***Py2,***Pz2;
   double px2[2],py2[2],pz2[2],x2[2],y2[2],w2[2],w[2],**min,**max;
   size_t N,stride;
   Particle ***particle;
   newP **newPart1,**newPart2;
   particle=D->particle;
   ptclList *p;

   W=0.0;
   pp0[0]=pp0[1]=pp0[2]=0.0;
   p=particle[i][j][k].head[s]->pt;
   minP[0]=1e5; maxP[0]=-1e5;
   minP[1]=1e5; maxP[1]=-1e5;
   minP[2]=1e5; maxP[2]=-1e5;
   while(p)  {
     weight=p->weight;
     px=p->px;     py=p->py;         pz=p->pz;
     W+=weight;
     pp0[0]+=weight*px; pp0[1]+=weight*py; pp0[2]+=weight*pz;
     if(px<minP[0]) minP[0]=px; else;
     if(px>maxP[0]) maxP[0]=px; else;
     if(py<minP[1]) minP[1]=py; else;
     if(py>maxP[1]) maxP[1]=py; else;
     if(pz<minP[2]) minP[2]=pz; else;
     if(pz>maxP[2]) maxP[2]=pz; else;
     p=p->next;
   }

   if(W>0) {
     for(l=0; l<3; l++) {
       pp0[l]/=W;
       rngP[l]=maxP[l]-minP[l];
     }
     if(rngP[0]>=rngP[1] && rngP[0]>=rngP[2]) divMode=0;
     else if(rngP[1]>=rngP[0] && rngP[1]>=rngP[2]) divMode=1;
     else if(rngP[2]>=rngP[0] && rngP[2]>=rngP[1]) divMode=2;
     else divMode=1;

     if(rngP[divMode]==0.0)  {
       flag=0;
     } else {
       dh=rngP[divMode]/8.0;
       for(l=0; l<9; l++) histo[l]=0.0;
       p=particle[i][j][k].head[s]->pt;
       while(p)  {
         weight=p->weight;
         pp[0]=p->px; pp[1]=p->py; pp[2]=p->pz;
         l=(int)((pp[divMode]-minP[divMode])/dh);
         histo[l]+=weight;
         p=p->next;
       }

       ind1=ind2=4;
       maxV1=maxV2=histo[ind1];
       for(l=1; l<5; l++) {
         if(maxV1<histo[4-l]) ind1=4-l; else;
         if(maxV2<histo[4+l]) ind2=4+l; else;
       }
       if(ind1<4 && ind2>4) flag=1; else flag=0;
     }
   }

   if(flag==0) pp0[divMode]=minP[divMode];
   else        pp0[divMode]=(ind1+ind2)*0.5*dh+minP[divMode];

   tmp2=D->targetW[s];
   cnt=(int)(1.0/tmp2);
   dx=1.0/((double)cnt);

   Q1=(double ***)malloc(cnt*sizeof(double **));
   Px1=(double ***)malloc(cnt*sizeof(double **));
   Py1=(double ***)malloc(cnt*sizeof(double **));
   Pz1=(double ***)malloc(cnt*sizeof(double **));
   Q2=(double ***)malloc(cnt*sizeof(double **));
   Px2=(double ***)malloc(cnt*sizeof(double **));
   Py2=(double ***)malloc(cnt*sizeof(double **));
   Pz2=(double ***)malloc(cnt*sizeof(double **));
   for(l=0; l<cnt; l++) {
     Q1[l]=(double **)malloc(1*sizeof(double *));
     Px1[l]=(double **)malloc(1*sizeof(double *));
     Py1[l]=(double **)malloc(1*sizeof(double *));
     Pz1[l]=(double **)malloc(1*sizeof(double *));
     Q2[l]=(double **)malloc(1*sizeof(double *));
     Px2[l]=(double **)malloc(1*sizeof(double *));
     Py2[l]=(double **)malloc(1*sizeof(double *));
     Pz2[l]=(double **)malloc(1*sizeof(double *));
     for(m=0; m<1; m++) {
       Q1[l][m]=(double *)malloc(2*sizeof(double *));
       Px1[l][m]=(double *)malloc(2*sizeof(double *));
       Py1[l][m]=(double *)malloc(2*sizeof(double *));
       Pz1[l][m]=(double *)malloc(2*sizeof(double *));
       Q2[l][m]=(double *)malloc(2*sizeof(double *));
       Px2[l][m]=(double *)malloc(2*sizeof(double *));
       Py2[l][m]=(double *)malloc(2*sizeof(double *));
       Pz2[l][m]=(double *)malloc(2*sizeof(double *));
     }
   }
   for(l=0; l<cnt; l++)
     for(m=0; m<1; m++)
       for(n=0; n<2; n++) {
         Q1[l][m][n]=0.0; Px1[l][m][n]=0.0;
         Py1[l][m][n]=0.0; Pz1[l][m][n]=0.0;
         Q2[l][m][n]=0.0; Px2[l][m][n]=0.0;
         Py2[l][m][n]=0.0; Pz2[l][m][n]=0.0;
       }
   min=(double **)malloc(2*sizeof(double *));
   max=(double **)malloc(2*sizeof(double *));
   for(l=0; l<2; l++) {
     min[l]=(double *)malloc(3*sizeof(double ));
     max[l]=(double *)malloc(3*sizeof(double ));
   }
   for(l=0; l<2; l++)
     for(m=0; m<3; m++) {
       min[l][m]=1e5;
       max[l][m]=-1e5;
     }


   p=particle[i][j][k].head[s]->pt;
   while(p)  {
     x=p->x; y=p->y; weight=p->weight;
     indexI=(int)(x/dx); indexJ=0;
     minX=indexI*dx;  maxX=(indexI+1)*dx;
     pp[0]=px=p->px; pp[1]=py=p->py; pp[2]=pz=p->pz;
     gamma=sqrt(1.0+px*px+py*py+pz*pz);
     tmp2=weight/dx;
     w[0]=(maxX-x);
     w[1]=(x-minX);
     if(pp[divMode]<pp0[divMode]) {
       for(l=0; l<3; l++) {
         if(pp[l]<min[0][l]) min[0][l]=pp[l]; else;
         if(pp[l]>max[0][l]) max[0][l]=pp[l]; else;
       }
       for(n=0; n<2; n++)  {
         Q1[indexI][indexJ][n]+=w[n]*tmp2;    
         Px1[indexI][indexJ][n]+=w[n]*tmp2*px;    
         Py1[indexI][indexJ][n]+=w[n]*tmp2*py;    
         Pz1[indexI][indexJ][n]+=w[n]*tmp2*pz;    
       }
       originE[0]+=weight*gamma;
     } else {
       for(l=0; l<3; l++) {
         if(pp[l]<min[1][l]) min[1][l]=pp[l]; else;
         if(pp[l]>max[1][l]) max[1][l]=pp[l]; else;
       }
       for(n=0; n<2; n++)  {
         Q2[indexI][indexJ][n]+=w[n]*tmp2;    
         Px2[indexI][indexJ][n]+=w[n]*tmp2*px;    
         Py2[indexI][indexJ][n]+=w[n]*tmp2*py;    
         Pz2[indexI][indexJ][n]+=w[n]*tmp2*pz;    
       }
       originE[1]+=weight*gamma;
     }
     p=p->next;
   }
   deleteParticle(particle,i,j,k,s);


     resolX=1.0/((double)D->resolX);
     resolY=1.0/((double)D->resolY);

   //lower division
   newPart1=(newP **)malloc(cnt*sizeof(newP *));
   for(l=0; l<cnt; l++)
     newPart1[l]=(newP *)malloc(1*sizeof(newP ));

   newE[0]=0.0;
   for(l=0; l<cnt; l++)  
     for(m=0; m<1; m++)  {
       W=0.0;
       for(n=0; n<2; n++)  W+=Q1[l][m][n];

       if(W>0)  {
         x0=Q1[l][m][1]/W;
         w2[0]=Q1[l][m][0];     w2[1]=Q1[l][m][1]; 
         x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
         Pr=Px1[l][m][1]; Pl=Px1[l][m][0];
         px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         Pr=Py1[l][m][1]; Pl=Py1[l][m][0];
         py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         Pr=Pz1[l][m][1]; Pl=Pz1[l][m][0];
         pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];

         x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX;
         newPart1[l][m].x1=x;
         newPart1[l][m].y1=0.0;
         x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX;
         newPart1[l][m].x2=x;
         newPart1[l][m].y2=0.0;
         newPart1[l][m].px1=px2[0];
         newPart1[l][m].px2=px2[1];
         newPart1[l][m].py1=py2[0];
         newPart1[l][m].py2=py2[1];
         newPart1[l][m].pz1=pz2[0];
         newPart1[l][m].pz2=pz2[1];
         newPart1[l][m].w1=w2[0]*resolX;
         newPart1[l][m].w2=w2[1]*resolX;
         newPart1[l][m].flag=1;
         for(n=0; n<2; n++)
           newE[0]+=w2[n]*sqrt(1.0+px2[n]*px2[n]+py2[n]*py2[n]+pz2[n]*pz2[n]);
       } else
         newPart1[l][m].flag=0;
     }	//Enf for l,m

   //upper division
   newPart2=(newP **)malloc(cnt*sizeof(newP *));
   for(l=0; l<cnt; l++)
     newPart2[l]=(newP *)malloc(1*sizeof(newP ));
 
   newE[1]=0.0;
   for(l=0; l<cnt; l++)  
     for(m=0; m<1; m++)  {
       W=0.0;
       for(n=0; n<2; n++)  W+=Q2[l][m][n];

       if(W>0)  {
         x0=Q2[l][m][1]/W;
         w2[0]=Q2[l][m][0];     w2[1]=Q2[l][m][1]; 
         x2[0]=x0*0.5;          x2[1]=(1.0+x0)*0.5;
         Pr=Px2[l][m][1]; Pl=Px2[l][m][0];
         px2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         px2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         Pr=Py2[l][m][1]; Pl=Py2[l][m][0];
         py2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         py2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];
         Pr=Pz2[l][m][1]; Pl=Pz2[l][m][0];
         pz2[0]=((Pl-Pr)+(Pr+Pl)*x0)/w2[0];
         pz2[1]=(2*Pr-(Pr+Pl)*x0)/w2[1];

         x=(l*dx+x2[0]*dx)+i-D->istart+D->minXSub*resolX;
         newPart2[l][m].x1=x;
         newPart2[l][m].y1=0.0;
         x=(l*dx+x2[1]*dx)+i-D->istart+D->minXSub*resolX;
         newPart2[l][m].x2=x;
         newPart2[l][m].y2=0.0;
         newPart2[l][m].px1=px2[0];
         newPart2[l][m].px2=px2[1];
         newPart2[l][m].py1=py2[0];
         newPart2[l][m].py2=py2[1];
         newPart2[l][m].pz1=pz2[0];
         newPart2[l][m].pz2=pz2[1];
         newPart2[l][m].w1=w2[0]*resolX;
         newPart2[l][m].w2=w2[1]*resolX;
         newPart2[l][m].flag=1;
         for(n=0; n<2; n++)
           newE[1]+=w2[n]*sqrt(1.0+px2[n]*px2[n]+py2[n]*py2[n]+pz2[n]*pz2[n]);
       } else
         newPart2[l][m].flag=0;
     }	//Enf for l,m

/*
   delE=originE[0]-newE[0];
   delE=sqrt(delE*delE);
   err=delE/originE[0];
   num=0;
   while(err>0.01 && num<10) {
     for(l=0; l<cnt; l+=2)
       for(m=0; m<cnt; m+=2)  {
         if(newPart1[l][m].flag==1) 
           energyConservation(newPart1,l,m,&delE,min,max,0);
         else ;
       }
     for(l=1; l<cnt; l+=2)
       for(m=1; m<cnt; m+=2)  {
         if(newPart1[l][m].flag==1) 
           energyConservation(newPart1,l,m,&delE,min,max,0);
         else ;
       }
     err=delE/originE[0];
     num++;
   }

   delE=originE[1]-newE[1];
   delE=sqrt(delE*delE);
   err=delE/originE[1];
   num=0;
   while(err>0.01 && num<10) {
     for(l=0; l<cnt; l+=2)
       for(m=0; m<cnt; m+=2)  {
         if(newPart2[l][m].flag==1) 
           energyConservation(newPart2,l,m,&delE,min,max,0);
         else ;
       }
     for(l=1; l<cnt; l+=2)
       for(m=1; m<cnt; m+=2)  {
         if(newPart2[l][m].flag==1) 
           energyConservation(newPart2,l,m,&delE,min,max,0);
         else ;
       }
     err=delE/originE[1];
     num++;
   }
*/

   for(l=0; l<cnt; l++)
     for(m=0; m<1; m++)  {
       if(newPart1[l][m].flag==1) storeParticle(particle,i,j,k,s,newPart1,l,m);
       else ;
     }
   for(l=0; l<cnt; l++)
     for(m=0; m<1; m++)  {
       if(newPart2[l][m].flag==1) storeParticle(particle,i,j,k,s,newPart2,l,m);
       else ;
     }



   for(l=0; l<cnt; l++) {
     for(m=0; m<1;m++) {
       free(Q1[l][m]);       free(Px1[l][m]);
       free(Py1[l][m]);      free(Pz1[l][m]);
       free(Q2[l][m]);       free(Px2[l][m]);
       free(Py2[l][m]);      free(Pz2[l][m]);
     }
     free(Q1[l]); free(Px1[l]); free(Py1[l]); free(Pz1[l]);
     free(Q2[l]); free(Px2[l]); free(Py2[l]); free(Pz2[l]);
     free(newPart1[l]); free(newPart2[l]);
   }
   free(Q1); free(Px1); free(Py1); free(Pz1);
   free(Q2); free(Px2); free(Py2); free(Pz2);
   free(newPart1); free(newPart2);
   for(l=0; l<2; l++) { free(min[l]); free(max[l]); }
   free(min); free(max);

}
/*
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
*/

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

void storeParticle(Particle ***particle,int i,int j,int k,int s,newP **newPart,int l,int m)
{
  ptclList *new;

  if(newPart[l][m].w1>0.0) {
  new = (ptclList *)malloc(sizeof(ptclList ));
  new->next = particle[i][j][k].head[s]->pt;
  particle[i][j][k].head[s]->pt=new;

  new->x=newPart[l][m].x1;
  new->y=newPart[l][m].y1;
  new->z=0.0;
  new->px=newPart[l][m].px1;
  new->py=newPart[l][m].py1;
  new->pz=newPart[l][m].pz1;
  new->weight=newPart[l][m].w1;
  } else ;

  if(newPart[l][m].w2>0.0) {
  new = (ptclList *)malloc(sizeof(ptclList ));
  new->next = particle[i][j][k].head[s]->pt;
  particle[i][j][k].head[s]->pt=new;

  new->x=newPart[l][m].x2;
  new->y=newPart[l][m].y2;
  new->z=0.0;
  new->px=newPart[l][m].px2;
  new->py=newPart[l][m].py2;
  new->pz=newPart[l][m].pz2;
  new->weight=newPart[l][m].w2;
  } else ;
}

void energyConservation(newP **newPart,int i,int j,double *delE,double **min,double **max,int l)
{
  int k,flag;
  double p1[3],p2[3],sig[3],w1,w2,gamma1,gamma2,sig1,sig2;
  double beta,dbeta,refE,tmp1,tmp2;

  w1=newPart[i][j].w1;
  p1[0]=newPart[i][j].px1;
  p1[1]=newPart[i][j].py1;
  p1[2]=newPart[i][j].pz1;
  gamma1=sqrt(1.0+p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);

  w2=newPart[i][j].w2;
  p2[0]=newPart[i][j].px2;
  p2[1]=newPart[i][j].py2;
  p2[2]=newPart[i][j].pz2;
  gamma2=sqrt(1.0+p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);

  for(k=0; k<3; k++) {
    if(p1[k]>p2[k]) {
      sig1=max[l][k]-p1[k]; sig2=min[l][k]-p2[k];
      if(sig1*sig1<sig2*sig2) sig[k]=sig1;
      else                    sig[k]=w2/w1*sig2;
    } else {
      sig1=min[l][k]-p1[k]; sig2=max[l][k]-p2[k];
      if(sig1*sig1<sig2*sig2) sig[k]=sig1;
      else                    sig[k]=w2/w1*sig2;
    }
  }

  refE=w1*gamma1+w2*gamma2;
  beta=calBeta(beta,dbeta,refE,delE,p1,p2,w1,w2,sig,l,min,max);

  newPart[i][j].px1=p1[0]+beta*sig[0];
  newPart[i][j].py1=p1[1]+beta*sig[1];
  newPart[i][j].pz1=p1[2]+beta*sig[2];
  newPart[i][j].px2=p2[0]-beta*sig[0]*w1/w2;
  newPart[i][j].py2=p2[1]-beta*sig[1]*w1/w2;
  newPart[i][j].pz2=p2[2]-beta*sig[2]*w1/w2;
}

double calBeta(double beta,double db,double refE,double *delE,double *p1,double *p2,double w1,double w2,double *sig,int l,double **min,double **max)
{
  int k,cnt,flag;
  double pp1[3],pp2[3],energy,dE,DelE;

  DelE=*delE;
  DelE=sqrt(DelE*DelE);
  cnt=0;
  flag=0;
  while(flag==0 )  {
    beta=randomValue(1.0);
    pp1[0]=p1[0]+beta*sig[0];
    pp1[1]=p1[1]+beta*sig[1];
    pp1[2]=p1[2]+beta*sig[2];
    energy=w1*sqrt(1.0+pp1[0]*pp1[0]+pp1[1]*pp1[1]+pp1[2]*pp1[2]);
    pp2[0]=p2[0]-w1/w2*beta*sig[0];
    pp2[1]=p2[1]-w1/w2*beta*sig[1];
    pp2[2]=p2[2]-w1/w2*beta*sig[2];
    energy+=w2*sqrt(1.0+pp2[0]*pp2[0]+pp2[1]*pp2[1]+pp2[2]*pp2[2]);
    dE=refE-energy;
    dE=sqrt(dE*dE);
    if(DelE>dE || dE==0.0) flag=1; else flag=0;
    cnt++;
  }

  flag=1;
  for(k=0; k<3; k++) {
    if(pp1[k]>max[l][k] || pp1[k]<min[l][k] ||
       pp2[k]>max[l][k] || pp2[k]<min[l][k]) {
       flag=0;
    } else ;
  }
//  if(cnt==10) flag=0; else ;

  if(flag==0) { beta=0.0; }
  else  *delE=DelE-dE;
  return beta;
}

