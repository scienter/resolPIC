#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


double maxwellianVelocity(double temperature);
void loadPlasma_crystal(Domain *D,LoadList *LL,int s);
void random1D_sobol(double *x,gsl_qrng *q1);
void random2D_sobol(double *x,double *y,gsl_qrng *q2);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q3);
void loadDefinedPlasma2D(Domain *D,LoadList *LL,int s);
void loadDefinedPlasma3D(Domain *D,LoadList *LL,int s);
void loadBeamPlasma1D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeamPlasma2D(Domain *D,LoadList *LL,int s,int iteration);
double gaussian_dist(double sig);


double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-(x-centerX)*(x-centerX)/gaussCoefX/gaussCoefX);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefX*(x-centerX)*(x-centerX);
    break;
  }
  return result;
}

double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ)
{
  double result;

  switch (mode)  {
  case 0 :	//Costant
    result=1.0;
    break;
  case 1 :	//Gaussian
    result=exp(-((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ))/gaussCoefYZ/gaussCoefYZ);
    break;
  case 2 :	//2nd polynomial
    result=1.0+polyCoefYZ*((y-centerY)*(y-centerY)+(z-centerZ)*(z-centerZ));
    break;
  }
  return result;
}

void loadPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
  void loadPolygonPlasma1D();
  void loadPolygonPlasma2D();
  void loadPolygonPlasma3D();
  void loadChannelPlasma();
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  

//  void loadPlasma_crystal(Domain *D,LoadList *LL,int s);


  switch((LL->type-1)*3+D->dimension)  {
  //1D
  case ((Polygon-1)*3+1):
    loadPolygonPlasma1D(D,LL,s,iteration); 
//    loadPlasma_crystal(D,LL,s);
    break;

  //2D
  case (Polygon-1)*3+2:
    loadPolygonPlasma2D(D,LL,s,iteration); 
    break;

  case (Defined-1)*3+2:
    loadDefinedPlasma2D(D,LL,s); 
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //3D
  case (Polygon-1)*3+3:
    loadPolygonPlasma3D(D,LL,s,iteration); 
    break;
  case (Defined-1)*3+3:
//    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime) {
    loadDefinedPlasma3D(D,LL,s); 
//    }
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  case (Beam-1)*3+1:
  case (Beam-1)*3+2:
  case (Beam-1)*3+3:
    ;
    break;
  default:
    ;
  }
}


void loadDefinedPlasma3D(Domain *D,LoadList *LL,int s)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,cnt;
   int l,n,myrank;
   double dx,dy,dz,v1,v2,v3,charge,weightCoef;
   double x,y,z,xPos,yPos,zPos;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New;
   DefPtcl *p;

   dx=D->dx;   dy=D->dy;   dz=D->dz;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   cnt=LL->numDefPtcls;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   //position define      
   p=LL->def;
   while(p) {
     cnt=p->numDefPtcls;
     for(n=0; n<cnt; n++)
     {
       xPos=p->define[n];
       i=(int)(xPos/dx-D->minXSub+istart);
       if(i<iend && i>=istart)
       {
         yPos=p->define[n+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         if(jstart<=j && j<jend)
         {
           zPos=p->define[n+2*cnt];
           k=(int)(zPos/dz-D->minZSub+kstart);
           if(kstart<=k && k<kend)
           {
             x=xPos/dx-D->minXSub+istart-i;
             y=yPos/dy-D->minYSub+jstart-j;
             z=zPos/dz-D->minZSub+kstart-k;
             New = (ptclList *)malloc(sizeof(ptclList));
             New->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = New;

             New->x = x;           New->oldX=i+x;
             New->y = y;           New->oldY=j+y;
             New->z = z;           New->oldZ=k+z;

             New->E1=New->E2=New->E3=0.0;
             New->B1=New->B2=New->B3=0.0;
             v1=maxwellianVelocity(LL->temperature)/velocityC;
             v2=maxwellianVelocity(LL->temperature)/velocityC;
             v3=maxwellianVelocity(LL->temperature)/velocityC;
             New->p1=-D->gamma*D->beta+v1;
             New->p2=v2;
             New->p3=v3;
             New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
             New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
             New->weight=1.0/LL->numberInCell*weightCoef;
             New->charge=charge;
             LL->index+=1;
             New->index=LL->index;
             New->core=myrank;
           }
         }
       }                //if(i,j,k<range)
     }                  //End of for(n<cnt)
     p=p->next;
   }                    //End of while(p)
}

void loadDefinedPlasma2D(Domain *D,LoadList *LL,int s)
{
   int ii,i,j,k,istart,iend,jstart,jend,kstart,kend,cnt;
   int l,n,myrank;
   double dx,dy,v1,v2,v3,charge,weightCoef;
   double ne,randTest,positionX,positionY,x,y,xPos,yPos;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New;
   DefPtcl *p;

   dx=D->dx;   dy=D->dy;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   cnt=LL->numDefPtcls;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;
   k=0;

   //position define     
   p=LL->def;
   while(p) {
     cnt=p->numDefPtcls;
     for(n=0; n<cnt; n++)
     {
       xPos=p->define[n];
       i=(int)(xPos/dx-D->minXSub+istart);
       if(i<iend && i>=istart)
       {
         yPos=p->define[n+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         if(jstart<=j && j<jend)
         {
           x=xPos/dx-D->minXSub+istart-i;
           y=yPos/dy-D->minYSub+jstart-j;
           New = (ptclList *)malloc(sizeof(ptclList));
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;

           New->x = x;           New->oldX=i+x;
           New->y = y;           New->oldY=j+y;
           New->z = 0;           New->oldZ=k +0;

           New->E1=New->E2=New->E3=0.0;
           New->B1=New->B2=New->B3=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
           New->p1=-D->gamma*D->beta+v1;
           New->p2=v2;
           New->p3=v3;
           New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
           New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
           New->weight=1.0/LL->numberInCell*weightCoef;
           New->charge=charge;
           LL->index+=1;
           New->index=LL->index;
           New->core=myrank;
         }
       }                //if(ii)
     }                  //End of for(n<cnt)
     p=p->next;
   }                    //End of while(p)
}




void loadPolygonPlasma1D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,n,istart,iend,intNum,cnt,l,t,modeX;
   double posX,v1,v2,v3,centerX,tmp,weight,charge;
   double ne,randTest,positionX,gaussCoefX,polyCoefX,weightCoef;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   centerX=LL->centerX;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   modeX=LL->modeX;
   weight=1.0/LL->numberInCell;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   srand(myrank-1);
   intNum=(int)LL->numberInCell;

   //position define      
   j=k=0;
   for(i=istart; i<iend; i++)
     {
       gsl_qrng *q1 = gsl_qrng_alloc (gsl_qrng_sobol,1);
       for(l=0; l<LL->xnodes-1; l++)
         {
           posX=(double)(i+D->minXSub-istart);
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
             ne*=tmp;
             weight=ne/LL->numberInCell*weightCoef;

             for(n=0; n<intNum; n++) {
//               positionX=randomValue(1.0);
               random1D_sobol(&positionX,q1);
               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = 0.0;
               New->oldY=j+0.0;
               New->z = 0.0;
               New->oldZ=k+0.0;

               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=v1;
               New->p2=v2;
               New->p3=v3;
               New->weight=weight;
               New->charge=charge;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            
             }		//end of for(i<intNum)

           } else;	
         } 		//end of for(lnodes)  

       gsl_qrng_free(q1);
     }			//End of for(i,j)         
}


void loadPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,t;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,tmp,weight,charge,weightCoef;
   double ne,randTest,positionX,positionY,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   k=0;
   centerX=LL->centerX; centerY=LL->centerY; centerZ=LL->centerZ;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   gaussCoefYZ=LL->gaussCoefYZ;
   polyCoefYZ=LL->polyCoefYZ;
   modeX=LL->modeX;
   modeYZ=LL->modeYZ;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       gsl_qrng *q2 = gsl_qrng_alloc (gsl_qrng_sobol,2);

       for(l=0; l<LL->xnodes-1; l++)
         for(t=0; t<LL->ynodes-1; t++)
         {
           posX=(double)(i+D->minXSub-istart);
           posY=(double)(j+D->minYSub-jstart);
           posZ=(double)(k+D->minZSub-kstart);
           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
              posY>=LL->ypoint[t] && posY<LL->ypoint[t+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             ne*=((LL->yn[t+1]-LL->yn[t])/(LL->ypoint[t+1]-LL->ypoint[t])*(posY-LL->ypoint[t])+LL->yn[t]);
//             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
//             ne*=tmp;
             tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
             ne*=tmp;
             weight=ne/LL->numberInCell*weightCoef;

             cnt=0;
             while(cnt<intNum)
             {      
//               positionX=randomValue(1.0);
//               positionY=randomValue(1.0);
               random2D_sobol(&positionX,&positionY,q2);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0.0;
               New->oldZ=0.0;

               New->E1=New->E2=New->E3=0.0;
               New->B1=New->B2=New->B3=0.0;
               v1=maxwellianVelocity(LL->temperature)/velocityC;
               v2=maxwellianVelocity(LL->temperature)/velocityC;
               v3=maxwellianVelocity(LL->temperature)/velocityC;
               New->p1=v1;      New->p2=v2;       New->p3=v3;
               New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
               New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
               New->weight=weight;
               New->charge=charge;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank;            

               cnt++; 
             }		//end of while(cnt)

           }	
         } 		//end of for(lnodes)  
       gsl_qrng_free(q2);

     }			//End of for(i,j)
         
}

void loadPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,m,n;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,weight,charge,centerX,centerY,centerZ,tmp,weightCoef;
   double ne,randTest,positionX,positionY,positionZ,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   centerX=LL->centerX;
   centerY=LL->centerY;
   centerZ=LL->centerZ;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   gaussCoefYZ=LL->gaussCoefYZ;
   polyCoefYZ=LL->polyCoefYZ;
   modeX=LL->modeX;
   modeYZ=LL->modeYZ;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define      
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
       for(k=kstart; k<kend; k++)
       {
         gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);
         for(l=0; l<LL->xnodes-1; l++)
           for(m=0; m<LL->ynodes-1; m++)
             for(n=0; n<LL->znodes-1; n++)
             {
               posX=(double)(i+D->minXSub-istart);
               posY=(double)(j+D->minYSub-jstart);
               posZ=(double)(k+D->minZSub-kstart);
 
               if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
                  posY>=LL->ypoint[m] && posY<LL->ypoint[m+1] &&
                  posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])
               {
                 ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
                 ne*=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
                 ne*=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
                 tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
                 ne*=tmp;
                 tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
                 ne*=tmp;
                 weight=ne/LL->numberInCell*weightCoef;
             
                 cnt=0;
                 while(cnt<intNum)
                 {               
//                   positionX=randomValue(1.0);
//                   positionY=randomValue(1.0);
//                   positionZ=randomValue(1.0);
                   random3D_sobol(&positionX,&positionY,&positionZ,q3);
  
                   New = (ptclList *)malloc(sizeof(ptclList)); 
                   New->next = particle[i][j][k].head[s]->pt;
                   particle[i][j][k].head[s]->pt = New;
 
                   New->x = positionX;
                   New->oldX=i+positionX;
                   New->y = positionY;
                   New->oldY=j+positionY;
                   New->z = positionZ;
                   New->oldZ=k +positionZ;
  
                   New->E1=New->E2=New->E3=0.0;
                   New->B1=New->B2=New->B3=0.0;
                   v1=maxwellianVelocity(LL->temperature)/velocityC;
                   v2=maxwellianVelocity(LL->temperature)/velocityC;
                   v3=maxwellianVelocity(LL->temperature)/velocityC;
                   New->p1=-D->gamma*D->beta+v1;
                   New->p2=v2;
                   New->p3=v3;
                   New->weight=weight;
                   New->charge=charge;
                   LL->index+=1;
                   New->index=LL->index;            
                   New->core=myrank;            
   
                   cnt++; 
                 }		//end of while(cnt)
               } else ;		//End of if (l,m,n)
             }			//End of for(l,m,n)

         gsl_qrng_free(q3);
       }		//End of for(i,j,k)         
}

double maxwellianVelocity(double temperature)
{
   double vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);
   
   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}

double gaussian_dist(double sig)
{
   double r,prob,v,random;
   int intRand,randRange=1e5;

   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 4.0*(random-0.5);
      prob=exp(-v*v);
   }
   return sig*v;
}
