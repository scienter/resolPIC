#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>

double randomValue(double beta);
double maxwellianVelocity(double temperature);
void loadMovingPlasma_crystal(Domain *D,LoadList *LL,int s);
double applyFunctionX(int mode,double centerX,double x,double gaussCoefX,double polyCoefX);
double applyFunctionYZ(int mode,double centerY,double y,double centerZ,double z,double gaussCoefYZ,double polyCoefYZ);
void random1D_sobol(double *x,gsl_qrng *q);
void random2D_sobol(double *x,double *y,gsl_qrng *q);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q);
void loadMovingDefinedPlasma2D(Domain *D,LoadList *LL,int s);
void loadMovingDefinedPlasma3D(Domain *D,LoadList *LL,int s);
void loadMovingPolygonPlasma1D(Domain *D,LoadList *LL,int s,int iteration);
void loadMovingPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration);
void loadMovingPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration);


void loadMovingPlasma(Domain *D,LoadList *LL,int s,int iteration)
{
  int rankX;
  int myrank, nTasks;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  rankX=myrank/(D->M*D->N);
  switch((LL->type-1)*3+D->dimension)  {
  //1D
  case 1:
    if(rankX==D->L-1)
      loadMovingPolygonPlasma1D(D,LL,s,iteration); 
    else	;
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //2D
  case (Polygon-1)*3+2:
    if(rankX==D->L-1)
      loadMovingPolygonPlasma2D(D,LL,s,iteration); 
    else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;
  case (Defined-1)*3+2:
    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime) {
      if(rankX==D->L-1) loadMovingDefinedPlasma2D(D,LL,s); 
      else ;
    }    else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  //3D
  case (Polygon-1)*3+3:
    if(rankX==D->L-1)
      loadMovingPolygonPlasma3D(D,LL,s,iteration); 
    else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;
  case (Defined-1)*3+3:
    if(LL->minLoadTime<=iteration && iteration<=LL->maxLoadTime) {
      if(rankX==D->L-1) loadMovingDefinedPlasma3D(D,LL,s); 
      else ;
    } else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  default:
    ;
  }
}

void loadMovingDefinedPlasma2D(Domain *D,LoadList *LL,int s)
{
   int ii,i,j,k,istart,iend,jstart,jend,kstart,kend;
   int n,myrank,cnt;
   double dx,dy,v1,v2,v3,weightCoef;
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
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;
   k=0;

   //position define      
   i=iend-1;
   p=LL->def;
   while(p) {
     cnt=p->numDefPtcls;
     for(n=0; n<cnt; n++) {
       xPos=p->define[n];
       ii=(int)(xPos/dx-D->minXSub+istart);
       if(i==ii)
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
           New->charge=LL->charge;
           LL->index+=1;
           New->index=LL->index;           
           New->core=myrank; 
         }
       }		//if(ii)
     }			//End of for(n<cnt)
     p=p->next;
   }			//End of while(p)
}

void loadMovingDefinedPlasma3D(Domain *D,LoadList *LL,int s)
{
   int ii,i,j,k,istart,iend,jstart,jend,kstart,kend;
   int n,myrank,cnt;
   double dx,dy,dz,v1,v2,v3,weightCoef;
   double ne,randTest,positionX,positionY,x,y,z,xPos,yPos,zPos;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New;   
   DefPtcl *p;

   dx=D->dx; dy=D->dy; dz=D->dz;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   cnt=LL->numDefPtcls;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;


   //position define      
   i=iend-1;
   p=LL->def;
   while(p) {
     cnt=p->numDefPtcls;
     for(n=0; n<cnt; n++)
     {
       xPos=p->define[n];
       ii=(int)(xPos/dx-D->minXSub+istart);
       if(i==ii)
       {
         yPos=p->define[n+cnt];
         j=(int)(yPos/dy-D->minYSub+jstart);
         zPos=p->define[n+2*cnt];
         k=(int)(zPos/dz-D->minZSub+kstart);
         if(jstart<=j && j<jend && kstart<=k && k<kend)
         {
           x=xPos/dx-D->minXSub+istart-i;
           y=yPos/dy-D->minYSub+jstart-j;
           z=zPos/dz-D->minZSub+kstart-k;
           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;
  
           New->x = x;           New->oldX=i+x;
           New->y = y;           New->oldY=j+y;
           New->z = z;           New->oldZ=k +z;
  
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
           New->charge=LL->charge;
           LL->index+=1;
           New->index=LL->index;            
           New->core=myrank; 
         }
       }		//if(ii)
     }			//End of for(n<cnt)
     p=p->next;
   }			//End of while(p)
}

void loadMovingPolygonPlasma1D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,intNum,cnt,l,myrank,modeX;
   double posX,v1,v2,v3,centerX,tmp,weight,charge,weightCoef;
   double ne,randTest,positionX,gaussCoefX,polyCoefX;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;
   iend=D->iend;
   j=k=0;
   centerX=LL->centerX;
   gaussCoefX=LL->gaussCoefX;
   polyCoefX=LL->polyCoefX;
   modeX=LL->modeX;
   charge=LL->charge;
   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   srand(iteration*(myrank+1));
   intNum=(int)LL->numberInCell;

   //position define      
   i=iend-1;
       gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,1);

       for(l=0; l<LL->xnodes-1; l++)
         {
           posX=(double)(i+D->minXSub-istart);

           if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])
           {
             ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
             tmp=applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
             ne*=tmp;

             weight=ne/LL->numberInCell*weightCoef;
             
             cnt=0;
             while(cnt<intNum)
             {               
               positionX=randomValue(1.0);
//               random1D_sobol(&positionX,q);
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
               New->p1=-D->gamma*D->beta+v1;
               New->p2=v2;
               New->p3=v3;
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
         } 		//end of for(lnodes)  i
       gsl_qrng_free(q);
         
}

void loadMovingPolygonPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,t,myrank;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,tmp,weight,charge,weightCoef;
   double ne,randTest,positionX,positionY,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   k=0;
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

   srand(iteration*(myrank+1));
   intNum=(int)LL->numberInCell;

   //position define      
   i=iend-1;
     for(j=jstart; j<jend; j++)
     {
       gsl_qrng *q = gsl_qrng_alloc (gsl_qrng_sobol,2);

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
               random2D_sobol(&positionX,&positionY,q);

               New = (ptclList *)malloc(sizeof(ptclList)); 
               New->next = particle[i][j][k].head[s]->pt;
               particle[i][j][k].head[s]->pt = New;
 
               New->x = positionX;
               New->oldX=i+positionX;
               New->y = positionY;
               New->oldY=j+positionY;
               New->z = 0;
               New->oldZ=k +0;

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
               New->weight=weight;
               New->charge=charge;
               LL->index+=1;
               New->index=LL->index;            
               New->core=myrank; 

               cnt++; 
             }		//end of while(cnt)
           }	
         } 		//end of for(lnodes)  i
       gsl_qrng_free(q);

     }			//End of for(i,j)
         
}

void loadMovingPolygonPlasma3D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,m,n,myrank;
   int modeX,modeYZ;
   double posX,posY,posZ,v1,v2,v3,centerX,centerY,centerZ,weight,charge,tmp,weightCoef;
   double ne,randTest,positionX,positionY,positionZ,gaussCoefX,polyCoefX,gaussCoefYZ,polyCoefYZ;
   Particle ***particle;
   particle=D->particle;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   ptclList *New,*p;   

   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   kstart=D->kstart; kend=D->kend;
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

   srand(iteration*(myrank+1));
   intNum=(int)LL->numberInCell;

   //position define      
   i=iend-1;
   for(j=jstart; j<jend; j++)
     for(k=kstart; k<kend; k++)
     {
       gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);

       for(l=0; l<LL->xnodes-1; l++)
         for(m=0; m<LL->ynodes-1; m++)
           for(n=0; n<LL->znodes-1; n++)
           {
             posX=i+D->minXSub-istart;
             posY=j+D->minYSub-jstart;
             posZ=k+D->minZSub-kstart;
 
             if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1] &&
                posY>=LL->ypoint[m] && posY<LL->ypoint[m+1] &&
                posZ>=LL->zpoint[n] && posZ<LL->zpoint[n+1])
             {
               ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
               ne*=((LL->yn[m+1]-LL->yn[m])/(LL->ypoint[m+1]-LL->ypoint[m])*(posY-LL->ypoint[m])+LL->yn[m]);
               ne*=((LL->zn[n+1]-LL->zn[n])/(LL->zpoint[n+1]-LL->zpoint[n])*(posZ-LL->zpoint[n])+LL->zn[n]);
//               applyFunctionX(modeX,centerX,posX,gaussCoefX,polyCoefX);
               tmp=applyFunctionYZ(modeYZ,centerY,posY,centerZ,posZ,gaussCoefYZ,polyCoefYZ);
               ne*=tmp;

               weight=ne/LL->numberInCell*weightCoef;	//it is the double number of superparticles.
             
               cnt=0;
               while(cnt<intNum)
               {               
//                 positionX=randomValue(1.0);
//                 positionY=randomValue(1.0);
//                 positionZ=randomValue(1.0);
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
                 New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
                 New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
                 New->weight=weight;
                 New->charge=charge;
                 LL->index+=1;
                 New->index=LL->index;            
                 New->core=myrank; 

                 cnt++; 
               }		//end of while(cnt)

             }			//End of if(pos...)
           } 			//end of for(l,m,n)  

       gsl_qrng_free(q3);
     }				//End of for(i,j,k)
         
}

void assignDefParticle(Domain *D)
{
  int n,cnt,flag;
  double tmp,tmpX,tmpY,tmpZ,tmpR,factY,factZ,randX,randY,randZ;
  LoadList *LL,*prevL;
  DefPtcl *p,*prevP;

  factY=factZ=0.0;
  if(D->dimension>1) factY=1.0; else;
  if(D->dimension>2) factZ=1.0; else;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      if(LL->pair==OFF) {
        p=LL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=LL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);
              for(n=0; n<cnt; n++)  {
                flag=0;
                while(flag==0) {
//                  random2D_sobol(&randX,&randY,q3); randZ=0.0;
                  random3D_sobol(&randX,&randY,&randZ,q3);
                  tmpX=-LL->RDef+randX*2.0*LL->RDef;
                  tmpY=(-LL->RDef+randY*2.0*LL->RDef)*factY;
                  tmpZ=(-LL->RDef+randZ*2.0*LL->RDef)*factZ;
//                  tmp=(double)(randomValue(1.0));
//                  tmpX=-LL->RDef+tmp*2.0*LL->RDef;
//                  tmp=(double)(randomValue(1.0));
//                  tmpY=(-LL->RDef+tmp*2.0*LL->RDef)*factY;
//                  tmp=(double)(randomValue(1.0));
//                  tmpZ=(-LL->RDef+tmp*2.0*LL->RDef)*factZ;
                  tmpR=sqrt(tmpX*tmpX+tmpY*tmpY+tmpZ*tmpZ);
                  if(tmpR<=LL->RDef) {
                    flag=1;
                    tmpX=p->xPos+tmpX;
                    tmpY=p->yPos+tmpY;
                    tmpZ=p->zPos+tmpZ;
                  } else  flag=0;
                }
                p->define[n]=tmpX;
                p->define[n+cnt]=tmpY;
                p->define[n+2*cnt]=tmpZ;
              }	//End of for(n<cnt)
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
        }
      } else { 	//if pair==ON
        p=LL->def; prevP=prevL->def;
        while(p) {
          if(p->flag==OFF) {
            if(p->xPos<=LL->createGuard) {
              cnt=p->numDefPtcls=prevL->numDefPtcls;           
              p->define=(double *)malloc((cnt*3)*sizeof(double ));
              for(n=0; n<cnt; n++)  {
                p->define[n]=prevP->define[n];
                p->define[n+cnt]=prevP->define[n+cnt];
                p->define[n+2*cnt]=prevP->define[n+2*cnt];
              }
              p->flag=ON;  
            } else ;
          } else ;	//End of if(flag==Off)
          p=p->next;
          prevP=prevP->next;	
        } 
      }
    } else ; 		//End of if(LL->type==Defined) 
    prevL=LL;
    LL=LL->next;
  }
}

void deleteDefParticle(Domain *D)
{
  int cnt;
  LoadList *LL;
  DefPtcl *p,*prev;

  LL=D->loadList;
  while(LL->next)      {
    if(LL->type==Defined) {
      p=LL->def;
      cnt=1;
      while(p) {
        if(cnt==1) prev=p; else ;
  
        if(p->xPos<=LL->deleteGuard) {
          if(cnt==1)  {
            LL->def=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=LL->def;
            cnt=1;
          } else {
            prev->next=p->next;
            p->next=NULL;
            free(p->define);
            free(p);
            p=prev->next;
          } 
        }
        else  {
          prev=p;
          p=p->next;
          cnt++;
        }
      }
    }	else ;	//End of if(LL->type==Defined)
    LL=LL->next;
  }
}

void saveDefParticle(Domain *D,int iteration)
{
  int n,cnt,num;
  double x,y,z;
  char name[100];
  FILE *out;

  LoadList *LL;
  DefPtcl *p;
  LL=D->loadList;
  num=0;
  while(LL->next)      {
    sprintf(name,"def%d_%d",iteration,num);
    out = fopen(name,"w");    

    p=LL->def;
    while(p) {
      cnt=p->numDefPtcls;
      for(n=0; n<cnt; n++)  {
        x=p->define[n]*D->lambda;
        y=p->define[n+cnt]*D->lambda;
        z=p->define[n+2*cnt]*D->lambda;
        fprintf(out,"%g %g %g\n",x,y,z);
      }	//End of for(n<cnt)
      p=p->next;
    }
    fclose(out);
    num++;
    LL=LL->next;
  }
}
