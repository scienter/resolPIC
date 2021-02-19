#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_qrng.h>


double maxwellianVelocity(double temperature);
void random1D_sobol(double *x,gsl_qrng *q1);
void random2D_sobol(double *x,double *y,gsl_qrng *q2);
void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q3);
void loadBeamPlasma1D(Domain *D,LoadList *LL,int s,int iteration);
void loadBeamPlasma2D(Domain *D,LoadList *LL,int s,int iteration);
double gaussian_dist(double sig);
void MPI_TransferRho_Xplus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Xminus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Yplus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_TransferRho_Yminus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_Transfer1F_Xplus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_Transfer1F_Xminus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_Transfer1F_Yplus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_Transfer1F_Yminus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_Transfer3F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void MPI_Transfer3F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void MPI_Transfer3F_Yplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share);
void MPI_Transfer3F_Yminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share);
void MPI_TransferJ_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void MPI_TransferJ_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void sor2D(Domain *D,double ***Phi,double ***f,double a,double b,double c,double d,double e,double rjac,double ratio,int nn);
void absorb_UD(Domain *D,double *upr,double *dnr,double *upd,double *dnd,double y,double upL,double downL,double LdU,double LdD,double rr,double rd);


void loadBeam(Domain *D,LoadList *LL,int s,int iteration)
{
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  switch((LL->type-1)*3+D->dimension)  {
  //1D
  case ((Polygon-1)*3+1):
	  ;
    break;

  //2D
  case (Polygon-1)*3+2:
          ;
    break;

  case (Defined-1)*3+2:
    ;
    break;

  //3D
  case (Polygon-1)*3+3:
    ;
    break;
  case (Defined-1)*3+3:
    ;
    break;

  case (Beam-1)*3+1:
    if(iteration==LL->loadingStep) {
      loadBeamPlasma1D(D,LL,s,iteration);
    }    else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;

  case (Beam-1)*3+2:
    if(iteration==LL->loadingStep) {
      loadBeamPlasma2D(D,LL,s,iteration);
    }    else ;
    MPI_Barrier(MPI_COMM_WORLD);
    break;


  default:
    ;
  }
}


void loadBeamPlasma1D(Domain *D,LoadList *LL,int s,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,jC,kC;
   int minYSub,maxYSub,totalCnt;
   double posX,posY,posZ,tmp,weight,charge,weightCoef,ne,nenergy,testY;
   double px,py,pz,positionX,positionY,y,yPrime,z,zPrime,sigYPrime,sigZPrime;
   double dt,dx,dy,dz,emitY,emitZ,betaY,betaZ,sigY,sigZ,gamma0,dGam,gamma,sigGam,delGam;
   double distance,vx,delT,sqrt2;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minYSub=D->minYSub;
   maxYSub=D->maxYSub;
   j=k=0;

   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;
   jC=0;
   kC=0;

   sqrt2=1.414214;
   gamma0=LL->energy/mc2;
   dGam=LL->spread*gamma0;
   emitY=LL->emitY/gamma0;  
   emitZ=LL->emitZ/gamma0;
   sigY=LL->sigY/D->lambda;	//normalized
   sigZ=LL->sigZ/D->lambda;	//normalized
   sigYPrime=emitY/sigY*0.25/D->lambda;
   sigZPrime=emitZ/sigZ*0.25/D->lambda;

   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   charge=LL->charge;
   distance=LL->focus/D->lambda-LL->xpoint[LL->xnodes-1]*D->dx;	//normalized
   vx=sqrt(gamma0*gamma0-1.0)/gamma0;				//normalized
   delT=distance/vx;						//normalized

   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define      
   for(i=istart; i<iend; i++)
   {
     gsl_qrng *q2 = gsl_qrng_alloc (gsl_qrng_sobol,2);
     for(l=0; l<LL->xnodes-1; l++) 
     {
       posX=(double)(i+D->minXSub-istart);
       if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])
       {
         ne=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l]);
         nenergy=((LL->xenergy[l+1]-LL->xenergy[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xenergy[l]);
	 totalCnt=intNum;
         weight=ne/LL->numberInCell*weightCoef;

         cnt=0;
         while(cnt<totalCnt)  {      
           positionX=randomValue(1.0);
           y=gaussian_dist(sigY*sqrt2);
           yPrime=gaussian_dist(sigYPrime*sqrt2);
           z=gaussian_dist(sigZ*sqrt2);
           zPrime=gaussian_dist(sigZPrime*sqrt2);
	   sigGam=dGam*0.5*1.201122;		// 1.201122=1/sqrt(ln(2))
	   delGam=gaussian_dist(sigGam);
           gamma=gamma0*nenergy+delGam;

	   px=sqrt((gamma*gamma-1.0)/(1.0+yPrime*yPrime+zPrime*zPrime));
	   py=yPrime*px;
	   pz=zPrime*px;
       
	   // transverse postion at start corresponding to the focus
	   y-=delT*py/gamma;
	   z-=delT*pz/gamma;

           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;
 
           New->x = positionX; New->y = 0.0; New->z = 0.0;
           New->oldX=i+positionX-px/gamma*dt/dx;
           New->oldY=0.0; New->oldZ=0.0;
  
           New->E1=New->E2=New->E3=0.0;
           New->B1=New->B2=New->B3=0.0;
           New->p1=px;      New->p2=py;       New->p3=pz;
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

void loadBeamPlasma2D(Domain *D,LoadList *LL,int s,int iteration)
{
   int ii,jj,i,j,k,istart,iend,jstart,jend,kstart,kend,intNum,cnt,l,jC,kC;
   int minYSub,maxYSub,totalCnt,i1,j1,intXc,intYc,oldI,oldJ;
   double posX,posY,posZ,tmp,weight,charge,weightCoef,ne,nenergy,testY,wTimesQ;
   double px,py,pz,positionX,positionY,yPrime,zPrime,sigYPrime,sigZPrime;
   double dt,dx,dy,dz,emitY,emitZ,sigY,sigZ,gamma0,dGam,gamma,sigGam,delGam;
   double distance,delT,sqrt2,wt,fVal,coeff,wx[2],wy[2],WX[2],WY[2];
   double rjac,ratio,x,y,z,a,b,c,d,e,vx,vy,vz,xc,yc,xcc,ycc,beta2,beta,oldX,oldY,x1,y1;
   FILE *out;
   char name[100];

   MPI_Status status;
   int myrank,nTasks,rankX,rank;
   double recv;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   Particle ***particle;
   particle=D->particle;

   ptclList *New,*p;   

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   kstart=D->kstart;   kend=D->kend;
   minYSub=D->minYSub;
   maxYSub=D->maxYSub;
   k=0;

   dx=D->dx; dy=D->dy; dz=D->dz; dt=D->dt;
   jC=0;
   kC=0;

   sqrt2=1.414214;
   gamma0=LL->energy/mc2+1.0;
   dGam=LL->spread*gamma0;
   emitY=LL->emitY/gamma0;  
   emitZ=LL->emitZ/gamma0;
   sigY=LL->sigY/D->lambda;	//normalized
   sigZ=LL->sigZ/D->lambda;	//normalized
   sigYPrime=emitY/sigY*0.25/D->lambda;
   sigZPrime=emitZ/sigZ*0.25/D->lambda;

   if(LL->species==Test) weightCoef=0.0; else weightCoef=1.0;

   charge=LL->charge;
   distance=LL->focus/D->lambda-LL->xpoint[LL->xnodes-1]*D->dx;	//normalized
   vx=sqrt(gamma0*gamma0-1.0)/gamma0;				//normalized
   if(vx==0.0) delT=0.0;
   else        delT=distance/vx;						//normalized
   coeff=LL->density/LL->criticalDensity;

   srand(myrank+1);
   intNum=(int)LL->numberInCell;

   //position define      
   for(i=istart; i<iend; i++)
   {
//     srand(i-istart+D->minXSub);
     gsl_qrng *q2 = gsl_qrng_alloc (gsl_qrng_sobol,2);
     for(l=0; l<LL->xnodes-1; l++) {
       posX=(double)(i+D->minXSub-istart);
       if(posX>=LL->xpoint[l] && posX<LL->xpoint[l+1])
       {
         ne=(LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xn[l];
         nenergy=(LL->xenergy[l+1]-LL->xenergy[l])/(LL->xpoint[l+1]-LL->xpoint[l])*(posX-LL->xpoint[l])+LL->xenergy[l];
	 totalCnt=ne*sigY/dy*intNum*sqrt(2*M_PI);
         weight=1.0/LL->numberInCell*weightCoef;

         cnt=0;
         while(cnt<totalCnt)  {      
           positionX=randomValue(1.0);
           y=gaussian_dist(sigY*sqrt2);
           yPrime=gaussian_dist(sigYPrime*sqrt2);
           z=gaussian_dist(sigZ*sqrt2);
           zPrime=gaussian_dist(sigZPrime*sqrt2);
	   sigGam=dGam*0.5*1.201122;		// 1.201122=1/sqrt(ln(2))
	   delGam=gaussian_dist(sigGam);
           gamma=gamma0*nenergy+delGam;

	   px=sqrt((gamma*gamma-1.0)/(1.0+yPrime*yPrime+zPrime*zPrime));
	   py=yPrime*px;
	   pz=zPrime*px;
       
	   // transverse postion at start corresponding to the focus
//	   printf("delT=%g, distance=%g, vx=%g,focus=%g,lambda=%g,pointX=%g,dx=%g\n",delT,distance,vx,LL->focus,D->lambda,LL->xpoint[LL->xnodes-1],D->dx);
	   y-=delT*py/gamma;
	   z-=delT*pz/gamma;

	   testY=y/dy;
           if(testY>=minYSub && testY<maxYSub) {
             j=(int)(y/dy-minYSub+jstart);
	     positionY=y/dy-minYSub-(int)(y/dy-minYSub);

	     New = (ptclList *)malloc(sizeof(ptclList)); 
             New->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = New;
 
             New->x = positionX; New->y = positionY; New->z = 0.0;
             New->oldX=i+positionX-px/gamma*dt/dx;
             New->oldY=j+positionY-py/gamma*dt/dy;
             New->oldZ=0.0;
  
             New->E1=New->E2=New->E3=0.0;
             New->B1=New->B2=New->B3=0.0;
             New->p1=px;      New->p2=py;       New->p3=pz;
             New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
             New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
             New->weight=weight;
             New->charge=charge;
             LL->index+=1;
             New->index=LL->index;            
             New->core=myrank;            

             wx[0]=1.0-positionX; wx[1]=positionX;
             wy[0]=1.0-positionY; wy[1]=positionY;
	     for(ii=0; ii<2; ii++)
	       for(jj=0; jj<2; jj++)
	         D->Den[i+ii][j+jj][k]-=coeff*charge*weight*wx[ii]*wy[jj]*4.0*M_PI*M_PI;
/*
	     wTimesQ=charge*weight;
             x1=(positionX+0.5)-(int)(positionX+0.5);
             y1=(positionY+0.5)-(int)(positionY+0.5);
             i1=i+(int)(positionX+0.5)-1;
             j1=j+(int)(positionY+0.5)-1;

             WX[1]=x1; WX[0]=1.0-x1;
             WY[1]=y1; WY[0]=1.0-y1;
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 D->Jx[i1+ii][j+jj][k]+=WX[ii]*wy[jj]*px/gamma*coeff*wTimesQ;		
		 D->Jy[i+ii][j1+jj][k]+=wx[ii]*WY[jj]*py/gamma*coeff*wTimesQ;
                 D->Jz[i+ii][j+jj][k]+=wx[ii]*wy[jj]*pz/gamma*coeff*wTimesQ;
               }
*/
/*	     
	     oldX=i+positionX+0.5*px/gamma*dt/dx;
             oldY=j+positionY+0.5*py/gamma*dt/dy;
             oldI=(int)oldX;   oldJ=(int)oldY;
	     wx[1]=oldX-oldI; wx[0]=1.0-wx[0];
	     wy[1]=oldY-oldJ; wy[0]=1.0-wy[0];
	     for(ii=0; ii<2; ii++)
	       for(jj=0; jj<2; jj++)
	         D->DenOld[oldI+ii][oldJ+jj][k]-=coeff*charge*wx[ii]*wy[jj]*weight*4.0*M_PI*M_PI;
*/
		 
	     xc=i+positionX+0.5*px/gamma*dt/dx;
             yc=j+positionY+0.5*py/gamma*dt/dy;
             intXc=(int)xc;   intYc=(int)yc;
             xcc=xc-intXc;    ycc=yc-intYc;
	     wx[1]=xcc;		           wx[0]=1.0-wx[1];
	     wy[1]=ycc;		           wy[0]=1.0-wy[1];
	     WX[1]=xcc+0.5-(int)(xcc+0.5); WX[0]=1.0-WX[1];
	     WY[1]=ycc+0.5-(int)(ycc+0.5); WY[0]=1.0-WY[1];
	     i1=(int)(xc+0.5)-1;       
	     j1=(int)(yc+0.5)-1;

	     vx=px/gamma;
	     vy=py/gamma;
	     vz=pz/gamma;
//	     for(ii=0; ii<2; ii++)
//	       for(jj=0; jj<2; jj++) {
//	         D->CurX[i1+ii][intYc+jj][k]-=coeff*charge*WX[ii]*wy[jj]*2.0*M_PI*vx;
//	         D->CurY[intXc+ii][j1+jj][k]-=coeff*charge*wx[ii]*WY[jj]*2.0*M_PI*vy;
//	         D->CurZ[intXc+ii][intYc+jj][k]-=coeff*charge*wx[ii]*wy[jj]*2.0*M_PI*vz;
//	         D->Jx[i1+ii][intYc+jj][k]+=coeff*charge*weight*WX[ii]*wy[jj]*vx;
//                 D->Jy[intXc+ii][j1+jj][k]+=coeff*charge*weight*wx[ii]*WY[jj]*vy;
//	         D->Jz[intXc+ii][intYc+jj][k]+=coeff*charge*weight*wx[ii]*wy[jj]*vz;
//	       }

	   } else ;
           cnt++; 	   
         }		//end of while(cnt)
       }	
     } 		//end of for(lnodes)  
     gsl_qrng_free(q2);
   }			//End of for(i,j)

   if(D->L>1)  {
     MPI_TransferRho_Xplus(D,D->Den,D->nySub+5,1,3);
     MPI_TransferRho_Xminus(D,D->Den,D->nySub+5,1,3);
   }  else     ;
   if(D->M>1)  {
     MPI_TransferRho_Yplus(D,D->Den,D->nxSub+5,1,3);
     MPI_TransferRho_Yminus(D,D->Den,D->nxSub+5,1,3);
   }  else     ;

   a=1.0/dx/dx/gamma0/gamma0;
   b=a;
   c=1.0/dy/dy;
   d=c;
   e=-2.0*(a+c);
   rjac=LL->rjac;
   ratio=LL->ratio;

   sor2D(D,D->Phi,D->Den,a,b,c,d,e,rjac,ratio,1);   
//   sor2D(D,D->PhiOld,D->DenOld,a,b,c,d,e,rjac,ratio,2);   
//   sor2D(D,D->Ax,D->CurX,a,b,c,d,e,rjac,ratio,2);   
//   sor2D(D,D->Ay,D->CurY,a,b,c,d,e,rjac,ratio,3);   
//   sor2D(D,D->Az,D->CurZ,a,b,c,d,e,rjac,ratio,4);   


   //field assign
   double upL,downL,LdU,LdD,rr,rd,tmpr,upr,upd,dnr,dnd;
   beta2=1.0-1.0/gamma0/gamma0;
   beta=sqrt(beta2);
   upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
   downL=(double)(D->minYDomain+D->pmlCellDown);
   LdU=D->pmlCellUp;
   LdD=D->pmlCellDown;
   rr=D->pmlr;
   rd=D->pmld;
   upr=upd=dnr=dnd=1.0;
   if(D->fieldType==Yee || D->fieldType==Pukhov) {
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++) {
         y=(j-jstart)+minYSub;
//       absorb_UD(D,&upr,&dnr,&upd,&dnd,y,upL,downL,LdU,LdD,rr,rd);
//       tmpr=upr*dnr;
         D->Ex[i][j][k]+=-0.5/M_PI/gamma0/gamma0*(D->Phi[i+1][j][k]-D->Phi[i][j][k])/dx;
         D->Ey[i][j][k]+=-0.5/M_PI*(D->Phi[i][j+1][k]-D->Phi[i][j][k])/dy;
         D->Bz[i][j][k]+=-0.25/M_PI*beta*(D->Phi[i][j+1][k]+D->Phi[i-1][j+1][k]-D->Phi[i][j][k]-D->Phi[i-1][j][k])/dy;
//       D->Bz[i][j][k]=-0.5/M_PI*beta*(D->Phi[i][j+1][k]-D->Phi[i][j][k])/dy;
//       D->Bz[i][j][k]=-0.25/M_PI*beta*(D->PhiOld[i+1][j+1][k]+D->PhiOld[i][j+1][k]-D->PhiOld[i+1][j][k]-D->PhiOld[i][j][k])/dy;
       }
     if(D->M>1)  {
       MPI_Transfer3F_Yminus(D,D->Ex,D->Ey,D->Bz,D->nxSub,1,3);
       MPI_Transfer3F_Yplus(D,D->Ex,D->Ey,D->Bz,D->nxSub,1,3);
     } else      ;
     if(D->L>1)  {
       MPI_Transfer3F_Xminus(D,D->Ex,D->Ey,D->Bz,D->nySub+5,1,3);
       MPI_Transfer3F_Xplus(D,D->Ex,D->Ey,D->Bz,D->nySub+5,1,3);
     } else      ;
   } else if(D->fieldType==Split) {
     for(i=istart; i<iend; i++)
       for(j=jstart; j<jend; j++) {
         y=(j-jstart)+minYSub;
         D->Ex[i][j][k]+=-0.25/M_PI/gamma0/gamma0*(D->Phi[i+1][j][k]-D->Phi[i-1][j][k])/dx;
         D->PrC[i][j][k]+=-0.5/M_PI*(D->Phi[i][j+1][k]-D->Phi[i][j][k])/dy;
         D->Pr[i][j][k]+=-0.25/M_PI*(D->Phi[i+1][j+1][k]+D->Phi[i][j+1][k]-D->Phi[i+1][j][k]-D->Phi[i][j][k])/dy;
       }
     if(D->M>1)  {
       MPI_Transfer3F_Yminus(D,D->Ex,D->Pr,D->PrC,D->nxSub,1,3);
       MPI_Transfer3F_Yplus(D,D->Ex,D->Pr,D->PrC,D->nxSub,1,3);
     } else      ;
     if(D->L>1)  {
       MPI_Transfer3F_Xminus(D,D->Ex,D->Pr,D->PrC,D->nySub+5,1,3);
       MPI_Transfer3F_Xplus(D,D->Ex,D->Pr,D->PrC,D->nySub+5,1,3);
     } else      ;
   } else ;

/*
   int rankY;
   rankY=myrank%D->M;
   if(rankY==0) {
     for(i=istart; i<iend; i++) {
       j=jstart-1;
       y=(j-jstart)+minYSub;
       D->Bz[i][j][k]=-0.25/M_PI*beta*(D->Phi[i][j+1][k]+D->Phi[i-1][j+1][k]-D->Phi[i][j][k]-D->Phi[i-1][j][k])/dy;
     }
   } else ;
*/

/*     
   sprintf(name,"poisson_%d",myrank);
   out = fopen(name,"w");
   for(i=istart; i<iend; i++) {
     x=(i-istart+D->minXSub)*dx*D->lambda;
     for(j=jstart-1; j<jend+1; j++) {
       y=(j-jstart+D->minYSub)*dy*D->lambda;
       fprintf(out,"%g %g %g %g %g %g %g\n",x,y,D->Ex[i][j][k],D->Ey[i][j][k],D->Bz[i][j][k],D->Phi[i][j][k],D->Den[i][j][k]);
     }
     fprintf(out,"\n");
   }
   fclose(out);
   printf("%s is made\n",name);
*/  
}

void sor2D(Domain *D,double ***Phi,double ***f,double a,double b,double c,double d,double e,double rjac,double ratio,int nn)
{
   int i,j,k=0,jsw,maxIter=10000,istart,iend,jstart,jend,cnt,minY;
   double x,y,refMin,min,error,Err,w=1.0;
   char name[100];
   FILE *out;

   MPI_Status status;
   int myrank,nTasks,rankX,rankY,rank,offX1,offX2,offY1,offY2;
   double recv;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;

   rankX=myrank/D->M;
   if(D->L>1) {
     if(rankX==0) { offX1=1; offX2=0; }
     else if(rankX==D->L-1) { offX1=0; offX2=1; }
     else { offX1=0; offX2=0; }
   } else { offX1=1; offX2=1; }

   rankY=myrank%D->M;
   if(D->M>1) {
     if(rankY==0) { offY1=1; offY2=0; }
     else if(rankY==D->M-1) { offY1=0; offY2=0; }
     else { offY1=0; offY2=0; }
   } else { offY1=1; offY2=0; }
//   if(rankY>0 && rankY<D->M-1) jstart-=1; else ;
//lala
/*
   sprintf(name,"rho_%d",myrank);
   out = fopen(name,"w");
   for(i=istart+off1; i<iend-off2; i++) {
     x=(i-istart+D->minXSub)*D->dx*D->lambda;
     for(j=jstart+1; j<jend-1; j++) {
       y=(j-jstart+D->minYSub)*D->dy*D->lambda;
       fprintf(out,"%g %g %g\n",x,y,D->Den[i][j][k]);
     }
     fprintf(out,"\n");
   }
   fclose(out);
   printf("%s is made\n",name);
*/
   // resid
   refMin=0.0;
   recv=0.0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++) {
       Err=-f[i][j][k];
       refMin+=fabs(Err);
     }
   if(myrank!=0)  MPI_Send(&refMin,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
   else {
     for(rank=1; rank<nTasks; rank++) {
       MPI_Recv(&recv,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
       refMin+=recv;
     }
   }
   MPI_Bcast(&refMin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   

//   rjac=0.9998;
//   ratio=1e-5;
   error=ratio*refMin;

   minY=D->minYSub-D->minYDomain;
   min=refMin;
   cnt=1;
   while(min>error && cnt<maxIter) {
     min=0.0;

     //first step
     if(D->minXSub%2==0 && minY%2==0) jsw=1;
     else if(D->minXSub%2==0 && minY%2==1) jsw=0;
     else if(D->minXSub%2==1 && minY%2==0) jsw=0;
     else jsw=1;
     for(i=istart+offX1; i<iend-offX2; i++) {
       for(j=jstart+jsw+offY1; j<jend-offY2; j+=2) {
         Err=a*Phi[i+1][j][k]+b*Phi[i-1][j][k]+c*Phi[i][j+1][k]+d*Phi[i][j-1][k]+e*Phi[i][j][k]-f[i][j][k];
	 min+=fabs(Err);
	 Phi[i][j][k]-=w*Err/e;
       }
       jsw=1-jsw;
     }
     w=(cnt==1 ? 1.0/(1.0-0.5*rjac*rjac) : 1.0/(1.0-0.25*rjac*rjac*w));

     if(D->L>1)  {
       MPI_Transfer1F_Xminus(D,Phi,D->nySub+5,1,3);
       MPI_Transfer1F_Xplus(D,Phi,D->nySub+5,1,3);
     }  else     ;
     if(D->M>1)  {
       MPI_Transfer1F_Yminus(D,Phi,D->nxSub+5,1,3);
       MPI_Transfer1F_Yplus(D,Phi,D->nxSub+5,1,3);
     }  else     ;


     //second step
     if(D->minXSub%2==0 && minY%2==0) jsw=0;
     else if(D->minXSub%2==0 && minY%2==1) jsw=1;
     else if(D->minXSub%2==1 && minY%2==0) jsw=1;
     else jsw=0;
     for(i=istart+offX1; i<iend-offX2; i++) {
       for(j=jstart+jsw+offY1; j<jend-offY2; j+=2) {
         Err=a*Phi[i+1][j][k]+b*Phi[i-1][j][k]+c*Phi[i][j+1][k]+d*Phi[i][j-1][k]+e*Phi[i][j][k]-f[i][j][k];
	 min+=fabs(Err);
	 Phi[i][j][k]-=w*Err/e;
       }
       jsw=1-jsw;
     }
     w=1.0/(1.0-0.25*rjac*rjac*w);
     if(D->L>1)  {
       MPI_Transfer1F_Xminus(D,Phi,D->nySub+5,1,3);
       MPI_Transfer1F_Xplus(D,Phi,D->nySub+5,1,3);
     }  else     ;
     if(D->M>1)  {
       MPI_Transfer1F_Yminus(D,Phi,D->nxSub+5,1,3);
       MPI_Transfer1F_Yplus(D,Phi,D->nxSub+5,1,3);
     }  else     ;

     recv=0.0;
     if(myrank!=0)  MPI_Send(&min,1,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
     else {
       for(rank=1; rank<nTasks; rank++) {
         MPI_Recv(&recv,1,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
         min+=recv;
       }
     }
     MPI_Bcast(&min,1,MPI_DOUBLE,0,MPI_COMM_WORLD);   

     cnt++;
     if(cnt%10==0 && myrank==0) 
        printf("cnt=%d, min=%g, error=%g,w=%g,rjac=%g,nn=%d\n",cnt,min,error,w,rjac,nn);

   }
}


