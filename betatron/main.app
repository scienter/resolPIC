#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "parameter.h"


void main(int argc, char *argv[])
{
   double thX,thY,th,dThX,dThY,energy,cosPhi,sinPhi;
   int id,core,i,j,k,n,cnt,step,iteration,numFreq;
   double alpha,f;
   int *idList,*coreList;
   FILE *in,*out;
   char fileName[100];
   void solveEnergy(Parameter *D,int idList,int coreList,int species,double th,double cosPhi,double sinPhi);
   void solveSpectrum(Parameter *D,int idList,int coreList,int species,double freq);
   void parameterSetting(Parameter *D,char *input);
   Parameter D;

   if(argc < 1)   { 
       printf("betatron input_filed\n");
       exit(0);  
   }

   parameterSetting(&D,argv[1]);

   double e=1.6e-19; 	//4.8e-10;
   double c=3e8;
   double pi=3.14;
   double hbar=6.582119514e-16;
   double invHbar=1.0/hbar;
   alpha=e*e/4.0/pi/c*D.superP;

   cnt=0;
   in = fopen("idList","r");
//     fgets(str,100,in);
   while(fscanf(in,"%d %d",&id,&core)!=EOF)
     cnt++;
   fclose(in);

   idList=(int *)malloc(cnt*sizeof(int ));
   coreList=(int *)malloc(cnt*sizeof(int ));
   D.B=(double *)malloc(3*sizeof(double ));
   D.C=(double *)malloc(3*sizeof(double ));


   in = fopen("idList","r");
   for(i=0; i<cnt; i++)
     fscanf(in,"%d %d",&idList[i],&coreList[i]);
   fclose(in);

   dThX=(D.maxThX-D.minThX)/((double)D.numThX);     
   dThY=(D.maxThY-D.minThY)/((double)D.numThY);     
   numFreq=(int)((D.maxE-D.minE)/D.dE);

   D.data=(double **)malloc(D.numThX*sizeof(double *));
   for(i=0; i<D.numThX; i++)
     D.data[i]=(double *)malloc(D.numThY*sizeof(double ));
   for(i=0; i<D.numThX; i++)
     for(j=0; j<D.numThY; j++)
       D.data[i][j]=0.0;
   D.freq=(double *)malloc(numFreq*sizeof(double ));
   for(i=0; i<numFreq; i++)
     D.freq[i]=0.0;

   switch (D.dimension)  {
   case 2:
/*
     //here is to calculate energy.
     for(i=0; i<D.numThX; i++)
     {
       thX=D.minThX+i*dThX;
       for(j=0; j<D.numThY; j++)
       {
         thY=D.minThY+j*dThY;
         th=sqrt(thX*thX+thY*thY);
         cosPhi=thX/th;
         sinPhi=thY/th;
         for(n=0; n<cnt; n++)
         {
           for(k=0; k<3; k++)
           {
             D.B[k]=0.0; 
             D.C[k]=0.0; 
           }
           solveEnergy(&D,idList[n],coreList[n],D.species,th,cosPhi,sinPhi);
           D.data[i][j]+=D.B[0];
         }
         printf("thX=%g, thY=%g, is done\n",thX,thY);
       }
     }
        
     sprintf(fileName,"energy");
     out = fopen(fileName,"w");
     for(i=0; i<D.numThX; i++)
     {
       thX=D.minThX+i*dThX;
       for(j=0; j<D.numThY; j++)
       {
         thY=D.minThY+j*dThY;
         fprintf(out,"%g %g %g\n",thX,thY,D.data[i][j]*alpha);
       }
       fprintf(out,"\n");
     }
     fclose(out);    
*/ 
     //here is to spectrum calculation.
     for(i=0; i<numFreq; i++)
     {
       f=(D.minE+i*D.dE)*invHbar;
       for(k=0; k<3; k++)
       {
         D.B[k]=0.0; 
         D.C[k]=0.0; 
       }
       for(n=0; n<cnt; n++)
         solveSpectrum(&D,idList[n],coreList[n],D.species,f);
       D.freq[i]=D.B[0]*D.B[0]+D.B[1]*D.B[1]+D.B[2]*D.B[2]
                +D.C[0]*D.C[0]+D.C[1]*D.C[1]+D.C[2]*D.C[2];           
       printf("%d/%d is done\n",i,numFreq); 
     }
     sprintf(fileName,"spectrum");
     out = fopen(fileName,"w");
     for(i=0; i<numFreq; i++)
     {
       energy=D.minE+i*D.dE;
       fprintf(out,"%g %g\n",energy,D.freq[i]);       
     }
     fclose(out);    
     break;
   }

   for(i=0; i<D.numThX; i++)
     free(D.data[i]);
   free(D.data);
   free(D.freq);

   free(idList);
   free(coreList);
   free(D.B);
   free(D.C);
}

void solveSpectrum(Parameter *D,int idList,int coreList,int species,double freq)
{
  double t,x,y,z,ux,uy,uz,gamma,phase,prevVx,prevVy,prevVz,step;
  double sinTh,cosTh,dt,n_dot_beta,n_dot_Dbeta,denomitor,t0;
  double vx,vy,vz,DbetaX,DbetaY,DbetaZ,betaX,betaY,betaZ;
  double RR,RRx,RRy,RRz,Rx,Ry,Rz,nx,ny,nz;
  double A[3];
  int iteration,i,id,core;
  char fileName[100];
  FILE *in=NULL;
  
  dt=D->dt;
  sprintf(fileName,"%dTrack%d_%d",species,idList,coreList);
  in = fopen(fileName,"r");
  if(in==NULL)
    exit(0);
  else
  {
    iteration=0;
    Rx=D->det_R;
    Ry=D->det_x;
    Rz=D->det_y;
    fscanf(in,"%lf %lf %lf %lf %lf %lf %d %d %lf",&x,&y,&z,&ux,&uy,&uz,&id,&core,&step);
    gamma=sqrt(1.0+ux*ux+uy*uy+uz*uz);
    prevVx=ux/gamma;
    prevVy=uy/gamma;
    prevVz=uz/gamma; //uz/gamma;
//    t0=step*dt;
    while(fscanf(in,"%lf %lf %lf %lf %lf %lf %d %d %lf",&x,&y,&z,&ux,&uy,&uz,&id,&core,&step)!=EOF)
    {
      t=step*dt;
      gamma=sqrt(1.0+ux*ux+uy*uy+uz*uz);
      vx=ux/gamma;
      vy=uy/gamma;
      vz=uz/gamma; //uz/gamma;
  
      DbetaX=(vx-prevVx)/dt;
      DbetaY=(vy-prevVy)/dt;
      DbetaZ=(vz-prevVz)/dt;
      betaX=(vx+prevVx)*0.5;
      betaY=(vy+prevVy)*0.5;
      betaZ=(vz+prevVz)*0.5;

      nx=1.0;
      ny=0.0;
      nz=0.0;
      n_dot_beta=nx*betaX+ny*betaY+nz*betaZ;
      n_dot_Dbeta=nx*DbetaX+ny*DbetaY+nz*DbetaZ;
      denomitor=(1.0-n_dot_beta)*(1.0-n_dot_beta);
      A[0]=n_dot_Dbeta*(nx-betaX)-DbetaX*(1.0-n_dot_beta)/denomitor*dt;
      A[1]=n_dot_Dbeta*(ny-betaY)-DbetaY*(1.0-n_dot_beta)/denomitor*dt;
      A[2]=n_dot_Dbeta*(nz-betaZ)-DbetaZ*(1.0-n_dot_beta)/denomitor*dt;
  
      phase=freq*(t-x/(3.0e8));

      for(i=0; i<3; i++)
      {
        D->B[i]+=A[i]*cos(phase);
        D->C[i]+=A[i]*sin(phase);
      }

      prevVx=vx;
      prevVy=vy;
      prevVz=vz;
      iteration++;
    }
  }
  fclose(in);
}


void solveEnergy(Parameter *D,int idList,int coreList,int species,double th,double cosPhi,double sinPhi)
{
  double t,x,y,z,ux,uy,uz,gamma,phase,prevVx,prevVy,prevVz,step,aveT,aveX;
  double sinTh,cosTh,dt,n_dot_beta,n_dot_Dbeta,denomitor5,tmp;
  double vx,vy,vz,DbetaX,DbetaY,DbetaZ,betaX,betaY,betaZ,denomitor,x0,invR,R;
  double A[3];
  int iteration,i,id,core;
  char fileName[100];
  FILE *in=NULL;
  
  dt=D->dt;
  cosTh=cos(th);
  sinTh=sin(th);
  x0=250e-6;
  R=3.0;
  sprintf(fileName,"%dTrack%d_%d",species,idList,coreList);
  in = fopen(fileName,"r");
  if(in==NULL)
    exit(0);
  else
  {
    prevVx=prevVy=prevVz=0.0;
    iteration=0;
//    aveT=0.5*(25500.0+40000.0);
//    aveX=0.5*(255e-6+378e-6);
    while(fscanf(in,"%lf %lf %lf %lf %lf %lf %d %d %lf",&x,&y,&z,&ux,&uy,&uz,&id,&core,&step)!=EOF)
    {
      t=step*dt;
      gamma=sqrt(1.0+ux*ux+uy*uy+uz*uz);
      vx=ux/gamma;
      vy=uy/gamma;
      vz=uz/gamma; //uz/gamma;
  
      DbetaX=(vx-prevVx)/dt;
      DbetaY=(vy-prevVy)/dt;
      DbetaZ=(vz-prevVz)/dt;
      betaX=(vx+prevVx)*0.5;
      betaY=(vy+prevVy)*0.5;
      betaZ=(vz+prevVz)*0.5;

      denomitor=(1.0-betaX*cosTh-betaY*sinTh*cosPhi-betaZ*sinTh*sinPhi);
      denomitor5=denomitor*denomitor*denomitor*denomitor*denomitor;
      n_dot_Dbeta=DbetaX*cosTh+DbetaY*sinTh*cosPhi+DbetaZ*sinTh*sinPhi;
      n_dot_beta=betaX*cosTh+betaY*sinTh*cosPhi+betaZ*sinTh*sinPhi;
      A[0]=n_dot_Dbeta*(cosTh-betaX)-DbetaX*(1.0-n_dot_beta);
      A[1]=n_dot_Dbeta*(sinTh*cosPhi-betaY)-DbetaY*(1.0-n_dot_beta);
      A[2]=n_dot_Dbeta*(sinTh*sinPhi-betaZ)-DbetaZ*(1.0-n_dot_beta);

      tmp=A[0]*A[0]+A[1]*A[1]+A[2]*A[2];
      D->B[0]+=tmp/denomitor5*dt;

      prevVx=vx;
      prevVy=vy;
      prevVz=vz;
      iteration++;
    }
  }
  fclose(in);
}
