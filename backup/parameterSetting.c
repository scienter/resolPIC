#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <gsl/gsl_qrng.h>

double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

int whatONOFF(char *str);
int whatLaserMode(char *str);
int whatLoadingMethod(char *str);

void parameterSetting(Domain *D,External *Ext, char *input)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *in=NULL;
   int FindParameters();
   int findLoadParameters();
   int findLaserParameters();
   int whatSaveMode();
   int whatFieldType();
   double x,y,z,px,py,pz,gamma;
   double positionX,factor,pMinX,pMaxX,pPosition;
   double normalB,normalE,Ex,Ey,Ez,Bx,By,Bz,dyoverdx,dzoverdx;
   char str[100],name[100],fileName[100];
   int rank,minT,maxT,tmpInt,fail=0,cnt;
   int i,j,k,n,numProbeX,numProbeY,numProbeZ,probeType,id,core,species;
   double lambda,tmpDouble,probeDx,probeDy,probeDz,maxProbeX,minProbeX,maxProbeY,minProbeY,maxProbeZ,minProbeZ,tmp;
   LoadList *LL, *New;
   LaserList *L, *LNew;


   //initially
   if(FindParameters("Domain",1,"dimension",input,str)) 
     D->dimension=atoi(str);
   else  {
      printf("in [Domain], dimension=?  (1:1D, 2:2D, 3:3D)\n");
      fail=1;
   }
   //Boost frame
   if(FindParameters("Domain",1,"boost_gamma",input,str)) D->gamma=atof(str);
   else D->gamma=1;
   if(FindParameters("Domain",1,"boost_ion",input,str)) D->boostIon=whatONOFF(str);
   else D->boostIon=ON;
   if(D->gamma>1)   D->boostOn=ON;
   else             D->boostOn=OFF;
   D->beta=sqrt(1-1.0/D->gamma/D->gamma);
   if(FindParameters("Domain",1,"lambda",input,str)) {
      D->lambda=atof(str);
      D->lambda*=D->gamma*(1+D->beta);
   } else  {
      printf("In [Domain], lambda=? [m]\n");
      fail=1;
   }

   if(FindParameters("Domain",1,"L",input,str)) D->L=atoi(str);
   else  {
      printf("in [Domain], L=?  (Sorry. Fix as L=1)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"M",input,str)) D->M=atoi(str);
   else  {
      printf("in [Domain], M=?  (y directionally dividing number)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"N",input,str)) D->N=atoi(str);
   else  {
      printf("in [Domain], N=?  (z directionally dividing number)\n");
      fail=1;
   }
   if(D->dimension==1)  {
     D->M=1;
     D->N=1;
   }  else if(D->dimension==2) 
     D->N=1;
   if(D->dimension==2)  {
     if(D->M*D->L!=nTasks)  {
       printf("L=%d, M=%d, check the values of L and M.\n",D->L,D->M);
       fail=1;
     }
   }  else if(D->dimension==3)  {
     if(D->M*D->N*D->L!=nTasks)  {
       printf("L=%d, M=%d, N=%d, check the values of L, M, and N.\n",D->L,D->M,D->N);
       fail=1;
     }
   }

   //Field Type
   if(FindParameters("Domain",1,"field_type",input,str)) 
     D->fieldType=whatFieldType(str);
   else  {
      printf("in [Domain], field_type=?  (Split/Yee/Pukhov)\n");
      fail=1;
   }
   //Current Type
   if(FindParameters("Domain",1,"current_order",input,str)) D->currentType=atoi(str);
   else  
      D->currentType=1;
   if(FindParameters("Domain",1,"interpolation_order",input,str)) D->interpolationType=atoi(str);
   else 
      D->interpolationType=1;

   //Domain parameter setting
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], maxStep=? \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_step",input,str)) D->saveStep=atoi(str);
   else  {
      printf("In [Domain], save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  {
      printf("In [Domain], save_start=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"center_save_step",input,str)) D->centerStep=atoi(str);
   else  {
      printf("In [Domain], center_save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"period_boundary_X",input,str)) D->periodX=whatONOFF(str);
   else  D->periodX=OFF;
   if(FindParameters("Domain",1,"period_boundary_Y",input,str)) D->periodY=whatONOFF(str);
   else  D->periodY=OFF;
   if(FindParameters("Domain",1,"period_boundary_Z",input,str)) D->periodZ=whatONOFF(str);
   else  D->periodZ=OFF;

   //save options
   if(FindParameters("Save",1,"field_format",input,str)) 
     D->saveFieldMode=whatSaveMode(str);
   else  
     D->saveFieldMode=TXT;
   if(FindParameters("Save",1,"particle_format",input,str)) 
     D->saveParticleMode=whatSaveMode(str);
   else  
     D->saveParticleMode=TXT;
   if(FindParameters("Save",1,"density_format",input,str)) 
     D->saveDensityMode=whatSaveMode(str);
   else  
     D->saveDensityMode=TXT;
   if(FindParameters("Save",1,"dump_format",input,str)) 
     D->saveDumpMode=whatSaveMode(str);
   else  
     D->saveDumpMode=TXT;
   if(FindParameters("Save",1,"dump_save",input,str)) 
     D->dumpSave=whatONOFF(str);
   else  
     D->dumpSave=OFF;
   if(FindParameters("Save",1,"dump_start",input,str)) 
     D->dumpStart=atoi(str);
   else  
     D->dumpStart=D->saveStart;
   if(FindParameters("Save",1,"dump_save_step",input,str)) 
     D->dumpSaveStep=atoi(str);
   else  
     D->dumpSaveStep=D->saveStep;
   D->dumpStep=0;
   if(FindParameters("Save",1,"field_save",input,str)) 
     D->fieldSave=whatONOFF(str);
   else  
     D->fieldSave=ON;
   if(FindParameters("Save",1,"raman_save",input,str)) 
     D->ramanSave=whatONOFF(str);
   else  
     D->ramanSave=ON;
   if(FindParameters("Save",1,"particle_save",input,str)) 
     D->particleSave=whatONOFF(str);
   else  
     D->particleSave=ON;
   if(FindParameters("Save",1,"density_save",input,str)) 
     D->densitySave=whatONOFF(str);
   else  
     D->densitySave=ON;
   if(FindParameters("Save",1,"test_save",input,str)) 
     D->testSave=whatONOFF(str);
   else  
     D->testSave=ON;

   if(FindParameters("Save",1,"resolution_change",input,str)) 
      D->resolChange=whatONOFF(str);
   else  D->resolChange=OFF;
   if(FindParameters("Save",1,"resolution_high",input,str)) 
      D->resolHigh=whatONOFF(str);
   else  D->resolHigh=OFF;
   if(FindParameters("Save",1,"resolution_low",input,str)) 
      D->resolLow=whatONOFF(str);
   else  D->resolLow=OFF;
   if(FindParameters("Save",1,"resolution_change_step",input,str)) 
      D->resolStep=atoi(str);
   else  D->resolStep=D->maxTime;
   if(FindParameters("Save",1,"resolution_rate_X",input,str)) 
      D->resolX=atoi(str);
   else  D->resolX=1;
   if(FindParameters("Save",1,"resolution_rate_Y",input,str)) 
      D->resolY=atoi(str);
   else  D->resolY=1;
   if(FindParameters("Save",1,"resolution_rate_Z",input,str)) 
      D->resolZ=atoi(str);
   else  D->resolZ=1;
   if(FindParameters("Domain",1,"minX",input,str)) 
   {
      D->minX=atof(str)/D->lambda;
      D->minX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("In [Domain], minX=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str)) 
   {
      D->maxX=atof(str)/D->lambda;
      D->maxX*=D->gamma*(1+D->beta);
   }   else  {
      printf("In [Domain], maxX=? [m].\n");
      fail=1;
   }

   if(D->dimension>1)
   {
     if(FindParameters("Domain",1,"minY",input,str)) 
       D->minY=atof(str)/D->lambda;
     else  {
       printf("In [Domain], minY=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"maxY",input,str)) 
       D->maxY=atof(str)/D->lambda;
     else  {
       printf("In [Domain], maxY=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"dy_over_dx",input,str)) 
       dyoverdx=atof(str);
     else  {
       printf("In [Domain], dy_over_dx=?  [dy/dx]\n");
       fail=1;
     }
   }
   else	;
   if(D->dimension>2)
   {
     if(FindParameters("Domain",1,"minZ",input,str)) 
       D->minZ=atof(str)/D->lambda;
     else  {
       printf("In [Domain], minZ=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"maxZ",input,str)) 
       D->maxZ=atof(str)/D->lambda;
     else  {
       printf("In [Domain], maxZ=? [m].\n");
       fail=1;
     }
     if(FindParameters("Domain",1,"dz_over_dx",input,str)) 
       dzoverdx=atof(str);
     else  
       dzoverdx=dyoverdx;
   }
   else	;
   if(FindParameters("Domain",1,"moving_domain",input,str)) D->moving=whatONOFF(str);
   else  {
      printf("In [Domain], moving_domain=? [ON/OFF].\n");
      fail=1;
   }  
   if(FindParameters("Domain",1,"division_lambda",input,str)) 
   {
      D->divisionLambda=atof(str);
   }
   else  {
      printf("In [Domain], divisionLambda=? [number of devided wavelength]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dt_ratio",input,str)) 
      D->dtRatio=atof(str);
   else  {
      printf("In [Domain], dt_ratio=? [<1.0]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"moving_velocity",input,str))
     D->movingV=atof(str);
   else D->movingV=D->dtRatio; 


   //External field parameter setting
   if(FindParameters("External",1,"Ex",input,str)) Ex=atof(str);
   else  {
      printf("In [External], Ex=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ey",input,str)) Ey=atof(str);
   else  {
      printf("In [External], Ey=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ez",input,str)) Ez=atof(str);
   else  {
      printf("In [External], Ez=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bx",input,str)) Bx=atof(str);
   else  {
      printf("In [External], Bx=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"By",input,str)) By=atof(str);
   else  {
      printf("In [External], By=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bz",input,str)) Bz=atof(str);
   else  {
      printf("In [External], Bz=? [Tesla]\n");
      fail=1;
   }

   //pml
   if(FindParameters("PML",1,"pml",input,str)) 
     D->pmlOn=whatONOFF(str);
   else  
     D->pmlOn=OFF;
   if(FindParameters("PML",1,"pml_start",input,str)) 
     D->pmlStart=atoi(str);
   else  {
      printf("In [PML], pml_start=? [what step]\n");
      fail=1;
   }
   if(FindParameters("PML",1,"right_pml_cells",input,str))
     D->pmlCellRight=atoi(str);
   else  D->pmlCellRight=1;
   if(FindParameters("PML",1,"left_pml_cells",input,str))
     D->pmlCellLeft=atoi(str);
   else  D->pmlCellLeft=1;
   if(FindParameters("PML",1,"up_pml_cells",input,str))
     D->pmlCellUp=atoi(str);
   else  D->pmlCellUp=1;
   if(FindParameters("PML",1,"down_pml_cells",input,str))
     D->pmlCellDown=atoi(str);
   else  D->pmlCellDown=1;
   if(FindParameters("PML",1,"front_pml_cells",input,str)) D->pmlCellFront=atoi(str);
   else  D->pmlCellFront=1;
   if(FindParameters("PML",1,"back_pml_cells",input,str)) D->pmlCellBack=atoi(str);
   else  D->pmlCellBack=1;
   if(FindParameters("PML",1,"pml_r",input,str)) D->pmlr=atof(str);
   else  { printf("In [PML], pml_r=? (retarding length)\n"); fail=1; }
   if(FindParameters("PML",1,"pml_d",input,str)) D->pmld=atof(str);
   else  { printf("In [PML], pml_d=? (damping length)\n"); fail=1; }

   //field ionization
   if(FindParameters("Domain",1,"field_ionization",input,str))
     D->fieldIonization=whatONOFF(str);
   else  D->fieldIonization=OFF;

   //additional Domain parameters  
   D->iteration=0;
   D->dx=1.0/D->divisionLambda;
   D->nx=((int)((D->maxX-D->minX)/D->dx));
   D->ny=1;
   D->nz=1;
   tmpDouble=D->dx/(D->gamma*(1+D->beta));
   D->dy=D->dz=1.0;
   if(D->fieldType==Split) {
     D->dt=D->dx; 
     D->dtRatio=1.0;
//     D->movingV=1.0;
     D->shiftDuration=D->maxStep+D->resolX;
   } else {
     D->dt=D->dx*D->dtRatio;
     D->shiftDuration=(int)(1.0/(1.0-D->movingV));
   }
   D->shiftStart=0;

   if(D->dimension>1)   {
     D->dy=D->dx*dyoverdx;
     D->ny=((int)((D->maxY-D->minY)/D->dy));
   }
   else	;
   if(D->dimension>2)   {
     D->dz=D->dx*dzoverdx;
//     D->dz=tmpDouble*dzoverdx;
//     if(D->dz<=D->dx*1.5)   {
//       printf("dz_over_dx is too low. It must be over than %g.\n", D->gamma*(1+D->beta)*1.5);
//       fail=1;
//     }
//     else	;
     D->nz=((int)((D->maxZ-D->minZ)/D->dz));
   }
   else	;  

   D->minXDomain=D->minYDomain=D->minZDomain=0;
   D->maxXDomain=D->nx;
   D->omega=2*pi*velocityC/D->lambda;
   if(D->boostOn==ON)   {
     D->minXSub=-D->nx;
     if(D->dimension>1)
       D->minYDomain=(int)(D->minY/D->dy);
     else	;
     if(D->dimension>2)
       D->minZDomain=(int)(D->minZ/D->dz);
     else	;
   }
   else   {
     D->minXSub=0;
     if(D->dimension>1)
       D->minYDomain=(int)(D->minY/D->dy);
     if(D->dimension>2)
       D->minZDomain=(int)(D->minZ/D->dz);
   }

   if(myrank==0)
   {
     printf("dx=%g,dt=%g\n",D->dx,D->dt);
     if(D->dimension==1)     {
       printf("gamma=%g, dx=%g, dt=%g, divisionLambda=%g\n",D->gamma,D->dx,D->dt,D->divisionLambda);
     } else if(D->dimension==2)     {
       printf("gamma=%g, dx=%g, dy=%g, dt=%g, divisionLambda=%g\n",D->gamma,D->dx,D->dy,D->dt,D->divisionLambda);
     } else if(D->dimension==3)     {
       printf("gamma=%g, dx=%g, dy=%g, dz=%g, dt=%g, divisionLambda=%g\n",D->gamma,D->dx,D->dy,D->dz,D->dt,D->divisionLambda);
     } else ;
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);

   //additional Boost parameters
   factor=D->gamma*(1+D->beta);
   D->minT=(int)(D->maxStep/factor/factor); 	//boost frame iteration
   D->maxT=(int)((D->maxStep+D->beta*D->nx)/(1+D->beta)-factor*D->gamma*D->minT*D->beta);	//boost frame iteration
   if(myrank==0)
     printf("maxT=%d, nx=%d, ny=%d, nz=%d\n",D->maxT,D->nx,D->ny,D->nz);
   else	;

   //additional external field parameters
   normalB=eMass*D->omega/(-eCharge);
   normalE=normalB*velocityC;
   Ext->E1=Ex/normalE;
   Ext->E2=Ey/normalE;
   Ext->E3=Ez/normalE;
   Ext->B1=Bx/normalB;
   Ext->B2=By/normalB;
   Ext->B3=Bz/normalB;

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input)) 
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

   //Plasma parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
      exit(0);
   else	;

}

int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   double positionX,positionY,positionZ,k0,tmp;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"loading_method",input,str)) L->loadMethod=whatLoadingMethod(str);
   else  L->loadMethod=0;

   if(L->loadMethod)
   {
     if(FindParameters("Laser",rank,"polarity",input,str)) L->polarity=atoi(str);
	  else  L->polarity=0;
     if(FindParameters("Laser",rank,"wavelength",input,str))  {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     }
     else  L->lambda=D->lambda;
     if(FindParameters("Laser",rank,"a0",input,str)) L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  { printf("in [Laser], focus=?  [m]\n"); fail=1; }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  positionX=0;
     if(FindParameters("Laser",rank,"loadPositionY",input,str)) positionY=atof(str);
     else  positionY=0;
     if(FindParameters("Laser",rank,"loadPositionZ",input,str)) positionZ=atof(str);
     else  positionZ=0;
     if(FindParameters("Laser",rank,"add",input,str)) L->add=whatONOFF(str);
     else  L->add=OFF;

     if(FindParameters("Laser",rank,"mode",input,str)) L->mode=whatLaserMode(str);
     else  L->mode=Gaussian;
     if(L->mode==Gaussian) {	 
       if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
       else  {
         printf("in [Laser], rU=? [# of basic wavelength]\n");
         fail=1;
       }
       if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
       else  {
         printf("in [Laser], rD=? [# of basic wavelength]\n");
         fail=1;
       }
       if(FindParameters("Laser",rank,"retard",input,str)) L->retard=atof(str);
       else  {
         printf("in [Laser], retard=? [# of basic wavelength]\n");
         fail=1;
       }
       if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
       else  L->flat=0.0;
       if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
       else  L->direction=1;
       if(FindParameters("Laser",rank,"gdd",input,str)) L->gdd=atof(str);
       else  L->gdd=0.0;
       if(D->dimension>1) {
         if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
         else  {
           printf("in [Laser], beamWaist=?  [m]\n");
           fail=1;
         }
       } else ;
       if(D->dimension==3) {
         if(FindParameters("Laser",rank,"elliptic",input,str)) L->elliptic=atof(str);
         else  L->elliptic=1.0;
       } else ;
     } 
     //SSTF setting
     else if(L->mode==SSTF) {
       if(FindParameters("Laser",rank,"lens_focus",input,str)) L->lensFocus=atof(str);
       else  { printf("in [Laser], lens_focus=?  [m]\n"); fail=1; }
       if(FindParameters("Laser",rank,"SSTF_focused_sigma_x",input,str)) L->sigX_t=atof(str);   //based on the field amplitude
       else  { printf("in [Laser], SSTF_focused_sigma_x=?  [m]\n"); fail=1; }
       if(FindParameters("Laser",rank,"SSTF_focused_sigma_y",input,str)) L->sigY_t=atof(str);
       else  { printf("in [Laser], SSTF_focused_sigma_y=?  [m]\n"); fail=1; }
       if(FindParameters("Laser",rank,"SSTF_linear_coef_y",input,str)) L->alphaYSSTF=atof(str);
       else  { printf("in [Laser], SSTF_linear_coef_y=?  [m/Hz]\n"); fail=1; }
       tmp=4.0*L->lensFocus*L->lensFocus/k0/k0/L->sigY_t/L->sigY_t-4.0*L->alphaYSSTF*L->alphaYSSTF*velocityC*velocityC/L->sigX_t/L->sigX_t;
       if(tmp<0.0) { printf("test value should be plus. test=%g\n",tmp); fail=1; }  else ;
       L->sigYSSTF=sqrt(tmp);
       if(D->dimension==3) {
         if(FindParameters("Laser",rank,"SSTF_focused_sigma_z",input,str)) L->sigZ_t=atof(str);
         else  { printf("in [Laser], SSTF_focused_sigma_z=?  [m]\n"); fail=1; }
         if(FindParameters("Laser",rank,"SSTF_linear_coef_z",input,str)) L->alphaZSSTF=atof(str);
         else  { printf("in [Laser], SSTF_linear_coef_z=?  [m/Hz]\n"); fail=1; }			
         tmp=4.0*L->lensFocus*L->lensFocus/k0/k0/L->sigZ_t/L->sigZ_t-4.0*L->alphaZSSTF*L->alphaZSSTF*velocityC*velocityC/L->sigX_t/L->sigX_t;
         if(tmp<0.0) { printf("test value should be plus. test=%g\n",tmp); fail=1; }  else ;
         L->sigZSSTF=sqrt(tmp);
         if(L->alphaZSSTF!=0.0 && L->alphaYSSTF!=0.0) fail=1; else ;
       } else ;
       L->sigOmegaSSTF=2.0*velocityC/L->sigX_t;
     } else ;

     //additional laser parameters
     k0=2.0*M_PI/L->lambda;
     L->omega=2*M_PI*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dx));   
     L->loadPointY=((int)(positionY/D->lambda/D->dy));   
     L->loadPointZ=((int)(positionZ/D->lambda/D->dz));   
     L->rayleighLength=pi/(L->lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;
     if(L->loadMethod==Shot) D->shiftStart=D->nx; else ;

     if(fail==1) exit(0); else ;
   }
   return L->loadMethod;
}

int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   LoadList *New;
   DefPtcl *Def;
   int whatSpecies();
   int whatPlasmaType();
   double whatMass(int species);
   double pointPosition,wp,pDt,tmpX,tmpY,tmpZ,tmpR;
   int whatCharge();
   int whatFunctionMode();
   int whatDefineMode();
   char name[100], str[100];
   int i,n,cnt,species,fail=0,tmpInt,flag;
   double tmp,max,min,randX,randY,randZ;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Plasma",rank,"type",input,name)) 
   {
     LL->type = whatPlasmaType(name);
     if(D->boostOn==ON)
       LL->type = BoostFrame;
   }
   else LL->type=0;

   if(LL->type>0)
   {
      if(FindParameters("Plasma",rank,"density",input,str)) 
      {
         LL->density=atof(str);
         LL->density*=D->gamma;
      }
      else  {
         printf("in [Plasma], density=? [m-3]\n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"pair",input,str)) LL->pair=whatONOFF(str);
      else  LL->pair=OFF;


      if(FindParameters("Plasma",rank,"species",input,name)) 
         species = whatSpecies(name);
      else  species = 0;
      LL->species=species;
      if(FindParameters("Plasma",rank,"numberInCell",input,str)) 
         LL->numberInCell=atoi(str);
      else  {
         printf("in [Plasma], numberInCell=? \n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"startIndex",input,str)) 
         LL->index=atoi(str);
      else  
         LL->index=0;
      if(FindParameters("Plasma",rank,"temperature",input,str)) LL->temperature=atof(str);
      else   LL->temperature=0.0;	
      if(FindParameters("Plasma",rank,"min_px",input,str)) LL->givenMinPx = atof(str);
      else  LL->givenMinPx = -1e9;
      LL->mass=whatMass(species);
      LL->charge=whatCharge(species);
      LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
//      LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy*D->lambda*D->dz/LL->numberInCell;

      //setting ionization
      ionizationSetup(LL,LL->species);
      if(FindParameters("Plasma",rank,"min_a0",input,str)) 
         LL->givenMinA0 = atof(str);
      else  LL->givenMinA0 = 0.0  ;

      srand(1*(myrank+1));
      switch (LL->type)  {
      case Polygon :
        if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
        else  {
          printf("in [Plasma], Xnodes=?\n");
          printf("Each nodes indicates the point of plasma density changing.\n");
          fail=1;
        }
        if(LL->xnodes>0)
        {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)
          {
            sprintf(name,"X%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xpoint[i] = atof(str)/D->gamma/D->lambda/D->dx;
            else 
            { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xn[i] = atof(str);
            else 
            { printf("Xn%d should be defined.\n",i);  fail=1; } 
          }
        }
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"Ynodes",input,str)) LL->ynodes=atoi(str);
          else  {
            printf("in [Plasma], Ynodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          if(LL->ynodes>0)
          {
            LL->ypoint = (double *)malloc(LL->ynodes*sizeof(double));
            LL->yn = (double *)malloc(LL->ynodes*sizeof(double));   
            for(i=0; i<LL->ynodes; i++)
            {
              sprintf(name,"Y%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) {
                LL->ypoint[i] = atof(str)/D->lambda/D->dy;
              }
              else 
              { printf("Y%d should be defined.\n",i);  fail=1; }
 
              sprintf(name,"Yn%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) 
                LL->yn[i] = atof(str);
              else 
              { printf("Yn%d should be defined.\n",i);  fail=1; } 
            }
          }
        }
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"Znodes",input,str)) LL->znodes=atoi(str);
          else  {
            printf("in [Plasma], Znodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          LL->zpoint = (double *)malloc(LL->znodes*sizeof(double));
          LL->zn = (double *)malloc(LL->znodes*sizeof(double));   
          for(i=0; i<LL->znodes; i++)
          {
            sprintf(name,"Z%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) {
              LL->zpoint[i] = atof(str)/D->lambda/D->dz;
            }
            else 
            { printf("Z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Zn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->zn[i] = atof(str);
            else 
            { printf("Zn%d should be defined.\n",i);  fail=1; } 
          }
        }
        if(FindParameters("Plasma",rank,"centerX",input,str))  
          LL->centerX=atof(str)/D->lambda/D->dx;
        else   LL->centerX=0.0;	
        LL->centerY=0.0;	
        LL->centerZ=0.0;	
        if(FindParameters("Plasma",rank,"gauss_coef_X",input,str))  
          LL->gaussCoefX=atof(str)/D->lambda/D->dx;
        else   LL->gaussCoefX=1.0;
        if(FindParameters("Plasma",rank,"poly_coef_X",input,str))  
          LL->polyCoefX=atof(str)/D->lambda/D->dx;
        else   LL->polyCoefX=0.0;	
        if(FindParameters("Plasma",rank,"function_mode_X",input,str))  
          LL->modeX=whatFunctionMode(str);
        else   LL->modeX=0;	
        if(FindParameters("Plasma",rank,"function_mode_YZ",input,str))  
          LL->modeYZ=whatFunctionMode(str);
        else   LL->modeYZ=0;	
        if(D->dimension>1)
        {	
          if(FindParameters("Plasma",rank,"centerY",input,str))  
            LL->centerY=atof(str)/D->lambda/D->dy;
          else   LL->centerY=0.0;	
          if(FindParameters("Plasma",rank,"gauss_coef_YZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/D->dy;
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"poly_coef_YZ",input,str))  
            LL->polyCoefYZ=atof(str)*D->lambda*D->dy*D->lambda*D->dy;
          else   LL->polyCoefYZ=0.0;	
        }
        else if(D->dimension>2)
        {	
          if(FindParameters("Plasma",rank,"centerZ",input,str))  
            LL->centerZ=atof(str)/D->lambda/D->dz;
          else   LL->centerZ=0.0;	
          if(FindParameters("Plasma",rank,"gauss_coef_YZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/sqrt(D->dy*D->dy+D->dz*D->dz);
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"poly_coef_YZ",input,str))  
            LL->polyCoefYZ=atof(str)*D->lambda*D->lambda*(D->dy*D->dy+D->dz*D->dz);
          else   LL->polyCoefYZ=0.0;	
        }
        else 	;
        break;

      case Defined :
//        srand(time(NULL)*(myrank+1));
        if(FindParameters("Plasma",rank,"define_mode",input,str))  
          LL->defineMode=whatDefineMode(str);
        else   LL->defineMode=byNumber;	
        if(FindParameters("Plasma",rank,"number_defined",input,str))  
          LL->numDefined=atoi(str);
        else   LL->numDefined=0;	
        if(LL->defineMode==byDensity)
        {
          if(FindParameters("Plasma",rank,"minX",input,str))  
            LL->minX=atof(str)/D->lambda;
          else   {   printf("minX=?  [m]\n");  exit(0);   }	
          if(FindParameters("Plasma",rank,"maxX",input,str))  
            LL->maxX=atof(str)/D->lambda;
          else   {   printf("maxX=?  [m]\n");  exit(0);   }	
          LL->minY=0.0; LL->maxY=0.0;
          LL->minZ=0.0; LL->maxZ=0.0;
          if(D->dimension>1)
          {
            if(FindParameters("Plasma",rank,"minY",input,str))  
              LL->minY=atof(str)/D->lambda;
            else   {   printf("minY=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxY",input,str))  
              LL->maxY=atof(str)/D->lambda;
            else   {   printf("maxY=?  [m]\n");  exit(0);   }	
          } else ;
          if(D->dimension>2)
          {
            if(FindParameters("Plasma",rank,"minZ",input,str))  
              LL->minZ=atof(str)/D->lambda;
            else   {   printf("minZ=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxZ",input,str))  
              LL->maxZ=atof(str)/D->lambda;
            else   {   printf("maxZ=?  [m]\n");  exit(0);   }	
          } else ;
        }
        else 	;

        if(FindParameters("Plasma",rank,"radius",input,str))  
          LL->RDef=atof(str)/D->lambda;
        else   LL->RDef=0.0;
        LL->createGuard=D->maxX;
        LL->deleteGuard=0.0;
//        LL->createGuard=D->maxX+LL->RDef*2.0;
//        LL->deleteGuard=D->maxX-LL->RDef*2.0;
        if(D->dimension==2) {
          LL->numDefPtcls=(int)(pi*LL->RDef*LL->RDef/D->dx/D->dy*LL->numberInCell);
        }
        else if(D->dimension==3)
          LL->numDefPtcls=(int)(4.0/3.0*pi*LL->RDef*LL->RDef*LL->RDef/D->dx/D->dy/D->dz*LL->numberInCell);
        else ;
  
        gsl_qrng *q1 = gsl_qrng_alloc (gsl_qrng_sobol,1);
        gsl_qrng *q2 = gsl_qrng_alloc (gsl_qrng_sobol,2);
        gsl_qrng *q3 = gsl_qrng_alloc (gsl_qrng_sobol,3);
        for(i=0; i<LL->numDefined; i++)
        {
          Def = (DefPtcl *)malloc(sizeof(DefPtcl ));
          Def->next = LL->def;
          LL->def = Def;
          if(i==0) Def->next=NULL; else ;

          if(LL->defineMode==byNumber)
          {
            if(myrank==0) {
              sprintf(name,"xPosition%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) 
                Def->xPos=atof(str)/D->lambda;
              else
              { printf("xPosition%d should be defined.\n",i); fail=1;}
            } else ;
            MPI_Bcast(&(Def->xPos),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            Def->yPos=0.0; Def->zPos=0.0;
            if(D->dimension>1) {
              if(myrank==0) {
                sprintf(name,"yPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str)) 
                  Def->yPos=atof(str)/D->lambda;
                else
                { printf("yPosition%d should be defined.\n",i); fail=1;}
              } else ;
              MPI_Bcast(&(Def->yPos),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
            } else ;
            if(D->dimension>2) {
              if(myrank==0) {
                sprintf(name,"zPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str)) 
                  Def->zPos=atof(str)/D->lambda;
                else
                { printf("zPosition%d should be defined.\n",i); fail=1;}
              } else ;
              MPI_Bcast(&(Def->zPos),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
            } else ;
          }
          else if(LL->defineMode==byDensity) {
            if(D->dimension==1) {
              random1D_sobol(&randX,q1);
              Def->xPos=LL->minX+randX*(LL->maxX-LL->minX);
              Def->yPos=0.0; Def->zPos=0.0;
            } else if(D->dimension==2) {
              random2D_sobol(&randX,&randY,q2);
              Def->xPos=LL->minX+randX*(LL->maxX-LL->minX);
              Def->yPos=LL->minY+randY*(LL->maxY-LL->minY);
              Def->zPos=0.0;
            } else  {
              random3D_sobol(&randX,&randY,&randZ,q3);
              Def->xPos=LL->minX+randX*(LL->maxX-LL->minX);
              Def->yPos=LL->minY+randY*(LL->maxY-LL->minY);
              Def->zPos=LL->minZ+randZ*(LL->maxZ-LL->minZ);
            }
          }
          else 	;
          Def->numDefPtcls=0;
          Def->flag=OFF;
        }
        gsl_qrng_free(q1);
        gsl_qrng_free(q2);
        gsl_qrng_free(q3);

        if(myrank==0) {
          if(FindParameters("Plasma",rank,"min_load_step",input,str))  
            LL->minLoadTime=atoi(str);
          else   LL->minLoadTime=0;
        } else ;
        MPI_Bcast(&(LL->minLoadTime),1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(myrank==0) {
          if(FindParameters("Plasma",rank,"max_load_step",input,str))  
            LL->maxLoadTime=atoi(str);
          else   LL->maxLoadTime=0;
        } else ;
        MPI_Bcast(&(LL->maxLoadTime),1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        break;

      case Beam :
        if(FindParameters("Plasma",rank,"loading_step",input,str)) LL->loadingStep=atoi(str);
        else                                                       LL->loadingStep=0;
	if(LL->loadingStep==0) D->shiftStart=D->nx;
	else                   D->shiftStart=0;
        if(FindParameters("Plasma",rank,"energy",input,str)) LL->energy=atof(str);
        else { printf("In [Plasma], energy=? [MeV] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"energy_spread",input,str)) LL->spread=atof(str)*1e-2;
        else { printf("In [Plasma], energy_spread=? [%%] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"norm_emittance_y",input,str)) LL->emitY=atof(str)*1e-6;
        else { printf("In [Plasma], norm_emittance_y=? [mm mrad] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"norm_emittance_z",input,str)) LL->emitZ=atof(str)*1e-6;
        else { printf("In [Plasma], norm_emittance_z=? [mm mrad] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"rms_sigma_y",input,str)) LL->sigY=atof(str);
        else { printf("In [Plasma], rms_sigma_y=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"rms_sigma_z",input,str)) LL->sigZ=atof(str);
        else { printf("In [Plasma], rms_sigma_z=? [m] \n"); fail=1;}
//        if(FindParameters("Plasma",rank,"position_x",input,str)) LL->posX=atof(str);
//        else { printf("In [Plasma], position_x=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"position_y",input,str)) LL->posY=atof(str);
        else { printf("In [Plasma], position_y=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"position_z",input,str)) LL->posZ=atof(str);
        else { printf("In [Plasma], position_z=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"focus",input,str)) LL->focus=atof(str);
        else { printf("In [Plasma], focus=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
        else  {
          printf("in [Plasma], Xnodes=?\n");
          printf("Each nodes indicates the point of plasma density changing.\n");
          fail=1;
        }
        if(LL->xnodes>0) {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          LL->xenergy = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)  {
            sprintf(name,"X%d",i);
            if(FindParameters("Plasma",rank,name,input,str))
	      LL->xpoint[i] = atof(str)/D->gamma/D->lambda/D->dx;
            else { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("Plasma",rank,name,input,str))LL->xn[i] = atof(str);
            else { printf("Xn%d should be defined.\n",i);  fail=1; } 
            sprintf(name,"Xenergy%d",i);
            if(FindParameters("Plasma",rank,name,input,str))LL->xenergy[i] = atof(str);
            else { printf("Xenergy%d should be defined.\n",i);  fail=1; } 
          }
        }
	// Poisson solver parameters
//        if(FindParameters("Plasma",rank,"rjac",input,str)) LL->rjac=atof(str);
//        else { printf("In [Plasma], rjac=? [<1.0] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"rjac_compensation",input,str)) LL->rjacRatio=atof(str);
        else { printf("In [Plasma], rjac_compensation=? [ex 0.999] \n"); fail=1;}
        LL->rjac=(cos(M_PI/(1.0*D->ny))+D->dx*D->dx/D->dy/D->dy*cos(M_PI/(1.0*D->nx)))/(1.0+D->dx*D->dx/D->dy/D->dy);
	LL->rjac*=LL->rjacRatio;

        if(FindParameters("Plasma",rank,"err_ratio",input,str)) LL->ratio=atof(str);
        else { printf("In [Plasma], err_ratio=? [(ex) 1e-5] \n"); fail=1;}

	LL->betaY=LL->sigY*LL->sigY/LL->emitY;
	LL->betaZ=LL->sigZ*LL->sigZ/LL->emitZ;
        break;

      }
   
   }	//end of if(species)

   if(fail==1)
      exit(0);

   return LL->type;
}

int whatDefineMode(char *str)
{
   if(strstr(str,"by_number")) 		return byNumber;
   else if(strstr(str,"by_density"))   	return byDensity;
   else 				return byNumber;
}

int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else 				return OFF;
}

int whatSaveMode(char *str)
{
   if(strstr(str,"TXT")) 		return TXT;
   else if(strstr(str,"HDF"))   	return HDF;
   else 				return TXT;
}

int whatFieldType(char *str)
{
   if(strstr(str,"Split")) 		return Split;
   else if(strstr(str,"Yee"))   	return Yee;
   else if(strstr(str,"Pukhov"))   	return Pukhov;
   else 				return 0;
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"Test")) 		return Test;
   else if(strstr(str,"Positron")) 		return Positron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else if(strstr(str,"NPlus0"))    return NPlus0;
   else if(strstr(str,"NPlus1"))    return NPlus1;
   else if(strstr(str,"NPlus2"))    return NPlus2;
   else if(strstr(str,"NPlus3"))    return NPlus3;
   else if(strstr(str,"NPlus4"))    return NPlus4;
   else if(strstr(str,"NPlus5"))    return NPlus5;
   else if(strstr(str,"NPlus6"))    return NPlus6;
   else if(strstr(str,"NPlus7"))    return NPlus7;
   else if(strstr(str,"AlPlus4"))   	return AlPlus4;
   else return 0;
}

int whatPlasmaType(char *str)
{
   if(strstr(str,"Polygon"))         return Polygon;
   else if(strstr(str,"Defined"))   	return Defined;
   else if(strstr(str,"Beam"))   	return Beam;
   else
   {
     printf("No Plasma type!\n"); 
     exit(0);
   }
   return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == Positron) 		return 30;
   else if(species == Test) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)           return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)           return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)           return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)           return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)           return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)           return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)           return (12.0111-6*eMassU)/eMassU;
   else if(species == NPlus0)           return (14.0064-0*eMassU)/eMassU;
   else if(species == NPlus1)           return (14.0064-1*eMassU)/eMassU;
   else if(species == NPlus2)           return (14.0064-2*eMassU)/eMassU;
   else if(species == NPlus3)           return (14.0064-3*eMassU)/eMassU;
   else if(species == NPlus4)           return (14.0064-4*eMassU)/eMassU;
   else if(species == NPlus5)           return (14.0064-5*eMassU)/eMassU;
   else if(species == NPlus6)           return (14.0064-6*eMassU)/eMassU;
   else if(species == NPlus7)           return (14.0064-7*eMassU)/eMassU;
   else if(species == AlPlus4)         return (26.9815-4*eMassU)/eMassU;
   else {  printf("Species (%d) mass not defined\n",species);  exit(0);  }
}

int whatCharge(int species)
{
   int fail;

   if(species == Electron) 		return -1;
   else if(species == Positron) 		return 1;
   else if(species == Test) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0)           return 0;
   else if(species == CPlus1)           return 1;
   else if(species == CPlus2)           return 2;
   else if(species == CPlus3)           return 3;
   else if(species == CPlus4)           return 4;
   else if(species == CPlus5)           return 5;
   else if(species == CPlus6)           return 6;
   else if(species == NPlus0)           return 0;
   else if(species == NPlus1)           return 1;
   else if(species == NPlus2)           return 2;
   else if(species == NPlus3)           return 3;
   else if(species == NPlus4)           return 4;
   else if(species == NPlus5)           return 5;
   else if(species == NPlus6)           return 6;
   else if(species == NPlus7)           return 7;
   else if(species == AlPlus4)           return 4;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

int whatFunctionMode(char *str)
{
   if(strstr(str,"Constant")) 		return Constant;
   else if(strstr(str,"Gaussian"))   	return Gaussian;
   else if(strstr(str,"Polynomial"))   	return Polynomial;
   else return 0;
}

void random1D_sobol(double *x,gsl_qrng *q)
{
   double v[1];

   gsl_qrng_get(q,v);
   *x=v[0];
}

void random2D_sobol(double *x,double *y,gsl_qrng *q)
{
   double v[2];

   gsl_qrng_get(q,v);
   *x=v[0];
   *y=v[1];
}

void random3D_sobol(double *x,double *y,double *z,gsl_qrng *q)
{
   double v[3];

   gsl_qrng_get(q,v);
   *x=v[0];
   *y=v[1];
   *z=v[2];
}

void ionizationSetup(LoadList *LL,int species)
{

  switch (species)  {
    case Electron :
    case Positron :
      LL->levels=0;
    case HPlus0 :
    case HPlus1 :
      LL->levels=1;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.0;
      if(species==HPlus1) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case HePlus0 :
    case HePlus1 :
    case HePlus2 :
      LL->levels=1;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.808;
      LL->ionEnergy[1]=4.002;
      if(species==HePlus2) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case CPlus0 :
    case CPlus1 :
    case CPlus2 :
    case CPlus3 :
    case CPlus4 :
    case CPlus5 :
    case CPlus6 :
      LL->levels=6;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=0.828;
      LL->ionEnergy[1]=1.793;
      LL->ionEnergy[2]=3.522;
      LL->ionEnergy[3]=4.743;
      LL->ionEnergy[4]=28.833;
      LL->ionEnergy[5]=36.033;
      if(species==CPlus6) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case NPlus0 :
    case NPlus1 :
    case NPlus2 :
    case NPlus3 :
    case NPlus4 :
    case NPlus5 :
    case NPlus6 :
    case NPlus7 :
      LL->levels=7;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.069;
      LL->ionEnergy[1]=2.177;
      LL->ionEnergy[2]=3.489;
      LL->ionEnergy[3]=5.697;
      LL->ionEnergy[4]=7.199;
      LL->ionEnergy[5]=40.598;
      LL->ionEnergy[6]=49.053;
      if(species==NPlus7) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
  }

}

int whatLaserMode(char *str)
{
  if(strstr(str,"Gaussian"))       return Gaussian;
  else if(strstr(str,"SSTF"))      return SSTF;
  else           return 0;
}

int whatLoadingMethod(char *str)
{
  if(strstr(str,"Boundary"))       return Boundary;
  else if(strstr(str,"Shot"))      return Shot;
  else           return 0;
}


