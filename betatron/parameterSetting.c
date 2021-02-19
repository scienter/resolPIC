#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "parameter.h"
#include "constants.h"

void parameterSetting(Parameter *D,char *input)
{
//   int myrank, nTasks;
//   MPI_Status status;

//   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
//   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     
   int fail=0;
   double lambda,divisionLambda;
   int FindParameters();
   char str[100];

   //initially
   if(FindParameters("Parameter",1,"dimension",input,str)) D->dimension=atoi(str);
   else  {
      printf("in [Parameter], dimension=?  (1:1D, 2:2D, 3:3D)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"minThX",input,str)) D->minThX=atof(str);
   else  {
      printf("in [Parameter], minThX=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"maxThX",input,str)) D->maxThX=atof(str);
   else  {
      printf("in [Parameter], maxThX=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"numThX",input,str)) D->numThX=atoi(str);
   else  {
      printf("in [Parameter], numThX=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"minThY",input,str)) D->minThY=atof(str);
   else  {
      printf("in [Parameter], minThY=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"maxThY",input,str)) D->maxThY=atof(str);
   else  {
      printf("in [Parameter], maxThY=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"numThY",input,str)) D->numThY=atoi(str);
   else  {
      printf("in [Parameter], numThY=?  (rad)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"species",input,str)) D->species=atoi(str);
   else  {
      printf("in [Parameter], species=?  (0 or 1 or ...)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"lambda",input,str)) lambda=atof(str);
   else  {
      printf("in [Parameter], lambda=?  [m]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"divisionLambda",input,str)) divisionLambda=atof(str);
   else  {
      printf("in [Parameter], divisionLambda=?  [lambda/dx]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"super_particle",input,str)) D->superP=atof(str);
   else  {
      printf("in [Parameter], super_particle=?  (in PIC code)\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"detector_r",input,str)) D->det_R=atof(str);
   else  {
      printf("in [Parameter], detector_r=?  [m]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"detector_x",input,str)) D->det_x=atof(str);
   else  {
      printf("in [Parameter], detector_x=?  [m]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"detector_y",input,str)) D->det_y=atof(str);
   else  {
      printf("in [Parameter], detector_y=?  [m]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"min_energy",input,str)) D->minE=atof(str);
   else  {
      printf("in [Parameter], min_energy=?  [eV]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"max_energy",input,str)) D->maxE=atof(str);
   else  {
      printf("in [Parameter], max_energy=?  [eV]\n");
      fail=1;
   }
   if(FindParameters("Parameter",1,"energy_gap",input,str)) D->dE=atof(str);
   else  {
      printf("in [Parameter], energy_gap=?  [eV]\n");
      fail=1;
   }

   D->dt=lambda/divisionLambda/velocityC;
   if(fail==1)
      exit(0);

}

