#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "mesh.h"


int main(int argc, char *argv[])
{
   int s,i,j,k;
   Domain D;
   char fileName[100];
   FILE *out;
   int myrank, nTasks;
   MPI_Status status; 

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if(argc < 9) 
   {  
     if(myrank==0) {
       printf("mpirun -np N ./resolParticle dimension mode step L M N resolX resolY resolZ targetW1 targetW2\n"); 
       printf("mode 1 -> high resolution change for particle\n"); 
     } else ;
     MPI_Barrier(MPI_COMM_WORLD); 
     MPI_Finalize();     
     exit(0); 
   }
   D.targetW = (double *)malloc(10*sizeof(double ));


   D.dimension=atoi(argv[1]);
   D.mode=whatMode(argv[2]);
   D.step=atoi(argv[3]);
   D.L=atoi(argv[4]);
   D.M=atoi(argv[5]);
   D.N=atoi(argv[6]);
   D.resolX=atoi(argv[7]);
   D.resolY=atoi(argv[8]);
   D.resolZ=atoi(argv[9]);
   D.targetW[0]=atof(argv[10]);
   D.targetW[1]=atof(argv[11]);
   if(D.L*D.M*D.N!=nTasks)  {
     printf("check nTasks!. Now L=%d,M=%d,N=%d\n",D.L,D.M,D.N);
     exit(0);
   }

   boundary(&D);
   MPI_Barrier(MPI_COMM_WORLD); 
 
   if(D.mode==HIGH)
     restoreParticleHDF(&D,D.step,D.particle);
   else if(D.mode==LOW)
     restoreParticleHDF(&D,D.step+1,D.particle);
   else ;
   printf("restored Particles\n");

   printf("On reCreate partidles\n");
//   reCreateParticle(D); 
   changeParticlePosition(&D); 

   for(s=0; s<D.nSpecies; s++)
     saveParticle(&D,s,D.particle);

   saveDumpParticleHDF(&D);
   cleanParticle(&D);

   //field save
   if(D.mode==HIGH) reHighCreateField(D);
   else if(D.mode==LOW) reLowCreateField(D);
   else ;

   free(D.targetW);

    MPI_Finalize();

    return 0;
}

