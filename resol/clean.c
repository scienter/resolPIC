#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <math.h>
#include <mpi.h>

void deleteField(double ***field,int nx,int ny,int nz);

void cleanParticle(Domain *D)
{
   int i,j,k,s,nxSub,nySub,nzSub;
   ptclList *p,*tmp;

    //remove field
    nxSub=D->reNxSub+5; nySub=1; nzSub=1;
    if(D->dimension>1) nySub=D->reNySub+5; else;
    if(D->dimension>2) nzSub=D->reNzSub+5; else;

   switch (D->dimension) {
   case 1 :
     j=k=0;
     for(i=D->istart; i<D->iend; i++)  
     {
       for(s=0; s<D->nSpecies; s++)  {
         p=D->particle[i][j][k].head[s]->pt;
         while(p)    {
           tmp=p->next;
           D->particle[i][j][k].head[s]->pt=tmp;
           p->next=NULL;  free(p);
           p=D->particle[i][j][k].head[s]->pt;
         }
         free(D->particle[i][j][k].head[s]);
       }
       free(D->particle[i][j][k].head);
     }
     break;
   case 2 :
     k=0;
     for(i=D->istart; i<D->reIend; i++)  
       for(j=D->jstart; j<D->reJend; j++)  
       {
         for(s=0; s<D->nSpecies; s++)  {
           p=D->particle[i][j][k].head[s]->pt;
           while(p)    {
             tmp=p->next;
             D->particle[i][j][k].head[s]->pt=tmp;
             p->next=NULL;  free(p);
             p=D->particle[i][j][k].head[s]->pt;
           }
           free(D->particle[i][j][k].head[s]);
         }
         free(D->particle[i][j][k].head);
       }
     break;

   }
   for(i=0; i<nxSub; i++)  {
     for(j=0; j<nySub; j++)  {     
       free(D->particle[i][j]);
     }
     free(D->particle[i]);
   }
   free(D->particle);

}


void deleteField(double ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)  
   {
     for(j=0; j<ny; j++)  
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}

void deleteFieldInt(int ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)
   {
     for(j=0; j<ny; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


