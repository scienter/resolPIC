#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include <mpi.h>

void movingDomain(Domain *D,int iteration)
{
  int shiftDuration;
  void movingDomain_Split();
  void movingDomain_Pukhov();
  void MPI_Transfer3F_Xminus();
  void MPI_Transfer3F_Xplus();
  void MPI_Transfer6F_Xminus();
  void MPI_Transfer6F_Xplus();

  switch((D->fieldType-1)*3+D->dimension)  {
  case ((Split-1)*3+1) :
  case ((Split-1)*3+2) :
  case ((Split-1)*3+3) :
    movingDomain_Split(D,iteration);
    break;
  case ((Yee-1)*3+1) :
  case ((Yee-1)*3+2) :
  case ((Yee-1)*3+3) :
  case ((Pukhov-1)*3+1) :
  case ((Pukhov-1)*3+2) :
  case ((Pukhov-1)*3+3) :
    movingDomain_Pukhov(D);
    break;

  default:
    printf("In movingDomain.c, what fieldType?(%d) and what dimension?(%d)\n",D->fieldType,D->dimension);
  }
}

void movingDomain_Split(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s;
    double x;
    ptclList *p,*tmp,*New;      
    Particle ***particle;
    particle=D->particle;
    LoadList *LL;
    int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   
    istart=0;    iend=D->iend+2;
    jstart=0;    jend=1;
    kstart=0;    kend=1;
    if(D->dimension>1) { jend=D->jend+3; } else;
    if(D->dimension>2) { kend=D->kend+3; } else;

    switch (D->dimension)   {
    case 1 :
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)    {
            D->Ex[i][j][k]=D->Ex[i+1][j][k];
            D->Bx[i][j][k]=D->Bx[i+1][j][k];
            D->Pr[i][j][k]=D->Pr[i+1][j][k];
            D->Pl[i][j][k]=D->Pl[i+1][j][k];
            D->Sr[i][j][k]=D->Sr[i+1][j][k];
            D->Sl[i][j][k]=D->Sl[i+1][j][k];
            D->Jx[i][j][k]=D->Jx[i+1][j][k];
            D->Jy[i][j][k]=D->Jy[i+1][j][k];
            D->Jz[i][j][k]=D->Jz[i+1][j][k];
          }     
      break;
    case 2 :
    case 3 :
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)    {
            D->Ex[i][j][k]=D->Ex[i+1][j][k];
            D->Bx[i][j][k]=D->Bx[i+1][j][k];
            D->Pr[i][j][k]=D->Pr[i+1][j][k];
            D->Pl[i][j][k]=D->Pl[i+1][j][k];
            D->Sr[i][j][k]=D->Sr[i+1][j][k];
            D->Sl[i][j][k]=D->Sl[i+1][j][k];
            D->Jx[i][j][k]=D->Jx[i+1][j][k];
            D->Jy[i][j][k]=D->Jy[i+1][j][k];
            D->Jz[i][j][k]=D->Jz[i+1][j][k];
            D->ExC[i][j][k]=D->ExC[i+1][j][k];
            D->BxC[i][j][k]=D->BxC[i+1][j][k];
            D->PrC[i][j][k]=D->PrC[i+1][j][k];
            D->PlC[i][j][k]=D->PlC[i+1][j][k];
            D->SrC[i][j][k]=D->SrC[i+1][j][k];
            D->SlC[i][j][k]=D->SlC[i+1][j][k];
            D->JxOld[i][j][k]=D->JxOld[i+1][j][k];
            D->JyOld[i][j][k]=D->JyOld[i+1][j][k];
            D->JzOld[i][j][k]=D->JzOld[i+1][j][k];
          }     
      break;
    }
    istart=D->istart; iend=D->iend;
    if(D->dimension>1) { jstart=D->jstart; jend=D->jend; } else;
    if(D->dimension>2) { kstart=D->kstart; kend=D->kend; } else;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)    
          for(s=0; s<D->nSpecies; s++)    {
            p=particle[i][j][k].head[s]->pt;
            while(p)     {
              p->x-=1.0;  
              p->oldX-=1.0;
              p=p->next;
            }
          } 

//    for(i=istart-1; i<=iend; i++)
//      for(j=0; j<jend+3; j++)
//      {
//        D->JxBoost[i][j][k]=D->JxBoost[i+1][j][k];
//        D->JyBoost[i][j][k]=D->JyBoost[i+1][j][k];
//        D->JzBoost[i][j][k]=D->JzBoost[i+1][j][k];
//      }     
//    i=iend+1;
//    for(j=0; j<jend+3; j++)
//    {
//      D->JxBoost[i][j][k]=0.0;
//      D->JyBoost[i][j][k]=0.0;
//      D->JzBoost[i][j][k]=0.0;
//    }     

    D->minXSub+=1;
    D->maxXSub+=1;
    D->minXDomain+=1;
    D->maxXDomain+=1;
    LL=D->loadList;
    while(LL->next)      {
      LL->createGuard+=D->dx;
      LL->deleteGuard+=D->dx;
      LL=LL->next;
    }
}

void movingDomain_Pukhov(Domain *D)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,shiftDuration;
    double x;
    ptclList *p,*tmp,*New;      
    Particle ***particle;
    particle=D->particle;
    LoadList *LL;

    istart=0;    iend=D->iend+2;
    jstart=0;    jend=1;
    kstart=0;    kend=1;
    if(D->dimension>1) { jend=D->jend+3; } else;
    if(D->dimension>2) { kend=D->kend+3; } else;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)          {
          D->Ex[i][j][k]=D->Ex[i+1][j][k];
          D->Ey[i][j][k]=D->Ey[i+1][j][k];
          D->Ez[i][j][k]=D->Ez[i+1][j][k];
          D->Bx[i][j][k]=D->Bx[i+1][j][k];
          D->By[i][j][k]=D->By[i+1][j][k];
          D->Bz[i][j][k]=D->Bz[i+1][j][k];
          D->Jx[i][j][k]=D->Jx[i+1][j][k];
          D->Jy[i][j][k]=D->Jy[i+1][j][k];
          D->Jz[i][j][k]=D->Jz[i+1][j][k];
        }    

    istart=D->istart; iend=D->iend;
    if(D->dimension>1) { jstart=D->jstart; jend=D->jend; } else;
    if(D->dimension>2) { kstart=D->kstart; kend=D->kend; } else;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)          
          for(s=0; s<D->nSpecies; s++)      {
            p=particle[i][j][k].head[s]->pt;
            while(p)  {
              p->x-=1.0;  
              p->oldX-=1.0;
	      if(p->y>=2.0 || p->y<=-2.0) printf("at moving, y=%g\n",p->y);
              p=p->next;
            }
          }	

    D->minXSub+=1;
    D->maxXSub+=1;
    D->minXDomain+=1;
    D->maxXDomain+=1;
    LL=D->loadList;
    while(LL->next)      {
      LL->createGuard+=D->dx;
      LL->deleteGuard+=D->dx;
      LL=LL->next;
    }
}  
