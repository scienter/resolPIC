#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"

void absorb_RL3(Domain *D);
void absorb_UD3(Domain *D);
void absorb_FB3(Domain *D);
void solveR1D_Split(Domain *D);
void solveL1D_Split(Domain *D);
void solve2DC_Split(Domain *D,int itertation);
void solve2D_Split(Domain *D,int iteration);
void solve3DC_Split(Domain *D,int itertation);
void solve3D_Split(Domain *D,int iteration);
void Bsolve2D_Pukhov(Domain *D,int iteration);
void Esolve2D_Pukhov(Domain *D,int iteration);
void Bsolve3D_Pukhov(Domain *D,int iteration);
void Esolve3D_Pukhov(Domain *D,int iteration);
void Bsolve3D_Yee(Domain *D,int iteration);
void Esolve3D_Yee(Domain *D,int iteration);
void Bsolve2D_Yee(Domain *D,int iteration);
void Esolve2D_Yee(Domain *D,int iteration);
void Bsolve1D_Yee_Pukhov(Domain *D);
void Esolve1D_Yee_Pukhov(Domain *D);
void MPI_Transfer3F_Period_Yminus(Domain *D,double ***f1,double ***f2,double ***f3,double sign1,double sign2,double sign3,int nx,int nz,int share);
void MPI_Transfer3F_Period_Yplus(Domain *D,double ***f1,double ***f2,double ***f3,double sign1,double sign2,double sign3,int nx,int nz,int share);
void Bfield_Treatment(Domain *D,int iteration);

void fieldSolve(Domain D,double t,int iteration)
{
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;
  void MPI_Transfer4F_Xminus();
  void MPI_Transfer4F_Xplus();
  void MPI_Transfer3F_Xminus();
  void MPI_Transfer3F_Xplus();
  void MPI_Transfer6F_Xminus();
  void MPI_Transfer6F_Xplus();
  void MPI_Transfer6F_Period_Yminus();
  void MPI_Transfer6F_Period_Yplus();
  void MPI_Transfer3F_Yminus();
  void MPI_Transfer3F_Yplus();
  void MPI_Transfer6F_Yminus();
  void MPI_Transfer6F_Yplus();
  void MPI_Transfer3F_Zminus();
  void MPI_Transfer3F_Zplus();
  void MPI_Transfer6F_Zminus();
  void MPI_Transfer6F_Zplus();

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch((D.fieldType-1)*3+D.dimension) {

  //1D field
  case (Split-1)*3+1:
    if(D.boostOn==OFF)    {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }
    solveR1D_Split(&D);
    solveL1D_Split(&D);
    if(D.L>1) {
      MPI_Transfer6F_Xplus(&D,D.Ex,D.Pr,D.Sr,D.Bx,D.Pl,D.Sl,1,1,3);
      MPI_Transfer6F_Xminus(&D,D.Ex,D.Pr,D.Sr,D.Bx,D.Pl,D.Sl,1,1,3);
    } else ;
    break;

  //split 2D
  case (Split-1)*3+2:
    if(D.boostOn==OFF)    {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }

    solve2DC_Split(&D,iteration);
    if(D.L>1)  {
      MPI_Transfer4F_Xminus(&D,D.ExC,D.BxC,D.PlC,D.SlC,D.nySub+5,1,3);
      MPI_Transfer4F_Xplus(&D,D.ExC,D.BxC,D.PrC,D.SrC,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.SrC,D.SlC,D.ExC,D.nxSub+5,1,3);
      MPI_Transfer3F_Yplus(&D,D.PrC,D.PlC,D.BxC,D.nxSub+5,1,3);
    } else	;

    solve2D_Split(&D,iteration);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,1,3);
      MPI_Transfer6F_Xplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,1,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
      MPI_Transfer6F_Yplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
    } else	;
//    if(D.period==ON)  {
//      MPI_Transfer6F_Period_Yminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
//      MPI_Transfer6F_Period_Yplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub+5,1,3);
//    } else	;

    break;

  //3D
  case (Split-1)*3+3:
    solve3DC_Split(&D,iteration);
    if(D.N>1)  {
      MPI_Transfer3F_Zminus(&D,D.PrC,D.PlC,D.ExC,D.nxSub,D.nySub,3);
      MPI_Transfer3F_Zplus(&D,D.SrC,D.SlC,D.BxC,D.nxSub,D.nySub,3);
    } else  ;  
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.SrC,D.SlC,D.ExC,D.nxSub,D.nzSub+5,3);
      MPI_Transfer3F_Yplus(&D,D.PrC,D.PlC,D.BxC,D.nxSub,D.nzSub+5,3);
    } else  ;
    if(D.L>1)  {
      MPI_Transfer4F_Xminus(&D,D.ExC,D.BxC,D.PlC,D.SlC,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer4F_Xplus(&D,D.ExC,D.BxC,D.PrC,D.SrC,D.nySub+5,D.nzSub+5,3);
    } else  ;

    solve3D_Split(&D,iteration);

    if(D.boostOn==OFF)    {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
        L=L->next;
      }
	 } else ;

    if(D.N>1)  {
      MPI_Transfer6F_Zminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub,D.nySub,3);
      MPI_Transfer6F_Zplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub,D.nySub,3);
    } else  ;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub,D.nzSub+5,3);
      MPI_Transfer6F_Yplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nxSub,D.nzSub+5,3);
    } else  ;
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer6F_Xplus(&D,D.Ex,D.Pr,D.Pl,D.Bx,D.Sr,D.Sl,D.nySub+5,D.nzSub+5,3);
    } else  ;
    break;

  //Yee, Pukhov 1D
  case (Yee-1)*3+1:
  case (Pukhov-1)*3+1:
    Esolve1D_Yee_Pukhov(&D);
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,1,1,3);
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,1,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)       {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    }  else ;

    Bsolve1D_Yee_Pukhov(&D);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,1,1,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,1,1,3);
    } else	;
    break;

  //Yee 2D
  case (Yee-1)*3+2:
    Esolve2D_Yee(&D,iteration);
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,1,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,1,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    }
    if(D.periodY==ON && D.M>1)  {
      MPI_Transfer3F_Period_Yminus(&D,D.Ex,D.Ey,D.Ez,0,1,1,D.nxSub+5,1,3);
      MPI_Transfer3F_Period_Yplus(&D,D.Ex,D.Ey,D.Ez,0,1,1,D.nxSub+5,1,3);
    } else ;

    Bsolve2D_Yee(&D,iteration);
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,1,3);
      MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,1,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
    } else	;
    if(D.periodY==ON && D.M>1)  {
      MPI_Transfer3F_Period_Yminus(&D,D.Bx,D.By,D.Bz,0,1,1,D.nxSub+5,1,3);
      MPI_Transfer3F_Period_Yplus(&D,D.Bx,D.By,D.Bz,0,1,1,D.nxSub+5,1,3);
    } else ;
    
    break;

  //Yee 3D
  case (Yee-1)*3+3:
    Esolve3D_Yee(&D,iteration);
    if(D.N>1)  {
      MPI_Transfer3F_Zminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nySub,3);
      MPI_Transfer3F_Zplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nySub,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nzSub+5,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nzSub+5,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,D.nzSub+5,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    }

    Bsolve3D_Yee(&D,iteration);
    if(D.N>1)  {
      MPI_Transfer6F_Zminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nySub,3);
      MPI_Transfer6F_Zplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nySub,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nzSub+5,3);
      MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nzSub+5,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,D.nzSub+5,3);
    } else	;
    break;

  //Pukhov 2D
  case (Pukhov-1)*3+2:
    Esolve2D_Pukhov(&D,iteration);
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,1,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,1,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)       {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }  else ;

    Bsolve2D_Pukhov(&D,iteration);
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,1,3);
      MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,1,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,1,3);
    } else	;
    break;

  //Pukhov 3D
  case (Pukhov-1)*3+3:
    Esolve3D_Pukhov(&D,iteration);
    if(D.N>1)  {
      MPI_Transfer3F_Zminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nySub,3);
      MPI_Transfer3F_Zplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nySub,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nzSub+5,3);
      MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub,D.nzSub+5,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer3F_Xplus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer3F_Xminus(&D,D.Ex,D.Ey,D.Ez,D.nySub+5,D.nzSub+5,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)       {
      L=D.laserList;
      while(L->next)  {
        if(L->loadMethod==Boundary) loadLaser(&D,L,t); else ;
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
        L=L->next;
      }
    }  else ;

    Bsolve3D_Pukhov(&D,iteration);
    if(D.N>1)  {
      MPI_Transfer6F_Zminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nySub,3);
      MPI_Transfer6F_Zplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nySub,3);
    } else	;
    if(D.M>1)  {
      MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nzSub+5,3);
      MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub,D.nzSub+5,3);
    } else	;
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,D.nzSub+5,3);
      MPI_Transfer6F_Xplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nySub+5,D.nzSub+5,3);
    } else	;
    break;

  default:
    printf("what fieldType? and what dimension?\n");
  }
}


void Esolve3D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,rankY,rankZ;
    double dtOverdx,dtOverdy,dtOverdz,dbPiDt;
    double oldEx,oldEy,oldEz;
    double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
    double tmp,tmpr,tmpd;

    int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    istart=D->istart; iend=D->iend;
    jstart=D->jstart; jend=D->jend;
    kstart=D->kstart; kend=D->kend;

    for(k=kstart; k<kend; k++) {
      D->frr[k]=1.0;
      D->frd[k]=1.0;
      D->bkr[k]=1.0;
      D->bkd[k]=1.0;
    }
    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy; 
	 dtOverdx=D->dt/D->dx; 
	 dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;	 

    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	 

          oldEx=D->Ex[i][j][k];
          oldEy=D->Ey[i][j][k];
          oldEz=D->Ez[i][j][k];

          tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dtOverdz*(D->By[i][j][k]-D->By[i][j][k-1])-dbPiDt*D->Jx[i][j][k];
          D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

          tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])+dtOverdz*(D->Bx[i][j][k]-D->Bx[i][j][k-1])-dbPiDt*D->Jy[i][j][k];
          D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

          tmp=dtOverdx*(D->By[i][j][k]-D->By[i-1][j][k])-dtOverdy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-dbPiDt*D->Jz[i][j][k];
          D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
        }
      }
    }

    rankY=(myrank%(D->M*D->N))%D->M;
    rankZ=(myrank%(D->M*D->N))/D->M;
	 if(rankY==0) {
		j=jstart;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
	 if(rankZ==0) {
		k=kstart;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
}

void Bsolve3D_Yee(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
	double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz;
	double oldBx,oldBy,oldBz;
	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
   double tmp,tmpr,tmpd;

   int nTasks,myrank;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

    for(k=kstart; k<kend; k++) {
      D->frr[k]=1.0;
      D->frd[k]=1.0;
      D->bkr[k]=1.0;
      D->bkd[k]=1.0;
    }
    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
      D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy; 
	 dtOverdx=D->dt/D->dx; 
	 dtOverdz=D->dt/D->dz;

    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	 

          oldBx=D->Bx[i][j][k];
          oldBy=D->By[i][j][k];
          oldBz=D->Bz[i][j][k];

          tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k])+dtOverdz*(D->Ey[i][j][k+1]-D->Ey[i][j][k]);
          D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

          tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k])-dtOverdz*(D->Ex[i][j][k+1]-D->Ex[i][j][k]);
          D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

          tmp=dtOverdy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])-dtOverdx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
          D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

          D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
          D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
          D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
        }
      }
    }
}

void solveR1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;
  nxSub=D->nxSub;

  j=k=0;
  for(i=iend-1; i>=istart; i--)
  {
    D->Ex[i][j][k]+=-2.0*M_PI*dt*D->Jx[i][j][k];
    D->Pr[i][j][k]=D->Pr[i-1][j][k]-M_PI*dt*D->Jy[i][j][k];
    D->Sr[i][j][k]=D->Sr[i-1][j][k]-M_PI*dt*D->Jz[i][j][k];
  }  
}

void solveL1D_Split(Domain *D)
{
  int i,j,k,istart,iend,nxSub;
  double dx,dt;

  dx=D->dx;          dt=D->dt;
  istart=D->istart;  iend=D->iend;
  nxSub=D->nxSub;

  j=k=0;
  for(i=istart; i<iend; i++)
  {
    D->Pl[i][j][k]=D->Pl[i+1][j][k]-M_PI*dt*D->Jy[i][j][k];
    D->Sl[i][j][k]=D->Sl[i+1][j][k]-M_PI*dt*D->Jz[i][j][k];
  }  
}

void Bsolve1D_Yee_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,oldBx,oldBy,oldBz;

    dx=D->dx;    dy=D->dy;    dt=D->dt; 
    nxSub=D->nxSub;    nySub=D->nySub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          oldBx=D->Bx[i][j][k];
          oldBy=D->By[i][j][k];
          oldBz=D->Bz[i][j][k];
          D->Bx[i][j][k]=0.0;
          D->By[i][j][k]+=dt/dx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
          D->Bz[i][j][k]+=-dt/dx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
          D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
          D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
          D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
        }
}

void Esolve1D_Yee_Pukhov(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt;

    dx=D->dx;    dy=D->dy;    dz=D->dz;    dt=D->dt;
    nxSub=D->nxSub;    nySub=D->nySub;    nzSub=D->nzSub;
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          D->Ex[i][j][k]+=-2*M_PI*dt*D->Jx[i][j][k];
          D->Ey[i][j][k]+=-dt/dx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-2*M_PI*dt*D->Jy[i][j][k];
          D->Ez[i][j][k]+=dt/dx*(D->By[i][j][k]-D->By[i-1][j][k])-2*M_PI*dt*D->Jz[i][j][k];
        }
}

void Bsolve3D_Pukhov(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
	double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz;
	double oldBx,oldBy,oldBz;
	double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
	double tmp,tmpr,tmpd;

   int nTasks,myrank;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
     D->ltd[i]=1.0;
   }

   ay=0.125*D->dx/D->dy;
   az=0.125*D->dx/D->dz;
   ax=ay+az;
   bx=1.0-2.0*ax;
   by=1.0-2.0*ay;
   bz=1.0-2.0*az;
   dtOverdx=D->dt/D->dx;
   dtOverdy=D->dt/D->dy;
   dtOverdz=D->dt/D->dz;

   for(i=istart; i<iend; i++)
   {
     rtr=D->rtr[i]; ltr=D->ltr[i];
     rtd=D->rtd[i]; ltd=D->ltd[i];
     for(j=jstart; j<jend; j++)
     {
       upr=D->upr[j]; dnr=D->dnr[j];
       upd=D->upd[j]; dnd=D->dnd[j];
       for(k=kstart; k<kend; k++)
       {
         tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
         tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

         oldBx=D->Bx[i][j][k];
         oldBy=D->By[i][j][k];
         oldBz=D->Bz[i][j][k];

         tmp=-dtOverdy*(bz*(D->Ez[i][j+1][k]-D->Ez[i][j][k])+az*(D->Ez[i][j+1][k+1]+D->Ez[i][j+1][k-1]-D->Ez[i][j][k+1]-D->Ez[i][j][k-1]))+dtOverdz*(by*(D->Ey[i][j][k+1]-D->Ey[i][j][k])+ay*(D->Ey[i][j+1][k+1]+D->Ey[i][j-1][k+1]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
         D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

         tmp=dtOverdx*(bz*(D->Ez[i+1][j][k]-D->Ez[i][j][k])+az*(D->Ez[i+1][j][k+1]+D->Ez[i+1][j][k-1]-D->Ez[i][j][k+1]-D->Ez[i][j][k-1]))-dtOverdz*(bx*(D->Ex[i][j][k+1]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]));
         D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

         tmp=dtOverdy*(bx*(D->Ex[i][j+1][k]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]))-dtOverdx*(by*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+ay*(D->Ey[i+1][j+1][k]+D->Ey[i+1][j-1][k]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
         D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

         D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
         D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
         D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
       }
     }
   }
}

void Esolve3D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    double ax,ay,az,bx,by,bz,dtOverdx,dtOverdy,dtOverdz,dbPiDt;
    double oldEx,oldEy,oldEz;
    double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
    double tmp,tmpr,tmpd;

    int nTasks,myrank,rankY,rankZ;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
     D->ltd[i]=1.0;
   }

    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=0.125*D->dx/D->dy;
    az=0.125*D->dx/D->dz;
    ax=ay+az;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    bz=1.0-2.0*az;
    dtOverdx=D->dt/D->dx;
    dtOverdy=D->dt/D->dy;
    dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;

    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_FB3(D);
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=0.125*D->dx/D->dy;
    az=0.125*D->dx/D->dz;
    ax=ay+az;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    bz=1.0-2.0*az;
    dtOverdx=D->dt/D->dx;
    dtOverdy=D->dt/D->dy;
    dtOverdz=D->dt/D->dz;
    dbPiDt=2.0*M_PI*D->dt;

    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

          oldEx=D->Ex[i][j][k];
          oldEy=D->Ey[i][j][k];
          oldEz=D->Ez[i][j][k];

          tmp=dtOverdy*(bz*(D->Bz[i][j][k]-D->Bz[i][j-1][k])+az*(D->Bz[i][j][k+1]+D->Bz[i][j][k-1]-D->Bz[i][j-1][k+1]-D->Bz[i][j-1][k-1]))-dtOverdz*(by*(D->By[i][j][k]-D->By[i][j][k-1])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i][j+1][k-1]-D->By[i][j-1][k-1]))-dbPiDt*D->Jx[i][j][k];
          D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

          tmp=-dtOverdx*(bz*(D->Bz[i][j][k]-D->Bz[i-1][j][k])+az*(D->Bz[i][j][k+1]+D->Bz[i][j][k-1]-D->Bz[i-1][j][k+1]-D->Bz[i-1][j][k-1]))+dtOverdz*(bx*(D->Bx[i][j][k]-D->Bx[i][j][k-1])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j][k-1]-D->Bx[i-1][j][k-1]))-dbPiDt*D->Jy[i][j][k];
          D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

          tmp=dtOverdx*(by*(D->By[i][j][k]-D->By[i-1][j][k])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i-1][j+1][k]-D->By[i-1][j-1][k]))-dtOverdy*(bx*(D->Bx[i][j][k]-D->Bx[i][j-1][k])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i-1][j-1][k]))-dbPiDt*D->Jz[i][j][k];
          D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
        }
      }
    }

	 rankY=(myrank%(D->M*D->N))%D->M;
    rankZ=(myrank%(D->M*D->N))/D->M;
	 if(rankY==0) {
		j=jstart;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
	 if(rankZ==0) {
		k=kstart;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
}

void Bsolve2D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double oldBx,oldBy,oldBz,ax,ay,bx,by,x,y;
    double dtOverdx,dtOverdy,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    k=0;
    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    //Solving B field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];

        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

        tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

        tmp=dtOverdy*(bx*(D->Ex[i][j+1][k]-D->Ex[i][j][k])+ax*(D->Ex[i+1][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i-1][j][k]))-dtOverdx*(by*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+ay*(D->Ey[i+1][j+1][k]+D->Ey[i+1][j-1][k]-D->Ey[i][j+1][k]-D->Ey[i][j-1][k]));
        D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
    }
}



void Esolve2D_Pukhov(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double dx,dy,dz,dt,ax,ay,bx,by,x1,x2,x,y;
    double oldEx,oldEy,oldEz;
    double dtOverdx,dtOverdy,dbPiDt,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

	 int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
	 dbPiDt=2.0*M_PI*D->dt;

    k=0;
    //Solving E field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldEx=D->Ex[i][j][k];
        oldEy=D->Ey[i][j][k];
        oldEz=D->Ez[i][j][k];

        tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dbPiDt*D->Jx[i][j][k];
        D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

        tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-dbPiDt*D->Jy[i][j][k];
        D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

        tmp=dtOverdx*(by*(D->By[i][j][k]-D->By[i-1][j][k])+ay*(D->By[i][j+1][k]+D->By[i][j-1][k]-D->By[i-1][j+1][k]-D->By[i-1][j-1][k]))-dtOverdy*(bx*(D->Bx[i][j][k]-D->Bx[i][j-1][k])+ax*(D->Bx[i+1][j][k]+D->Bx[i-1][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i-1][j-1][k]))-dbPiDt*D->Jz[i][j][k];
        D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
      }

//		j=jstart-1;
//      D->Ey[i][j][k]-=dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k]);
    }

 	 if(myrank%D->M==0) {
     for(i=istart; i<iend; i++) D->Ex[i][jstart][k]=0.0;
   } else ;
}

void Bsolve2D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double oldBx,oldBy,oldBz,ax,ay,bx,by;
    double dtOverdx,dtOverdy,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    k=0;
    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    //Solving B field
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];

        tmp=-dtOverdy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->Bx[i][j][k]=tmpd*(oldBx+tmp*tmpr);

        tmp=dtOverdx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->By[i][j][k]=tmpd*(oldBy+tmp*tmpr);

        tmp=dtOverdy*(D->Ex[i][j+1][k]-D->Ex[i][j][k])-dtOverdx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
        D->Bz[i][j][k]=tmpd*(oldBz+tmp*tmpr);

        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
    }

}

void Esolve2D_Yee(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double dx,dy,dz,dt,ax,ay,bx,by;
    double oldEx,oldEy,oldEz;
    double dtOverdx,dtOverdy,dbPiDt,rtr,ltr,rtd,ltd;
    double tmp,tmpr,tmpd;
    int nTasks,myrank;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;
		D->upd[j]=1.0;
		D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    ay=ax=0.125*D->dx/D->dy;
    bx=1.0-2.0*ax;
    by=1.0-2.0*ay;
    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
	 dbPiDt=2.0*M_PI*D->dt;

    k=0;
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        oldEx=D->Ex[i][j][k];
        oldEy=D->Ey[i][j][k];
        oldEz=D->Ez[i][j][k];

        tmp=dtOverdy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-dbPiDt*D->Jx[i][j][k];
        D->Ex[i][j][k]=tmpd*(oldEx+tmp*tmpr);

        tmp=-dtOverdx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-dbPiDt*D->Jy[i][j][k];
        D->Ey[i][j][k]=tmpd*(oldEy+tmp*tmpr);

        tmp=dtOverdx*(D->By[i][j][k]-D->By[i-1][j][k])-dtOverdy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-dbPiDt*D->Jz[i][j][k];
        D->Ez[i][j][k]=tmpd*(oldEz+tmp*tmpr);
     }

   }

   if(myrank%D->M==0) {
     for(i=istart; i<iend; i++) D->Ex[i][jstart][k]=0.0;
   } else ;
}

void solve2DC_Split(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub;  
    double dtOverdx,dtOverdy,piDt,qtDtByDy,hfPiDt;
    double tmp,tmpd,tmpr,ltr,rtr,ltd,rtd;

    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    
    
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend; k=0;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy;
    dtOverdx=D->dt/D->dx;
	 piDt=M_PI*D->dt;
	 qtDtByDy=0.25*D->dt/D->dy;
	 hfPiDt=0.5*M_PI*D->dt;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=dtOverdy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])-piDt*(D->JxOld[i][j][k]+D->Jx[i][j][k]);
        D->ExC[i][j][k]=tmpd*(D->ExC[i][j][k]+tmp*tmpr);

        tmp=-dtOverdy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
        D->BxC[i][j][k]=tmpd*(D->BxC[i][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->Ex[i+1][j+1][k]+D->Ex[i][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i][j][k])-hfPiDt*(D->JyOld[i+1][j][k]+D->Jy[i+1][j][k]);
        D->PlC[i][j][k]=tmpd*(D->PlC[i+1][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i][j-1][k])-hfPiDt*(D->JzOld[i+1][j][k]+D->Jz[i+1][j][k]);
        D->SlC[i][j][k]=tmpd*(D->SlC[i+1][j][k]+tmp*tmpr);
      }	
    }

    // PrC,SrC
    for(i=iend-1; i>=istart; i--)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=qtDtByDy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-hfPiDt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
        D->PrC[i][j][k]=tmpd*(D->PrC[i-1][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-hfPiDt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
        D->SrC[i][j][k]=tmpd*(D->SrC[i-1][j][k]+tmp*tmpr);
      }	
	 }

    if(myrank%D->M==0) {
      for(i=istart; i<iend; i++) {
        D->ExC[i][jstart][k]=0.0;
//        D->PrC[i][jstart-1][k]=D->PrC[i][jstart][k];
      }
    } else ;
}

void solve2D_Split(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;  
    double dtOverdx,dtOverdy,piDt,qtDtByDy,dbPiDt;
    double tmp,tmpr,tmpd,ltr,ltd,rtr,rtd;

    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend; k=0;

    for(j=jstart; j<jend; j++) {
      D->upr[j]=1.0;                                                                                                       D->upd[j]=1.0;                                                                                                       D->dnr[j]=1.0;
      D->dnd[j]=1.0;
    }
    for(i=istart; i<iend; i++) {
      D->rtr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltr[i]=1.0;
      D->ltd[i]=1.0;
    }
    if(D->pmlOn==ON && iteration>D->pmlStart) {
      absorb_UD3(D);
      absorb_RL3(D);
    } else ;

    dtOverdy=D->dt/D->dy;
	 piDt=M_PI*D->dt;
	 qtDtByDy=0.25*D->dt/D->dy;
	 dbPiDt=2.0*M_PI*D->dt;

    // PrC,PlC,E1C,SrC,SlC,B1C
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=dtOverdy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])-dbPiDt*D->Jx[i][j][k];
        D->Ex[i][j][k]=tmpd*(D->Ex[i][j][k]+tmp*tmpr);

        tmp=-dtOverdy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]);
        D->Bx[i][j][k]=tmpd*(D->Bx[i][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->ExC[i+1][j+1][k]+D->ExC[i][j+1][k]-D->ExC[i+1][j][k]-D->ExC[i][j][k])-piDt*D->Jy[i+1][j][k];
        D->Pl[i][j][k]=tmpd*(D->Pl[i+1][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->BxC[i+1][j][k]+D->BxC[i][j][k]-D->BxC[i+1][j-1][k]-D->BxC[i][j-1][k])-piDt*D->Jz[i+1][j][k];
        D->Sl[i][j][k]=tmpd*(D->Sl[i+1][j][k]+tmp*tmpr);
      }
	 }

    // Pr,Sr
    for(i=iend-1; i>=istart; i--)
	 {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        tmpr=rtr*ltr*D->upr[j]*D->dnr[j];
		  tmpd=rtd*ltd*D->upd[j]*D->dnd[j];

        tmp=qtDtByDy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-piDt*D->Jy[i][j][k];
        D->Pr[i][j][k]=tmpd*(D->Pr[i-1][j][k]+tmp*tmpr);

        tmp=-qtDtByDy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-piDt*D->Jz[i][j][k];
        D->Sr[i][j][k]=tmpd*(D->Sr[i-1][j][k]+tmp*tmpr);
      }	
	 }

    if(myrank%D->M==0) {
      for(i=istart; i<iend; i++) {
        D->Ex[i][jstart][k]=0.0;
//        D->Pr[i][jstart-1][k]=D->Pr[i][jstart][k];
      }
    } else ;
}


void solve3DC_Split(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
   double dtOverdy,dtOverdz,qtDtByDy,qtDtByDz,hfPiDt,piDt;
   double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
   double tmp,tmpr,tmpd;

   int nTasks,myrank,rankY,rankZ;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;
 
   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
	  D->ltd[i]=1.0;
	}
	if(D->pmlOn==ON && iteration>D->pmlStart) {
	  absorb_FB3(D);
	  absorb_UD3(D);
	  absorb_RL3(D);
	} else ;

	dtOverdy=D->dt/D->dy;
	dtOverdz=D->dt/D->dz;
	piDt=M_PI*D->dt;
	qtDtByDy=0.25*D->dt/D->dy;
	qtDtByDz=0.25*D->dt/D->dz;
	hfPiDt=0.5*M_PI*D->dt;

    // PlC,ExC,SlC,BxC
    for(i=istart; i<iend; i++)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
          tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
          tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];

          tmp=dtOverdy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])+dtOverdz*(D->Sr[i][j][k]-D->Sr[i][j][k-1]-D->Sl[i][j][k]+D->Sl[i][j][k-1])-piDt*(D->JxOld[i][j][k]+D->Jx[i][j][k]);
          D->ExC[i][j][k]=tmpd*(D->ExC[i][j][k]+tmp*tmpr);

          tmp=dtOverdz*(D->Pr[i][j][k+1]-D->Pr[i][j][k]+D->Pl[i][j][k+1]-D->Pl[i][j][k])-dtOverdy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
          D->BxC[i][j][k]=tmpd*(D->BxC[i][j][k]+tmp*tmpr);

          tmp=qtDtByDz*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j][k-1]-D->Bx[i][j][k-1])-qtDtByDy*(D->Ex[i+1][j+1][k]+D->Ex[i][j+1][k]-D->Ex[i+1][j][k]-D->Ex[i][j][k])-hfPiDt*(D->JyOld[i+1][j][k]+D->Jy[i+1][j][k]);
          D->PlC[i][j][k]=tmpd*(D->PlC[i+1][j][k]+tmp*tmpr);

          tmp=-qtDtByDy*(D->Bx[i+1][j][k]+D->Bx[i][j][k]-D->Bx[i+1][j-1][k]-D->Bx[i][j-1][k])-qtDtByDz*(D->Ex[i+1][j][k+1]+D->Ex[i][j][k+1]-D->Ex[i+1][j][k]-D->Ex[i][j][k])-hfPiDt*(D->JzOld[i+1][j][k]+D->Jz[i+1][j][k]);
          D->SlC[i][j][k]=tmpd*(D->SlC[i+1][j][k]+tmp*tmpr);
        }
      }
    }

    // PrC,SrC
    for(i=iend-1; i>=istart; i--)
    {
      rtr=D->rtr[i]; ltr=D->ltr[i];
      rtd=D->rtd[i]; ltd=D->ltd[i];
      for(j=jstart; j<jend; j++)
      {
        upr=D->upr[j]; dnr=D->dnr[j];
        upd=D->upd[j]; dnd=D->dnd[j];
        for(k=kstart; k<kend; k++)
        {
           tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
           tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];
																																			 
           tmp=qtDtByDz*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j][k-1]-D->Bx[i-1][j][k-1])+qtDtByDy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-hfPiDt*(D->JyOld[i][j][k]+D->Jy[i][j][k]);
           D->PrC[i][j][k]=tmpd*(D->PrC[i-1][j][k]+tmp*tmpr);
																																										
           tmp=-qtDtByDy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])+qtDtByDz*(D->Ex[i][j][k+1]+D->Ex[i-1][j][k+1]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-hfPiDt*(D->JzOld[i][j][k]+D->Jz[i][j][k]);
           D->SrC[i][j][k]=tmpd*(D->SrC[i-1][j][k]+tmp*tmpr);
        }
      }
	 }	 

	 rankY=(myrank%(D->M*D->N))%D->M;
    rankZ=(myrank%(D->M*D->N))/D->M;
	 if(rankY==0) {
		j=jstart;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
	       D->ExC[i][j][k]=0.0;
	 } else ;
	 if(rankZ==0) {
		k=kstart;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
	       D->ExC[i][j][k]=0.0;
	 } else ;
}

void solve3D_Split(Domain *D,int iteration)
{
   int i,j,k,istart,iend,jstart,jend,kstart,kend;
   double dtOverdy,dtOverdz,qtDtByDy,qtDtByDz,dbPiDt,piDt;
   double rtr,rtd,ltr,ltd,upr,upd,dnr,dnd;
   double tmp,tmpr,tmpd;

   int nTasks,myrank,rankY,rankZ;
   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;    iend=D->iend;
   jstart=D->jstart;    jend=D->jend;
   kstart=D->kstart;    kend=D->kend;

   for(k=kstart; k<kend; k++) {
     D->frr[k]=1.0;
     D->frd[k]=1.0;
     D->bkr[k]=1.0;
     D->bkd[k]=1.0;
   }
   for(j=jstart; j<jend; j++) {
     D->upr[j]=1.0;
     D->upd[j]=1.0;
     D->dnr[j]=1.0;
     D->dnd[j]=1.0;
   }
   for(i=istart; i<iend; i++) {
     D->rtr[i]=1.0;
     D->rtd[i]=1.0;
     D->ltr[i]=1.0;
	  D->ltd[i]=1.0;
	}
	if(D->pmlOn==ON && iteration>D->pmlStart) {
	  absorb_FB3(D);
	  absorb_UD3(D);
	  absorb_RL3(D);
	} else ;


   dtOverdy=D->dt/D->dy;
   dtOverdz=D->dt/D->dz;
   piDt=M_PI*D->dt;
   qtDtByDy=0.25*D->dt/D->dy;
   qtDtByDz=0.25*D->dt/D->dz;
   dbPiDt=2.0*M_PI*D->dt;
								
   // Pl,Ex,Sl,Bx
   for(i=istart; i<iend; i++)
   {
     rtr=D->rtr[i]; ltr=D->ltr[i];
     rtd=D->rtd[i]; ltd=D->ltd[i];
     for(j=jstart; j<jend; j++)
     {
       upr=D->upr[j]; dnr=D->dnr[j];
       upd=D->upd[j]; dnd=D->dnd[j];
       for(k=kstart; k<kend; k++)
       {
         tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
         tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];	

         tmp=dtOverdy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])+dtOverdz*(D->SrC[i][j][k]-D->SrC[i][j][k-1]-D->SlC[i][j][k]+D->SlC[i][j][k-1])-dbPiDt*D->Jx[i][j][k];
         D->Ex[i][j][k]=tmpd*(D->Ex[i][j][k]+tmp*tmpr);
						  
         tmp=dtOverdz*(D->PrC[i][j][k+1]-D->PrC[i][j][k]+D->PlC[i][j][k+1]-D->PlC[i][j][k])-dtOverdy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]);
         D->Bx[i][j][k]=tmpd*(D->Bx[i][j][k]+tmp*tmpr);
													 
         tmp=qtDtByDz*(D->BxC[i+1][j][k]+D->BxC[i][j][k]-D->BxC[i+1][j][k-1]-D->BxC[i][j][k-1])-qtDtByDy*(D->ExC[i+1][j+1][k]+D->ExC[i][j+1][k]-D->ExC[i+1][j][k]-D->ExC[i][j][k])-piDt*D->Jy[i+1][j][k];
         D->Pl[i][j][k]=tmpd*(D->Pl[i+1][j][k]+tmp*tmpr);
																	
         tmp=-qtDtByDy*(D->BxC[i+1][j][k]+D->BxC[i][j][k]-D->BxC[i+1][j-1][k]-D->BxC[i][j-1][k])-qtDtByDz*(D->ExC[i+1][j][k+1]+D->ExC[i][j][k+1]-D->ExC[i+1][j][k]-D->ExC[i][j][k])-piDt*D->Jz[i+1][j][k];
         D->Sl[i][j][k]=tmpd*(D->Sl[i+1][j][k]+tmp*tmpr);
       }
     }
   }		
	
   //Pr,Sr
   for(i=iend-1; i>=istart; i--)
   {
     rtr=D->rtr[i]; ltr=D->ltr[i];
     rtd=D->rtd[i]; ltd=D->ltd[i];
     for(j=jstart; j<jend; j++)
     {
       upr=D->upr[j]; dnr=D->dnr[j];
       upd=D->upd[j]; dnd=D->dnd[j];
       for(k=kstart; k<kend; k++)
       {
         tmpr=rtr*ltr*upr*dnr*D->frr[k]*D->bkr[k];
         tmpd=rtd*ltd*upd*dnd*D->frd[k]*D->bkd[k];
																																			 
         tmp=qtDtByDz*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j][k-1]-D->BxC[i-1][j][k-1])+qtDtByDy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-piDt*D->Jy[i][j][k];
         D->Pr[i][j][k]=tmpd*(D->Pr[i-1][j][k]+tmp*tmpr);
			 			
         tmp=-qtDtByDy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])+qtDtByDz*(D->ExC[i][j][k+1]+D->ExC[i-1][j][k+1]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-piDt*D->Jz[i][j][k];
			D->Sr[i][j][k]=tmpd*(D->Sr[i-1][j][k]+tmp*tmpr);
       }
     }
   }	

	rankY=(myrank%(D->M*D->N))%D->M;
    rankZ=(myrank%(D->M*D->N))/D->M;
	 if(rankY==0) {
		j=jstart;
      for(i=istart; i<iend; i++)
        for(k=kstart; k<kend; k++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
	 if(rankZ==0) {
		k=kstart;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
	       D->Ex[i][j][k]=0.0;
	 } else ;
		
}
