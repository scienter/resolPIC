#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"

double dampingU(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x<=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

double dampingD(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x>=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

void absorb_L(Domain *D,double *leftr,double *leftd,double x,double lL,double LdL,double rr,double rd)
{
  double tmp;

  if(x<=lL && LdL>0)  {
    tmp=(lL-x)/LdL;
    *leftr=1.0-rr*rr*tmp*tmp;
    *leftd=1.0-rd*rd*tmp*tmp;
  } else  {  
    *leftr=1.0;
    *leftd=1.0;
  }
}

void absorb_R(Domain *D,double *rtr,double *rtd,double x,double rL,double LdR,double rr,double rd)
{
  double tmp;

  if(x>=rL && LdR>0)  {
    tmp=(x-rL)/LdR;
    *rtr=1.0-rr*rr*tmp*tmp;
    *rtd=1.0-rd*rd*tmp*tmp;
  } else  {  
    *rtr=1.0;
    *rtd=1.0;
  }
}

void absorb_RL3(Domain *D)
{
  double tmp,x,LdR,LdL,rr,rd,rtL,ltL;
  int i,istart,iend,minXSub;

  istart=D->istart;    iend=D->iend;
  minXSub=D->minXSub;

  rtL=(double)(D->nx-D->pmlCellRight);
  ltL=(double)(D->pmlCellLeft);
  LdR=D->pmlCellRight;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;
  rd=D->pmld;

  for(i=istart; i<iend; i++) {
    x=(i-istart)+minXSub;
    if(x>=rtL && LdR>0)  {
      tmp=(x-rtL)/LdR;
      D->rtr[i]=1.0-rr*rr*tmp*tmp;
      D->ltr[i]=1.0;
      D->rtd[i]=1.0-rd*rd*tmp*tmp;
      D->ltd[i]=1.0;
    } else if(x<=ltL && LdL>0)  {
      tmp=(ltL-x)/LdL;
      D->rtr[i]=1.0;
      D->ltr[i]=1.0-rr*rr*tmp*tmp;
      D->rtd[i]=1.0;
      D->ltd[i]=1.0-rd*rd*tmp*tmp;
    } else  {  
      D->rtr[i]=1.0;
      D->ltr[i]=1.0;
      D->rtd[i]=1.0;
      D->ltd[i]=1.0;
    }
  } 
}

void absorb_UD3(Domain *D)
{
  double tmp,y,LdU,LdD,rr,rd,upL,downL;
  int j,jstart,jend,minYSub;

  jstart=D->jstart;    jend=D->jend;
  minYSub=D->minYSub;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  downL=(double)(D->minYDomain+D->pmlCellDown);
  LdU=D->pmlCellUp;
  LdD=D->pmlCellDown;
  rr=D->pmlr;
  rd=D->pmld;

  for(j=jstart; j<jend; j++) {
    y=(j-jstart)+minYSub;
    if(y>=upL && LdU>0)  {
      tmp=(y-upL)/LdU;
      D->upr[j]=1.0-rr*rr*tmp*tmp;
      D->dnr[j]=1.0;
      D->upd[j]=1.0-rd*rd*tmp*tmp;
      D->dnd[j]=1.0;
    } else if(y<=downL && LdD>0)  {
      tmp=(downL-y)/LdD;
      D->upr[j]=1.0;
      D->dnr[j]=1.0-rr*rr*tmp*tmp;
      D->upd[j]=1.0;
      D->dnd[j]=1.0-rd*rd*tmp*tmp;
    } else  {  
      D->upr[j]=1.0;
      D->dnr[j]=1.0;
      D->upd[j]=1.0;
      D->dnd[j]=1.0;
    }
  } 
}

void absorb_FB3(Domain *D)
{
  double tmp,z,LdF,LdB,rr,rd,frontL,backL;
  int k,kstart,kend,minZSub;

  kstart=D->kstart;    kend=D->kend;
  minZSub=D->minZSub;

  frontL=(double)(D->nz+D->minZDomain-D->pmlCellFront);
  backL=(double)(D->minZDomain+D->pmlCellBack);
  LdF=D->pmlCellFront;
  LdB=D->pmlCellBack;
  rr=D->pmlr;
  rd=D->pmld;

  for(k=kstart; k<kend; k++) {
    z=(k-kstart)+minZSub;
    if(z>=frontL && LdF>0)  {
      tmp=(z-frontL)/LdF;
      D->frr[k]=1.0-rr*rr*tmp*tmp;
      D->bkr[k]=1.0;
      D->frd[k]=1.0-rd*rd*tmp*tmp;
      D->bkd[k]=1.0;
    } else if(z<=backL && LdB>0)  {
      tmp=(backL-z)/LdB;
      D->frr[k]=1.0;
      D->bkr[k]=1.0-rr*rr*tmp*tmp;
      D->frd[k]=1.0;
      D->bkd[k]=1.0-rd*rd*tmp*tmp;
    } else  {  
      D->frr[k]=1.0;
      D->bkr[k]=1.0;
      D->frd[k]=1.0;
      D->bkd[k]=1.0;
    }
  } 
}

void absorb_UD(Domain *D,double *upr,double *btr,double *upd,double *btd,double y,double upL,double bottomL,double LdU,double LdB,double rr,double rd)
{
  double tmp;

  if(y>=upL && LdU>0)  {
    tmp=(y-upL)/LdU;
    *upr=1.0-rr*rr*tmp*tmp;
    *btr=1.0;
    *upd=1.0-rd*rd*tmp*tmp;
    *btd=1.0;
  } else if(y<=bottomL && LdB>0)  {
    tmp=(bottomL-y)/LdB;
    *upr=1.0;
    *btr=1.0-rr*rr*tmp*tmp;
    *upd=1.0;
    *btd=1.0-rd*rd*tmp*tmp;
  } else  {
    *upr=1.0;
    *btr=1.0;
    *upd=1.0;
    *btd=1.0;
  }
}


void absorb_UD2(Domain *D,double *upr,double *btr,double *upd,double *btd,double y,double upL,double bottomL,double LdU,double LdB,double rr,double rd,int flagX)
{
  double tmp;

  if(y>=upL && LdU>0 && flagX==ON)  {
    tmp=(y-upL)/LdU;
    *upr=1.0-rr*rr*tmp*tmp;
    *btr=1.0;
    *upd=1.0-rd*rd*tmp*tmp;
    *btd=1.0;
  } else if(y<=bottomL && LdB>0 && flagX==ON)  {
    tmp=(bottomL-y)/LdB;
    *upr=1.0;
    *btr=1.0-rr*rr*tmp*tmp;
    *upd=1.0;
    *btd=1.0-rd*rd*tmp*tmp;
  } else  {  
    *upr=1.0;
    *btr=1.0;
    *upd=1.0;
    *btd=1.0;
  }
}
/*
void absorb3DC(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,kinit,kfinal,pmlCell;
  double prevSrC,nowSrC,prevPrC,nowPrC,tmp;
  double dx,dt,dy,dz,rr,y,z;
  double damp11,damp12,damp21,damp22,damp31,damp32,damp41,damp42;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  dx=D->dx;
  dy=D->dy;
  dz=D->dz;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
//  rd=D->pmld;

  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        }
      }
    }
    break;
  case UPFRONT :    
    jinit=jend-D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case UPBACK :    
    jinit=jend-D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWNBACK :
    jfinal=jstart+D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case DOWNFRONT :
    jfinal=jstart+D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*damp31*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*damp41*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp21*damp31*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp21*damp31*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp11*damp41*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp11*damp41*D->SlC[i][j][k];
        }
      }
    }
    break;
  case FRONT :    
    kinit=kend-D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kinit; k<kend; k++)
    {
      z=k-kinit+1;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp11*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp11*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp21*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp21*D->SlC[i][j][k];
        }
      }
    }
    break;
  case BACK :
    kfinal=kstart+D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kstart; k<kfinal; k++)
    {
      z=kfinal-k;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
          D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
          D->PrC[i][j][k]=damp11*D->PrC[i][j][k];
          D->PlC[i][j][k]=damp11*D->PlC[i][j][k];
          D->SrC[i][j][k]=damp21*D->SrC[i][j][k];
          D->SlC[i][j][k]=damp21*D->SlC[i][j][k];
        }
      }
    }
    break;
  }
}

void absorb3D(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,kinit,kfinal,pmlCell;
  double prevSr,nowSr,prevPr,nowPr;
  double dx,dt,dy,dz,rr,y,z;
  double damp11,damp12,damp21,damp22,damp31,damp32,damp41,damp42;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;
  dx=D->dx;
  dy=D->dy;
  dz=D->dz;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
//  rd=D->pmld;

  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case UPFRONT :    
    jinit=jend-D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case UPBACK :    
    jinit=jend-D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=jinit; j<jend+3; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Pr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }  
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    kinit=kstart+D->pmlCell;
    kfinal=kend-D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kstart; k<kend; k++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        }
      }
    } 
    break;
  case DOWNBACK :
    jfinal=jstart+D->pmlCell;
    kinit=kstart;
    kfinal=kstart+D->pmlCell;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=0; k<kfinal; k++)
      {
        z=kfinal-k;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }
    }
    break;
  case DOWNFRONT :
    jfinal=jstart+D->pmlCell;
    kinit=kend-D->pmlCell;
    kfinal=kend;
    for(j=0; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(k=kinit; k<kend+3; k++)
      {
        z=k-kinit+1;
        damp31=func(z,rr,pmlCell);
        damp41=func(z-0.5,rr,pmlCell);
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*damp31*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*damp41*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp21*damp31*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp21*damp31*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp11*damp41*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp11*damp41*D->Sl[i][j][k];
        }
      }
    }
    break;
  case FRONT :    
    kinit=kend-D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kinit; k<kend; k++)
    {
      z=k-kinit+1;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp11*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp11*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp21*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp21*D->Sl[i][j][k];
        }
      }
    }
    break;
  case BACK :
    kfinal=kstart+D->pmlCell;
    jinit=jstart+D->pmlCell;
    jfinal=jend-D->pmlCell;
    for(k=kstart; k<kfinal; k++)
    {
      z=kfinal-k;
      damp11=func(z,rr,pmlCell);
      damp21=func(z-0.5,rr,pmlCell);
      for(j=jstart; j<jend; j++)
      {
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
          D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
          D->Pr[i][j][k]=damp11*D->Pr[i][j][k];
          D->Pl[i][j][k]=damp11*D->Pl[i][j][k];
          D->Sr[i][j][k]=damp21*D->Sr[i][j][k];
          D->Sl[i][j][k]=damp21*D->Sl[i][j][k];
        }
      }
    }
    break;
  }
}

void absorbC(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,pmlCell;
  double prevSrC,nowSrC,prevPrC,nowPrC;
  double dx,dt,dy,dz,rr,rd,y,damp11,damp12,damp21,damp22;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  dx=D->dx;
  dy=D->dy;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;
  rd=D->pmld;

  
  switch(position) {
  case UP :    
    jinit=jend-pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
        D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
        D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
        D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
        D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
        D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->ExC[i][j][k]=damp11*D->ExC[i][j][k];
        D->SrC[i][j][k]=damp11*D->SrC[i][j][k];
        D->SlC[i][j][k]=damp11*D->SlC[i][j][k];
        D->BxC[i][j][k]=damp21*D->BxC[i][j][k];
        D->PrC[i][j][k]=damp21*D->PrC[i][j][k];
        D->PlC[i][j][k]=damp21*D->PlC[i][j][k];
      }
    }
    break;
  }
}

void absorb(Domain *D,int position)
{
  int istart,iend,jstart,jend,kstart,kend;
  int i,j,k,nxSub,nySub,nzSub,jinit,jfinal,pmlCell;
  double prevSr,nowSr,prevPr,nowPr;
  double dx,dt,dy,dz,rr,y,damp11,damp12,damp21,damp22;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  dx=D->dx;
  dy=D->dy;
  dt=D->dt;

  pmlCell=D->pmlCell;
  rr=D->pmlr;

  k=0;
  switch(position) {
  case UP :    
    jinit=jend-D->pmlCell;
    for(j=jinit; j<jend; j++)
    {
      y=j-jinit+1;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
        D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
        D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
        D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
        D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
        D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
      }
    }
    break;
  case DOWN :
    jfinal=jstart+D->pmlCell;
    for(j=jstart; j<jfinal; j++)
    {
      y=jfinal-j;
      damp11=func(y,rr,pmlCell);
      damp21=func(y-0.5,rr,pmlCell);
      for(i=istart; i<iend; i++)
      {
        D->Ex[i][j][k]=damp11*D->Ex[i][j][k];
        D->Sr[i][j][k]=damp11*D->Sr[i][j][k];
        D->Sl[i][j][k]=damp11*D->Sl[i][j][k];
        D->Bx[i][j][k]=damp21*D->Bx[i][j][k];
        D->Pr[i][j][k]=damp21*D->Pr[i][j][k];
        D->Pl[i][j][k]=damp21*D->Pl[i][j][k];
      }
    }
    break;
  }

}
*/
