#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void fieldIonization1D(Domain *D);
void fieldIonization2D(Domain *D);
void fieldIonization3D(Domain *D);
void RungeKutta(double *W,double *prob,int iter,double dt,int start,int end,int flag);

void fieldIonization(Domain D,int iteration)
{

  switch(D.dimension)  {
    case 1 :
      fieldIonization1D(&D);
      break;
    case 2 :
      fieldIonization2D(&D);
      break;
    case 3 :
      fieldIonization3D(&D);
      break;
  default :
    printf("In fieldIonization(), what dimension?\n");
  }

}

void fieldIonization3D(Domain *D)
{
  int ii,jj,kk,i,j,k,s,Z,i1,j1,k1,n;
  int *levels,*species,*ionFinal,initZ,cnt,ionLevel;
  int istart,iend,jstart,jend,kstart,kend,nSpecies;
  double Ea,Edc,wa,Efield,c2,n_eff,ionE,phase,beforeW,nowW,testProb,dt;
  double **ionEnergy,*rho,wx[2],WX[2],wy[2],WY[2],wz[2],WZ[2];
  double x,y,z,prob,unitEBydt,value,weight,x1,y1,z1;
  Particle ***particle;
  particle=D->particle;
  LoadList *LL;

  ptclList *p,*New;
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  dt=D->dt; nSpecies=D->nSpecies;

  Ea=eCharge*5.1e11/eMass/D->omega/velocityC;
  wa=4.13e16/D->omega;
  unitEBydt=13.59844/0.511*1.0e-6/dt;

  rho=(double *)malloc(nSpecies*sizeof(double ));
  levels=(int *)malloc(nSpecies*sizeof(int ));
  species=(int *)malloc(nSpecies*sizeof(int ));
  ionFinal=(int *)malloc(nSpecies*sizeof(int ));
  ionEnergy=(double **)malloc(nSpecies*sizeof(double *));

  LL=D->loadList;
  s=0;
  while(LL->next)
  {
    rho[s]=LL->density/LL->criticalDensity;
    levels[s]=LL->levels;
    species[s]=LL->species;
    ionFinal[s]=LL->ionFinal;
    //ionEnergy Z number : 0,1,2,...,Z-1
    ionEnergy[s]=(double *)malloc(LL->levels*sizeof(double ));
    
    for(i=0; i<LL->levels; i++) 
      ionEnergy[s][i]=LL->ionEnergy[i];
    LL=LL->next;
    s++;
  }

  //initialize J
  if(D->fieldType==Split)  {
    for(i=0; i<iend+3; i++)  
      for(j=0; j<jend+3; j++)  
        for(k=0; k<kend+3; k++)  {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  } else {
    for(i=0; i<iend+3; i++)  
      for(j=0; j<jend+3; j++)  
        for(k=0; k<kend+3; k++)  {
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  }


  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)  
      for(k=kstart; k<kend; k++)  
      {
        for(s=0; s<nSpecies; s++) {
          //when except cases
          if(species[s]==Electron || species[s]==Test || ionFinal[s]==ON) ;
          //when ionization process

          else {
            p=particle[i][j][k].head[s]->pt;

            while(p)  {
              Edc=sqrt(p->E1*p->E1+p->E2*p->E2+p->E3*p->E3);
              initZ=p->charge;
              weight=p->weight;

              if(initZ<levels[s] && Edc>0.0) {
                //test probability
                testProb=randomValue(1.0);

                Efield=Ea/Edc;
                Z=levels[s]-1; ionLevel=-1;               
                while(Z>=initZ) {
                  ionE=ionEnergy[s][Z];
                  n_eff=(Z+1)*sqrt(1.0/ionE);
                  c2=pow(5.43656/n_eff,2.0*n_eff)/n_eff;
                  phase=2.0*Efield*pow(ionE,1.5);
                  nowW=wa*c2*(0.5*ionE)*pow(phase,2.0*n_eff-1)*exp(-phase/3.0)*dt;
                  prob=1.0-exp(-nowW);
                  if(prob>testProb) { ionLevel=Z; break; }
                  else Z--;
                }

                if(ionLevel>=initZ) { 
                  value=ionEnergy[s][ionLevel]*unitEBydt*rho[s]/Edc/Edc*weight;
                  p->charge=ionLevel+1;
                  //creating electrons
                  cnt=p->charge-initZ;
                  for(n=0; n<cnt; n++) {
                    New = (ptclList *)malloc(sizeof(ptclList));
                    New->next = particle[i][j][k].head[s-1]->pt;
                    particle[i][j][k].head[s-1]->pt = New;
                    x=randomValue(1.0); y=randomValue(1.0); z=randomValue(1.0);
                    New->x = x; New->oldX=i+x;
                    New->y = y; New->oldY=j+y;
                    New->z = z; New->oldZ=k+z;
                    New->p1=p->p1; New->p2=p->p2; New->p3=p->p3;
                    New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
                    New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
                    New->weight=weight;
                    New->charge=-1;
                    New->index=p->index;
                    New->core=p->core;

                    x=p->x; y=p->y; z=p->z;
                    x1=(x+0.5)-(int)(x+0.5);
                    y1=(y+0.5)-(int)(y+0.5);
                    z1=(z+0.5)-(int)(z+0.5);
                    i1=i+(int)(x+0.5)-1;
                    j1=j+(int)(y+0.5)-1;
                    k1=k+(int)(z+0.5)-1;
 
                    wx[1]=x; wx[0]=1.0-x;
                    wy[1]=y; wy[0]=1.0-y;
                    wz[1]=z; wz[0]=1.0-z;
                    WX[1]=x1; WX[0]=1.0-x1;
                    WY[1]=y1; WY[0]=1.0-y1;
                    WZ[1]=z1; WZ[0]=1.0-z1;

                    for(ii=0; ii<2; ii++)
                      for(jj=0; jj<2; jj++) 
                        for(kk=0; kk<2; kk++) {
                          D->Jx[i1+ii][j+jj][k+kk]+=WX[ii]*wy[jj]*wz[kk]*p->E1*value;
                          D->Jy[i+ii][j1+jj][k+kk]+=wx[ii]*WY[jj]*wz[kk]*p->E2*value;
                          D->Jz[i+ii][j+jj][k1+kk]+=wx[ii]*wy[jj]*WZ[kk]*p->E3*value;
                        }
                  }  	//End for(n<cnt)
                } else ;

              } else ;  //End of if(initZ<levels[s])

              p=p->next;

            }  //End of while(p)

          }    //End of if(speces test) else;

        }      //End of for(s)
      }        //End of for(i,j,k)
          

  free(rho); free(levels); free(species); free(ionFinal);
  for(n=0; n<nSpecies; n++) free(ionEnergy[n]); free(ionEnergy);
}

/*
void fieldIonization2D(Domain *D,int iteration)
{
  int ii,jj,i,j,k,s,Z,i1,j1,n;
  int *levels,*species,*ionFinal,initZ,cnt,ionLevel;
  int istart,iend,jstart,jend,kstart,kend,nSpecies;
  double Ea,Edc,wa,Efield,c2,n_eff,ionE,phase,beforeW,nowW,testProb,dt;
  double **ionEnergy,*rho,wx[2],WX[2],wy[2],WY[2];
  double x,y,z,prob,unitEBydt,value,weight,x1,y1;
  Particle ***particle;
  particle=D->particle;
  LoadList *LL;

  ptclList *p,*New;
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  dt=D->dt; nSpecies=D->nSpecies;

  Ea=eCharge*5.1e11/eMass/D->omega/velocityC;
  wa=4.13e16/D->omega;
  unitEBydt=13.59844/0.511*1.0e-6/dt;

  rho=(double *)malloc(nSpecies*sizeof(double ));
  levels=(int *)malloc(nSpecies*sizeof(int ));
  species=(int *)malloc(nSpecies*sizeof(int ));
  ionFinal=(int *)malloc(nSpecies*sizeof(int ));
  ionEnergy=(double **)malloc(nSpecies*sizeof(double *));

  LL=D->loadList;
  s=0;
  while(LL->next)
  {
    rho[s]=LL->density/LL->criticalDensity;
    levels[s]=LL->levels;
    species[s]=LL->species;
    ionFinal[s]=LL->ionFinal;
    //ionEnergy Z number : 0,1,2,...,Z-1
    ionEnergy[s]=(double *)malloc(LL->levels*sizeof(double ));
    //prob for Z : 0,1,2,...,Z-1,Z
//    prob[s]=(double *)malloc((LL->levels+1)*sizeof(double ));
    
    for(i=0; i<LL->levels; i++) 
      ionEnergy[s][i]=LL->ionEnergy[i];
    LL=LL->next;
    s++;
  }

  //initialize J
  k=0;
  if(D->fieldType==Split)  {
    for(i=0; i<iend+3; i++)  
      for(j=0; j<jend+3; j++)  {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  } else {
    for(i=0; i<iend+3; i++)  
      for(j=0; j<jend+3; j++)  {
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  }


  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)  
      {
        for(s=0; s<nSpecies; s++) {
          //when except cases
          if(species[s]==Electron || species[s]==Test || ionFinal[s]==ON) ;
          //when ionization process

          else {
            p=particle[i][j][k].head[s]->pt;

            while(p)  {
              Edc=sqrt(p->E1*p->E1+p->E2*p->E2+p->E3*p->E3);
              initZ=p->charge;
              weight=p->weight;

              if(initZ<levels[s] && Edc>0.0) {
                //test probability
                testProb=randomValue(1.0);

                Efield=Ea/Edc;
                Z=levels[s]-1; ionLevel=-1;               
                while(Z>=initZ) {
                  ionE=ionEnergy[s][Z];
                  n_eff=(Z+1)*sqrt(1.0/ionE);
                  c2=pow(5.43656/n_eff,2.0*n_eff)/n_eff;
                  phase=2.0*Efield*pow(ionE,1.5);
                  nowW=wa*c2*(0.5*ionE)*pow(phase,2.0*n_eff-1)*exp(-phase/3.0)*dt;
                  prob=1.0-exp(-nowW);
                  if(prob>testProb) { ionLevel=Z; break; }
                  else Z--;
                }

                if(ionLevel>=initZ) { 
                  value=ionEnergy[s][ionLevel]*unitEBydt*rho[s]/Edc/Edc*weight;
                  p->charge=ionLevel+1;
                  //creating electrons
                  cnt=p->charge-initZ;
                  for(n=0; n<cnt; n++) {
                    New = (ptclList *)malloc(sizeof(ptclList));
                    New->next = particle[i][j][k].head[s-1]->pt;
                    particle[i][j][k].head[s-1]->pt = New;
                    x=randomValue(1.0); y=randomValue(1.0); z=0.0;
                    New->x = x; New->oldX=i+x;
                    New->y = y; New->oldY=j+y;
                    New->z = z; New->oldZ=z;
                    New->p1=p->p1; New->p2=p->p2; New->p3=p->p3;
                    New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
                    New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
                    New->weight=weight;
                    New->charge=-1;
                    New->index=p->index;
                    New->core=p->core;

                    x=p->x; y=p->y;
                    x1=(x+0.5)-(int)(x+0.5);
                    y1=(y+0.5)-(int)(y+0.5);
                    i1=i+(int)(x+0.5)-1;
                    j1=j+(int)(y+0.5)-1;
 
                    wx[1]=x; wx[0]=1.0-x;
                    wy[1]=y; wy[0]=1.0-y;
                    WX[1]=x1; WX[0]=1.0-x1;
                    WY[1]=y1; WY[0]=1.0-y1;
                    for(ii=0; ii<2; ii++)
                      for(jj=0; jj<2; jj++) {
                        D->Jx[i1+ii][j+jj][k]+=WX[ii]*wy[jj]*p->E1*value;
                        D->Jy[i+ii][j1+jj][k]+=wx[ii]*WY[jj]*p->E2*value;
                        D->Jz[i+ii][j+jj][k]+=wx[ii]*wy[jj]*p->E3*value;
                      }
                  }  	//End for(n<cnt)
                } else ;

              } else ;  //End of if(initZ<levels[s])

              p=p->next;

            }  //End of while(p)

          }    //End of if(speces test) else;

        }      //End of for(s)
      }        //End of for(i,j,k)
          

  free(rho); free(levels); free(species); free(ionFinal);
  for(n=0; n<nSpecies; n++) free(ionEnergy[n]); free(ionEnergy);
}
*/


//lala
void fieldIonization2D(Domain *D)
{
	
  int n,i,j,k,s,Z,i1,j1,ii,jj,flag;
  int initZ,cnt,ionLevel,index,sCnt,*sList;
  int istart,iend,jstart,jend,kstart,kend,nSpecies;
  double Ea,Edc,wa,Efield,c2,n_eff,l_eff,ionE,phase,dt;
  double **ionEnergy,**W,**prob,wx[2],wy[2],WX[2],WY[2];
  double x,y,z,unitEBydt,value,weight,x1,y1,testProb;
  Particle ***particle;
  particle=D->particle;
  LoadList *LL;

  ptclList *p,*New;
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  dt=D->dt; nSpecies=D->nSpecies;
 
  
  Ea=eCharge*5.1e11/eMass/D->omega/velocityC;
  wa=4.13e16/D->omega;

  unitEBydt=13.59844/0.511*1.0e-6/dt;

  double rho[nSpecies];
  double minA0[nSpecies];
  int levels[nSpecies],species[nSpecies],ionFinal[nSpecies];
  ionEnergy=(double **)malloc(nSpecies*sizeof(double *));
  W=(double **)malloc(nSpecies*sizeof(double *));
  prob=(double **)malloc(nSpecies*sizeof(double *));


  LL=D->loadList;
  s=0; sCnt=0;
  while(LL->next)
  {
    rho[s]=LL->density/LL->criticalDensity;
    levels[s]=LL->levels;
    species[s]=LL->species;
    ionFinal[s]=LL->ionFinal;
    minA0[s]=LL->givenMinA0;
    //ionEnergy Z number : 0,1,2,...,Z-1
    ionEnergy[s]=(double *)malloc(LL->levels*sizeof(double ));
    W[s]=(double *)malloc(LL->levels*sizeof(double ));
    prob[s]=(double *)malloc((LL->levels+1)*sizeof(double ));
    prob[s][LL->levels]=1.0;

    for(i=0; i<LL->levels; i++) 
      ionEnergy[s][i]=LL->ionEnergy[i];

    if(LL->species!=Electron && LL->ionFinal==OFF) sCnt++; else ;
    LL=LL->next;
    s++;
  }

  sList=(int *)malloc(sCnt*sizeof(int ));
  LL=D->loadList;
  s=0; index=0;
  while(LL->next)
  {
    if(LL->species!=Electron && LL->ionFinal==OFF) {
      sList[index]=s; index++;
    }  else ;
    LL=LL->next;
    s++;
  }


  //initialize J
  k=0;
  if(D->fieldType==Split)  {
    for(i=0; i<iend+3; i++)  
      for(j=0; j<jend+3; j++)  {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  } else {
    for(i=0; i<iend+3; i++) 
      for(j=0; j<jend+3; j++)  {
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  }
 
  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      {
        for(index=0; index<sCnt; index++) {
          s=sList[index];

            p=particle[i][j][k].head[s]->pt;

            while(p)  {
              Edc=sqrt(p->E1*p->E1+p->E2*p->E2+p->E3*p->E3);
              initZ=p->charge;
              weight=p->weight;
 
              if(initZ<levels[s] && Edc>minA0[s])  {
                testProb=randomValue(1.0);
  
                for(Z=initZ; Z<levels[s]; Z++) {
                  Efield=Ea/Edc;
                  ionE=ionEnergy[s][Z];
                  n_eff=(Z+1)*sqrt(1.0/ionE);
                  c2=pow(5.43656/n_eff,2.0*n_eff)/(n_eff);

                  phase=2.0*Efield*pow(ionE,1.5);
                  W[s][Z]=wa*c2*(0.5*ionE)*pow(phase,2.0*n_eff-1)*exp(-1.0*phase/3.0);
                }
                RungeKutta(&W[s][0],&prob[s][0],10,dt,initZ,levels[s],flag);
flag=0;
                Z=initZ; ionLevel=-1;
                while(Z<levels[s]) {
                  if(prob[s][Z]<testProb && testProb<prob[s][Z+1]) {
                    ionLevel=Z;
                    break;
                  }
                  Z++;
                }

                if(ionLevel>=initZ) { 
                  value=ionEnergy[s][ionLevel]*unitEBydt*rho[s]/Edc/Edc*weight;
                  p->charge=ionLevel+1;
                  //creating electrons
                  cnt=p->charge-initZ;
                  for(n=0; n<cnt; n++) {
                    New = (ptclList *)malloc(sizeof(ptclList));
                    New->next = particle[i][j][k].head[s-1]->pt;
                    particle[i][j][k].head[s-1]->pt = New;
                    x=randomValue(1.0); y=randomValue(1.0);  z=0.0;
                    New->x = x; New->oldX=i+x;
                    New->y = y; New->oldY=j+y;
                    New->z = z; New->oldZ=k+z;
                    New->p1=p->p1; New->p2=p->p2; New->p3=p->p3;
                    New->E1=p->E1; New->E2=p->E2; New->E3=p->E3;
                    New->B1=p->B1; New->B2=p->B2; New->B3=p->B3;
                    New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
                    New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
                    New->weight=weight;
                    New->charge=-1;
                    New->index=p->index;
                    New->core=p->core;

                    x=p->x; y=p->y;
                    x1=(x+0.5)-(int)(x+0.5);
                    y1=(y+0.5)-(int)(y+0.5);
                    i1=i+(int)(x+0.5)-1;
                    j1=j+(int)(y+0.5)-1;
 
                    wx[1]=x; wx[0]=1.0-x;
                    wy[1]=y; wy[0]=1.0-y;
                    WX[1]=x1; WX[0]=1.0-x1;
                    WY[1]=y1; WY[0]=1.0-y1;
                    for(ii=0; ii<2; ii++)
                      for(jj=0; jj<2; jj++) {
                        D->Jx[i1+ii][j+jj][k]+=WX[ii]*wy[jj]*p->E1*value;
                        D->Jy[i+ii][j1+jj][k]+=wx[ii]*WY[jj]*p->E2*value;
                        D->Jz[i+ii][j+jj][k]+=wx[ii]*wy[jj]*p->E3*value;
                      }
                  }
                } else ;	//End of electron generation 

              }  	//End of if(initZ<levels[s])

              p=p->next;

            }  //End of while(p)


        }      //End of for(s)
      }        //End of for(i,j,k)
         

    free(sList);
    for(n=0; n<nSpecies; n++) {
      free(ionEnergy[n]); free(W[n]); free(prob[n]);
    } free(ionEnergy); free(W); free(prob);

}



void fieldIonization1D(Domain *D)
{
  int n,i,j,k,s,Z,i1,flag;
  int initZ,cnt,ionLevel,index,sCnt,*sList;
  int istart,iend,jstart,jend,kstart,kend,nSpecies;
  double Ea,Edc,wa,Efield,c2,n_eff,l_eff,ionE,phase,dt;
  double **ionEnergy,**W,**prob;
  double x,y,z,unitEBydt,value,weight,x1,testProb;
  Particle ***particle;
  particle=D->particle;
  LoadList *LL;

  ptclList *p,*New;
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  kstart=D->kstart;   kend=D->kend;
  dt=D->dt; nSpecies=D->nSpecies;

  Ea=eCharge*5.1e11/eMass/D->omega/velocityC;
  wa=4.13e16/D->omega;

  unitEBydt=13.59844/0.511*1.0e-6/dt;

  double rho[nSpecies];
  double minA0[nSpecies];
  int levels[nSpecies],species[nSpecies],ionFinal[nSpecies];
  ionEnergy=(double **)malloc(nSpecies*sizeof(double *));
  W=(double **)malloc(nSpecies*sizeof(double *));
  prob=(double **)malloc(nSpecies*sizeof(double *));

  LL=D->loadList;
  s=0; sCnt=0;
  while(LL->next)
  {
    rho[s]=LL->density/LL->criticalDensity;
    levels[s]=LL->levels;
    species[s]=LL->species;
    ionFinal[s]=LL->ionFinal;
    minA0[s]=LL->givenMinA0;
    //ionEnergy Z number : 0,1,2,...,Z-1
    ionEnergy[s]=(double *)malloc(LL->levels*sizeof(double ));
    W[s]=(double *)malloc(LL->levels*sizeof(double ));
    prob[s]=(double *)malloc((LL->levels+1)*sizeof(double ));
    prob[s][LL->levels]=1.0;

    for(i=0; i<LL->levels; i++) 
      ionEnergy[s][i]=LL->ionEnergy[i];

    if(LL->species!=Electron && LL->ionFinal==OFF) sCnt++; else ;
    LL=LL->next;
    s++;
  }

  sList=(int *)malloc(sCnt*sizeof(int ));
  LL=D->loadList;
  s=0; index=0;
  while(LL->next)
  {
    if(LL->species!=Electron && LL->ionFinal==OFF) {
      sList[index]=s; index++;
    }  else ;
    LL=LL->next;
    s++;
  }

  //initialize J
  j=k=0;
  if(D->fieldType==Split)  {
    for(i=0; i<iend+3; i++)  {
          D->JxOld[i][j][k]=D->Jx[i][j][k];
          D->JyOld[i][j][k]=D->Jy[i][j][k];
          D->JzOld[i][j][k]=D->Jz[i][j][k];
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  } else {
    for(i=0; i<iend+3; i++)  {  
          D->Jx[i][j][k]=0.0;
          D->Jy[i][j][k]=0.0;
          D->Jz[i][j][k]=0.0;
        }
  }

  for(i=istart; i<iend; i++)
      {
        for(index=0; index<sCnt; index++) {
          s=sList[index];

            p=particle[i][j][k].head[s]->pt;

            while(p)  {
              Edc=sqrt(p->E1*p->E1+p->E2*p->E2+p->E3*p->E3);
              initZ=p->charge;
              weight=p->weight;
 
              if(initZ<levels[s] && Edc>minA0[s])  {
                testProb=randomValue(1.0);
  
                for(Z=initZ; Z<levels[s]; Z++) {
                  Efield=Ea/Edc;
                  ionE=ionEnergy[s][Z];
                  n_eff=(Z+1)*sqrt(1.0/ionE);
                  c2=pow(5.43656/n_eff,2.0*n_eff)/(n_eff);

                  phase=2.0*Efield*pow(ionE,1.5);
                  W[s][Z]=wa*c2*(0.5*ionE)*pow(phase,2.0*n_eff-1)*exp(-1.0*phase/3.0);
                }
                RungeKutta(&W[s][0],&prob[s][0],10,dt,initZ,levels[s],flag);
flag=0;
                Z=initZ; ionLevel=-1;
                while(Z<levels[s]) {
                  if(prob[s][Z]<testProb && testProb<prob[s][Z+1]) {
                    ionLevel=Z;
                    break;
                  }
                  Z++;
                }
//lala
//if(ionLevel>=5) printf("prob=%g, testProb=%g, ionLevel=%d\n",prob,testProb,ionLevel);

                if(ionLevel>=initZ) { 
                  value=ionEnergy[s][ionLevel]*unitEBydt*rho[s]/Edc/Edc*weight;
                  p->charge=ionLevel+1;
                  //creating electrons
                  cnt=p->charge-initZ;
                  for(n=0; n<cnt; n++) {
                    New = (ptclList *)malloc(sizeof(ptclList));
                    New->next = particle[i][j][k].head[s-1]->pt;
                    particle[i][j][k].head[s-1]->pt = New;
                    x=randomValue(1.0); y=0.0; z=0.0;
                    New->x = x; New->oldX=i+x;
                    New->y = y; New->oldY=j+y;
                    New->z = z; New->oldZ=k+z;
                    New->p1=p->p1; New->p2=p->p2; New->p3=p->p3;
                    New->p1Old1=0.0; New->p2Old1=0.0; New->p3Old1=0.0;
                    New->p1Old2=0.0; New->p2Old2=0.0; New->p3Old2=0.0;
                    New->weight=weight;
                    New->charge=-1;
                    New->index=p->index;
                    New->core=p->core;

                    x=p->x;
                    x1=(x+0.5)-(int)(x+0.5);
                    i1=i+(int)(x+0.5)-1;
                    D->Jx[i1][j][k]+=(1.0-x1)*p->E1*value;
                    D->Jx[i1+1][j][k]+=x1*p->E1*value;
                    D->Jy[i][j][k]+=(1.0-x)*p->E2*value;
                    D->Jy[i+1][j][k]+=x*p->E2*value;
                    D->Jz[i][j][k]+=(1.0-x)*p->E3*value;
                    D->Jz[i+1][j][k]+=x*p->E3*value;

                  }
                } else ;	//End of electron generation 

              }  	//End of if(initZ<levels[s])

              p=p->next;

            }  //End of while(p)


        }      //End of for(s)
      }        //End of for(i,j,k)
          
  free(sList);
  for(n=0; n<nSpecies; n++) {
    free(ionEnergy[n]); free(W[n]); free(prob[n]);
  } free(ionEnergy); free(W); free(prob);
}



double calW(double ionEnergy,double charge,double omega)
{
  double c2,n;
  double Ea=5.1e11, wa=4.13e16, euler=2.71828;
  
  Ea=eCharge*Ea/eMass/omega/velocityC;
  n=charge*sqrt(1.0/ionEnergy);

  c2=0.5/M_PI/n*pow(2*euler/n,2*n);
}
