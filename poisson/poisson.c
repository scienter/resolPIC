#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void sor(double **phi,double **f,double a,double b,double c,double d,double e,int nx,int ny,double rjac);

int main(int argc, char *argv[])
{
   int i,j,nx,ny;
   double dx,dy,x,y,r,cenX,cenY,rho,rjac;
   double a,b,c,d,e;
   double **den,**phi,**f;
   char name[100];
   FILE *out;

   nx=400; ny=80;
   dx=0.025; dy=0.25;
   cenX=nx*dx*0.5;
   cenY=ny*dy*0.5;
   phi=(double **)malloc(nx*sizeof(double *));
   den=(double **)malloc(nx*sizeof(double *));
   f=(double **)malloc(nx*sizeof(double *));
   for(i=0; i<nx; i++) {
     phi[i]=(double *)malloc(ny*sizeof(double ));
     den[i]=(double *)malloc(ny*sizeof(double ));
     f[i]=(double *)malloc(ny*sizeof(double ));
   }
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++) {
       phi[i][j]=0.0;
       den[i][j]=0.0;
       f[i][j]=0.0;
     }
  
   //source
   for(i=0; i<nx; i++) {
     x=i*dx;
     for(j=0; j<ny; j++) {
       y=j*dy;
       r=sqrt((x-cenX)*(x-cenX)+(y-cenY)*(y-cenY));
       if(r<3) den[i][j]=-0.1; else ;
     }
   }
   for(i=0; i<nx; i++) 
     for(j=0; j<ny; j++) 
       f[i][j]=-2*M_PI*den[i][j];

   a=1.0/dx/dx;
   b=1.0/dx/dx;
   c=1.0/dy/dy;
   d=1.0/dy/dy;
   e=-2.0*(1.0/dx/dx+1.0/dy/dy);

//   rjac=(cos(M_PI/(1.0*nx))+(dx*dx/dy/dy)*cos(M_PI/(1.0*ny)))/(1.0*dx*dx/dy/dy);
//   rjac=1.0-M_PI/(1.0*nx);
   rjac=0.99999;
   sor(phi,f,a,b,c,d,e,nx,ny,rjac);

   //save file
   sprintf(name,"den");
   out=fopen(name,"w");
   for(i=0; i<nx; i++) {
     x=i*dx;
     for(j=0; j<ny; j++) {
       y=j*dy;
       rho=den[i][j];
       fprintf(out,"%g %g %g\n",x,y,rho);
     }
     fprintf(out,"\n");
   }
   fclose(out);
   printf("%s is made.\n",name);

   sprintf(name,"phi");
   out=fopen(name,"w");
   for(i=0; i<nx; i++) {
     x=i*dx;
     for(j=0; j<ny; j++) {
       y=j*dy;
       rho=phi[i][j];
       fprintf(out,"%g %g %g\n",x,y,rho);
     }
     fprintf(out,"\n");
   }
   fclose(out);
   printf("%s is made.\n",name);



   for(i=0; i<nx; i++) { free(phi[i]); free(den[i]); free(f[i]); } free(phi); free(den); free(f);
   return 0;
}


void sor(double **phi,double **f,double a,double b,double c,double d,double e,int nx,int ny,double rjac)
{
   int ipass,i,isw,j,jsw,n,maxit=10000;
   double anorm,anormf=0.0,omega=1.0,resid,ratio=1e-5;

   for (i=1; i<nx-1; i++)
     for (j=1; j<ny-1; j++)
       anormf+=fabs(f[i][j]);

   anorm=anormf;
   n=1;
   while(anorm>anormf*ratio && n<maxit)
   {
     anorm=0.0;
     jsw=1;
     //first step
     for(i=1; i<nx-1; i++)  {
       for(j=jsw; j<ny-1; j+=2)  {
   	 resid=a*phi[i+1][j]
	        +b*phi[i-1][j]
	        +c*phi[i][j+1]
	        +d*phi[i][j-1]
	        +e*phi[i][j]
	        -f[i][j];
	 anorm+=fabs(resid);
	 phi[i][j]-=omega*resid/e;
       }
       jsw=3-jsw;
     }
     omega=(n==1 ? 1.0/(1.0-0.5*rjac*rjac) :
		       1.0/(1.0-0.25*rjac*rjac*omega));

     jsw=2;
     //first step
     for(i=1; i<nx-1; i++)  {
       for(j=jsw; j<ny-1; j+=2)  {
   	 resid=a*phi[i+1][j]
	        +b*phi[i-1][j]
	        +c*phi[i][j+1]
	        +d*phi[i][j-1]
	        +e*phi[i][j]
	        -f[i][j];
	 anorm+=fabs(resid);
	 phi[i][j]-=omega*resid/e;
       }
       jsw=3-jsw;
     }
     omega=1.0/(1.0-0.25*rjac*rjac*omega);

     n++;

     if(n%10==0) 
       printf("cnt=%d, anorm=%g, error=%g, rjac=%g, w=%g,dt=%g\n",n,anorm,ratio*anormf,rjac,omega,-omega*resid/e);
   }
}

//Classic method
/*
void sor(double **phi,double **f,double a,double b,double c,double d,double e,int nx,int ny,double rjac)
{
   int ipass,i,isw,j,jsw,n,maxit=1000;
   double anorm,anormf=0.0,omega=1.0,resid,ratio=1e-5,dt;
   double **Err,err;

   Err=(double **)malloc(nx*sizeof(double *));
   for (i=0; i<nx; i++)
     Err[i]=(double *)malloc(ny*sizeof(double ));
   for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
       Err[i][j]=0.0;

   dt=0.5/(a+c);
//   /(1.0+1.0/a+1.0/c);
   printf("a=%g, c=%g,dt=%g\n",a,c,dt);

   for (i=1; i<nx-1; i++)
     for (j=1; j<ny-1; j++)
       anormf+=fabs(f[i][j]);

   for(n=1; n<=maxit; n++) {
     anorm=0.0;
     for(i=1; i<nx-1; i++) 
       for(j=1; j<ny-1; j++)  {
         err=a*phi[i+1][j]+b*phi[i-1][j]+c*phi[i][j+1]+d*phi[i][j-1]+e*phi[i][j]-f[i][j];
         anorm+=fabs(err);
	 Err[i][j]=err;
       }

     for(i=1; i<nx-1; i++) 
       for(j=1; j<ny-1; j++)  
         phi[i][j]+=dt*Err[i][j];

     if(anorm<ratio*anormf) return;

     if(n%10==0) 
       printf("cnt=%d, anorm=%g, error=%g, rjac=%g, w=%g,dt=%g\n",n,anorm,ratio*anormf,rjac,omega,dt);
   }
}
*/

//original
/*
void sor(double **phi,double **f,double a,double b,double c,double d,double e,int nx,int ny,double rjac)
{
   int ipass,i,isw,j,jsw,n,maxit=10000;
   double anorm,anormf=0.0,omega=1.0,resid,ratio=1e-5;

   for (i=1; i<nx-1; i++)
     for (j=1; j<ny-1; j++)
       anormf+=fabs(f[i][j]);

   for(n=1; n<=maxit; n++) {
     anorm=0.0;
     isw=0;
     for(ipass=1; ipass<=2; ipass++) {
       jsw=isw;
       for(i=1; i<nx-1; i++)  {
         for(j=jsw+1; j<ny-1; j+=2)  {
   	   resid=a*phi[i+1][j]
		   +b*phi[i-1][j]
		   +c*phi[i][j+1]
		   +d*phi[i][j-1]
		   +e*phi[i][j]
		   -f[i][j];
	   anorm+=fabs(resid);
	   phi[i][j]-=omega*resid/e;
	 }
	 jsw=3-jsw;
       }
       isw=3-isw;
       omega=(n==1 && ipass==1 ? 1.0/(1.0-0.5*rjac*rjac) :
		       1.0/(1.0-0.25*rjac*rjac*omega));
     }
     if(anorm<ratio*anormf) return;

     if(n%10==0) 
       printf("cnt=%d, anorm=%g, error=%g, rjac=%g, w=%g,dt=%g\n",n,anorm,ratio*anormf,rjac,omega,-omega*resid/e);
   }
}
*/
