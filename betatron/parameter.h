
typedef struct _Parameter 
{
   int dimension;

   double minThX;
   double maxThX;
   int numThX;
   double minThY;
   double maxThY;
   int numThY;

   int species;

   double dt;
   double superP;
   double *B;
   double *C;
   double **data;  
   double *share;
   double *shareBC;

   //here is detector information for spectrums.
   double det_R; 
   double det_x; 
   double det_y; 
   double minE;
   double maxE;
   double dE;
   double *freq;

   //mpi parameter
   int minIdList;
   int maxIdList;
}  Parameter; 

