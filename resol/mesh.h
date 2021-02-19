#define HIGH   1
#define LOW    2

typedef struct _Domain
{
  int dimension;
  int mode;
  int step;
  int L,M,N;
  int nextXrank,prevXrank,nextYrank,prevYrank,nextZrank,prevZrank;

  int resolX,resolY,resolZ;
  int rankX,rankY,rankZ;
  double *targetW;

  int nSpecies;
  int nx,ny,nz;
  int nxSub,nySub,nzSub;
  int minXSub,minYSub,minZSub,maxXSub,maxYSub,maxZSub;
  int reMinXSub,reMinYSub,reMinZSub,reMaxXSub,reMaxYSub,reMaxZSub;
  int minXDomain,minYDomain,minZDomain;
  int maxXDomain,maxYDomain,maxZDomain;

  int istart,jstart,kstart,iend,jend,kend;

  //resize
  int reNx,reNy,reNz;
  int reNxSub,reNySub,reNzSub;
  int reIend,reJend,reKend;
  int minX,maxX,minY,maxY,minZ,maxZ;
  int reMinX,reMaxX,reMinY,reMaxY,reMinZ,reMaxZ;
  int saveIstart,saveIend,saveNxSub,biasX;
  int saveJstart,saveJend,saveNySub,biasY;
  int saveKstart,saveKend,saveNzSub,biasZ;
  int reSaveIstart,reSaveIend,reSaveNxSub;
  int reSaveJstart,reSaveJend,reSaveNySub;
  int reSaveKstart,reSaveKend,reSaveNzSub;

  int *recv,*recvDataCnt,*sharePNum,*coreCnt;
  int *minXSubList,*minYSubList,*minZSubList;
  int *maxXSubList,*maxYSubList,*maxZSubList;
  int *reMinXSubList,*reMinYSubList,*reMinZSubList;
  int *reMaxXSubList,*reMaxYSubList,*reMaxZSubList;
  double **recvData,**sendData;

  //original field memory
  double ***fieldE,***fieldOld,***fieldNow,***fieldNext;
   
  //particle momory
  struct _Particle ***particle;
//  struct _Particle ***befoParticle;
//  struct _Particle ***nextParticle;

} Domain;

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _Particle
{
   ptclHead **head;
}  Particle;

typedef struct _ptclList  {
    double x;
    double y;
    double z;
    double px;    //momentum  
    double py;
    double pz;
    double pxSam;
    double pySam;
    double pzSam;
    int index;
    int core;
    double weight;
    double charge;
    struct _ptclList *next;
} ptclList;

typedef struct _queList  {
  double x;    double y;    double z;
  double px;   double py;   double pz;
  double gapX; double gapY; double gapZ;
  double weight;
  struct _queList *next;
} queList;

typedef struct _newP  {
  int flag;
  double x1; double y1;
  double px1; double py1; double pz1; double w1;
  double x2; double y2;
  double px2; double py2; double pz2; double w2;
} newP;

void boundary(Domain *D);
void saveParticle(Domain *D,int s,Particle ***particle);
void saveDumpParticleHDF(Domain *D);
void restoreParticleHDF(Domain *D,int iteration,Particle ***particle);
void cleanParticle(Domain *D);
void reCreateParticle(Domain D);
void reHighCreateField(Domain D);
void reLowCreateField(Domain D);
int whatMode(char *str);
void changeParticlePosition(Domain *D);



