#include "libkdtree.h"

#define deg2rad M_PI/180

//-----------------------------------------------------------------------

void comp_kdt_pairs_loglin_smu(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			       double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,double linbinsz,double **C,double *norm) 
{
  Point *point1,*point2;

  // convert to cartesian coordinates. 
  fillPoint3d(ra,dec,dc,ngal,&point1);
  fillPoint3d(ra_rand,dec_rand,dc_rand,nrand,&point2);

  Node *point1Tree = buildTree(point1,ngal);
  Node *point2Tree = buildTree(point2,nrand);

  // counts pairs. for a range in r. 
  rangeCount2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,point1Tree,point2Tree);
  // binning. 2D two bins. make a 3D version of this. input as trees, 2 can point to the same tree.

  if (ngal==nrand) *norm=(double)ngal*((double)ngal-1.0);
  else *norm=(double)ngal*(double)nrand;

  free(point1); 
  free(point2);
  free_tree(point1Tree);
  free_tree(point2Tree);
}

// same as above, but with weight. 
void comp_kdt_pairs_loglin_smu_w(double *ra,double *dec,double *dc,double *wei,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,
			     double *wei_rand,int nrand, double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,
			     double linbinsz,double **C,double *norm) 
{
    int i,j;
    Point *point1,*point2;

    if (fabs(wei[1]-dc[1])<1e-30) fillPoint3d(ra,dec,dc,ngal,&point1);
    else fillPoint3dw(ra,dec,dc,wei,ngal,&point1);
    
    if (fabs(wei_rand[1]-dc_rand[1])<1e-30) fillPoint3d(ra_rand,dec_rand,dc_rand,nrand,&point2);
    else fillPoint3dw(ra_rand,dec_rand,dc_rand,wei_rand,nrand,&point2);
    
    Node *point1Tree = buildTree(point1,ngal);
    Node *point2Tree = buildTree(point2,nrand);

    rangeCount2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,point1Tree,point2Tree);
    // C stores the number of pairs in each bin, nlogbins x nlinbins. 
    
    
    *norm=0;
    if (point1Tree->N!=point2Tree->N) for(i=0;i<point1Tree->N;i++) for(j=0;j<point2Tree->N;j++) *norm+=point1Tree->point[i].weight*point2Tree->point[j].weight;
    else for(i=0;i<point1Tree->N;i++) for(j=i+1;j<point2Tree->N;j++) *norm+=2*point1Tree->point[i].weight*point2Tree->point[j].weight;
    
    free(point1); 
    free(point2);
    free_tree(point1Tree);
    free_tree(point2Tree);
}

//-----------------------------------------------------------------------

size_t rangeCount2d_loglin_smu(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
			       double linbinsz, Node *node1, Node *node2) 
{
  int i,j; //,ii,jj;

    i=1+(int)floor((0.5*log10(distNodeNodeMin(node1,node2))-zerolog)/logbinsz);
    j=1+(int)floor((0.5*log10(distNodeNodeMax(node1,node2))-zerolog)/logbinsz);
    //ii=1+(int)floor((sqrt(distNodeNodeMin(node1,node2))-zerolin)/linbinsz);
    //jj=1+(int)floor((sqrt(distNodeNodeMax(node1,node2))-zerolin)/linbinsz);

    
    if (i>nlogbins || j<1) return 0;
    /*
    else if (i==j && ii==jj) {
      C[i][ii]+=node1->N*node2->N;
      return 0;
    }
    */

    if(node1->Left == NULL && node1->Right == NULL && node2->Left == NULL && 
       node2->Right == NULL) return countInNodes2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1,node2);
    else {
	if(node2->Left == NULL && node2->Right == NULL) {
	    rangeCount2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1->Left,node2);
	    rangeCount2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1->Right,node2);
	    return 0;
	} else {
	    rangeCount2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node2->Left,node1);
	    rangeCount2d_loglin_smu(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node2->Right,node1);
	    return 0;
	}
    }
}

size_t countInNodes2d_loglin_smu(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
				 double linbinsz, Node *node1, Node *node2) 
{
    size_t i,j,indi,indj;
    double dTmps,dTmpmu;
    logbinsz=1.0/logbinsz;
    linbinsz=1.0/linbinsz;
    
    for(i=0;i<node1->N;i++) {
	for(j=0;j<node2->N;j++) {
	    dTmps = distPoint(node1->point[i],node2->point[j]);
	    dTmpmu = distPointmu(node1->point[i],node2->point[j]);
	    indi=1+(int)floor((0.5*log10(dTmps)-zerolog)*logbinsz);
	    indj=1+(int)floor((sqrt(dTmpmu)-zerolin)*linbinsz);
	    if (indi>=1 && indi<=nlogbins && indj>=1 && indj<=nlinbins) C[indi][indj]++;
	} 
    }
    return 0;
}

//-----------------------------------------------------------------------

size_t rangeCount2d_loglin_smu_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
				 double linbinsz, Node *node1, Node *node2) 
{
  int i,j; //,ii,jj;

    // This will need to be changed. 

    i=1+(int)floor((0.5*log10(distNodeNodeMin(node1,node2))-zerolog)/logbinsz);
    j=1+(int)floor((0.5*log10(distNodeNodeMax(node1,node2))-zerolog)/logbinsz);
    //ii=1+(int)floor((sqrt(distNodeNodeMin(node1,node2))-zerolin)/linbinsz);
    //jj=1+(int)floor((sqrt(distNodeNodeMax(node1,node2))-zerolin)/linbinsz);

    if (i>nlogbins || j<1) return 0;
    /*
    else if (i==j && ii==jj) {
      C[i][ii]+=countWeightsInNodes_smu(node1,node2);
      return 0;
    }
    */

    // IF BOTH ARE LEAVES **
    if(node1->Left == NULL && node1->Right == NULL && node2->Left == NULL && 
       node2->Right == NULL) return countInNodes2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1,node2);
    else {
	if(node2->Left == NULL && node2->Right == NULL) {
	    rangeCount2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1->Left,node2);
	    rangeCount2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node1->Right,node2);
	    return 0;
	} else {
	    rangeCount2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node2->Left,node1);
	    rangeCount2d_loglin_smu_w(C,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,node2->Right,node1);
	    return 0;
	}
    }
}

size_t countInNodes2d_loglin_smu_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
 				   double linbinsz, Node *node1, Node *node2) 
{   // ** Brute force calculation for leaves **.  This will need changed. 
    size_t i,j,indi,indj;
    double dTmps,dTmpmu;
    logbinsz=1.0/logbinsz;
    linbinsz=1.0/linbinsz;
    
    for(i=0;i<node1->N;i++) {
	for(j=0;j<node2->N;j++) {
	    dTmps = distPoint(node1->point[i],node2->point[j]);
	    dTmpmu = distPointmu(node1->point[i],node2->point[j]);    
	    indi=1+(int)floor((0.5*log10(dTmps)-zerolog)*logbinsz);
	    indj=1+(int)floor((sqrt(dTmpmu)-zerolin)*linbinsz);
	    if (indi>=1 && indi<=nlogbins && indj>=1 && indj<=nlinbins) C[indi][indj]+=node1->point[i].weight*node2->point[j].weight;
	} 
    }
    return 0;
}

//-----------------------------------------------------------------------

double distPoint(Point a, Point b){
  int i;
  
  double sum = 0.0;
  for (i=0;i<NDIM;i++) sum += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
  return sum;
}

double distPointmu(Point a, Point b){
  int i;
  
  double sum = 0.0, suma = 0.0,sumb = 0.0;

  for (i=0;i<NDIM;i++) {
    sum  += (a.x[i]-b.x[i])*(a.x[i]-b.x[i]);
    suma += a.x[i]*a.x[i]-b.x[i]*b.x[i];
    sumb += (a.x[i]+b.x[i])*(a.x[i]+b.x[i]);
  }
  return suma*suma/sumb/sum;
}

//-----------------------------------------------------------------------

double distNodeNodeMin(Node *node1, Node *node2){
  int i;
  
  double sum = 0.0;
  for (i=0;i<NDIM;i++) {

    if (node1->xmin[i] > node2->xmin[i] && node1->xmin[i] > node2->xmax[i]) 
      sum += (node1->xmin[i]-node2->xmax[i])*(node1->xmin[i]-node2->xmax[i]);

    if (node2->xmin[i] > node1->xmax[i] && node2->xmax[i] > node1->xmax[i])
      sum += (node1->xmax[i]-node2->xmin[i])*(node1->xmax[i]-node2->xmin[i]);
  }  
  return sum;
}

double distNodeNodeMax(Node *node1, Node *node2){
  int i;
  
  double sum = 0.0,d_lo,d_hi,dx;
  for (i=0;i<NDIM;i++) {
    
    d_hi = (node1->xmax[i]-node2->xmin[i]);
    d_lo = (node2->xmax[i]-node1->xmin[i]);
    dx=MAX(d_lo,d_hi);

    sum += dx*dx;
  }
  return sum;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*
 * Tree initialization routines                                              *
 *---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------- *
 *  Creates a tree for the set of points "point" and returns a pointer       *
 *  on the tree root.                                                        *
 * ------------------------------------------------------------------------- */

Node *buildTree(Point *point,size_t N)
{
  size_t i,j;
  int SplitDim;
  double xmin[NDIM],xmax[NDIM];
  Point *pointTmp = (Point *)malloc(N*sizeof(Point));
  
  for(i=0;i<N;i++) {
    for(j=0;j<NDIM;j++) pointTmp[i].x[j] = point[i].x[j];
    //pointTmp[i].index = i;
    pointTmp[i].index = point[i].index;
    pointTmp[i].weight = point[i].weight;
  }
  
  //Computes limits and fills indexes
  for(i=0;i<N;i++) {
    if(i==0) {
      for(j=0;j<NDIM;j++) {
	xmin[j]=pointTmp[i].x[j];
	xmax[j]=pointTmp[i].x[j];
      }
    } else {
      for(j=0;j<NDIM;j++) {
	xmin[j]=MIN(pointTmp[i].x[j],xmin[j]);
	xmax[j]=MAX(pointTmp[i].x[j],xmax[j]);
      }
    }
  }
  
  //Takes the larger coordinates for the first cut
  SplitDim = 0;
  for(j=1;j<NDIM;j++) if( (xmax[j] - xmin[j]) > (xmax[0] - xmin[0])) SplitDim = j;
  
  //Necessary for comparePoints()
  for(i=0;i<N;i++) pointTmp[i].SplitDim = SplitDim;
  
  return createNode(pointTmp,N,SplitDim,xmin,xmax);
}

/* ------------------------------------------------------------------------- *
 *    Computes SplitValue and returns a pointer on the node, recursively.    *
 *    If N < Nmin returns a node with                                        *
 *    no child => node.Left = NULL and node.Right = NULL <=> leaf.           *
 *    To build a tree, just gives the whole point, SplitDim, xmin and xmax.  *
 * ------------------------------------------------------------------------- */

Node *createNode(Point *point, size_t N, int SplitDim, double xmin[NDIM], double xmax[NDIM])
{
  int NewSplitDim;
  size_t i,j,Nmin = 200; //Nmin should be > 1
  size_t NLeft,NRight;
  double SplitValue;
  double xminLeft[NDIM],xmaxLeft[NDIM];
  double xminRight[NDIM],xmaxRight[NDIM];
  
  //Sorts the values
  qsort(point,N,sizeof(Point),comparePoints);

  //Finds SplitValue [using the median along the maxspread dimension]
  switch(PARITY(N)){
  case EVEN:
    NLeft = N/2;
    NRight = N/2;
    break;
  case ODD:
    NLeft = (N+1)/2;
    NRight = (N-1)/2;
    break;
  }
  SplitValue = (point[NLeft-1].x[SplitDim] + point[NLeft].x[SplitDim])/2.0;
 
  //Allocates memory for THIS node
  Node *result = (Node *)malloc(sizeof(Node));

  //Writes N, SplitDim, and SplitValue
  result->N = N;
  result->SplitDim = SplitDim;
  result->SplitValue = SplitValue;

  /*printf("NLeft= %d,NRight= %d\n",NLeft,NRight);*/

  //Limits
  for(j=0;j<NDIM;j++){
    result->xmin[j] = xmin[j];
    result->xmax[j] = xmax[j];
    xminLeft[j] = xmin[j];
    xmaxLeft[j] = xmax[j];
    xminRight[j] = xmin[j];
    xmaxRight[j] = xmax[j];
  }

  result->point = point;
  
  if(N < Nmin) {
      result->Left = NULL;
      result->Right = NULL;
  } else {
      //Sets the next coordinate
      NewSplitDim = SplitDim + 1 ;
      if(NewSplitDim > NDIM-1)  NewSplitDim = 0;
      
      //Left
      Point *pointLeft = (Point *)malloc(NLeft*sizeof(Point));
      for(i=0;i<NLeft;i++){
	  for(j=0;j<NDIM;j++) pointLeft[i].x[j] = point[i].x[j];
	  pointLeft[i].index = point[i].index;
	  pointLeft[i].SplitDim = NewSplitDim;
	  pointLeft[i].weight = point[i].weight;
      }
      xmaxLeft[SplitDim] = SplitValue;
      
      //Pointer to the left node
      result->Left = createNode(pointLeft,NLeft,NewSplitDim,xminLeft,xmaxLeft);
      
      //Right
      Point *pointRight = (Point *)malloc(NRight*sizeof(Point));
      for(i=NLeft;i<N;i++){
	  for(j=0;j<NDIM;j++) pointRight[i-NLeft].x[j] = point[i].x[j];
	  pointRight[i-NLeft].index = point[i].index;
	  pointRight[i-NLeft].SplitDim = NewSplitDim;
	  pointRight[i-NLeft].weight = point[i].weight;
      }
      xminRight[SplitDim] = SplitValue;
      
      //Pointer to the right node
      result->Right = createNode(pointRight,NRight,NewSplitDim,xminRight,xmaxRight);
  }
  
  return result;
}

int free_tree(Node *t)
{
  free(t->point);
  if (t->Left!=NULL) free_tree(t->Left);
  if (t->Right!=NULL) free_tree(t->Right);
  free(t);
  return 0;
}

/*----------------------------------------------------------------*
 *Initialization routines                                         *
 *----------------------------------------------------------------*/

int readPars(int argc, char **argv, FILE **fileIn1, int *x0col1, int *x1col1, int *x2col1, FILE **fileIn2, int *x0col2, int *x1col2, int *x2col2) 
{
    int i;
    
    for(i=0;i<argc;i++){
	if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc < 9){
	    printf("Usage:  %s fileIn1 id_x id_y id_z fileIn2 id_x id_y id_z\n",argv[0]);
	    exit(-1);
	}
    }
    
    *fileIn1 = fopenAndCheck(argv[1],"r");   
    *x0col1 = atoi(argv[2]);
    *x1col1 = atoi(argv[3]);
    *x2col1 = atoi(argv[4]);
    
    *fileIn2 = fopenAndCheck(argv[5],"r");
    *x0col2 = atoi(argv[6]);
    *x1col2 = atoi(argv[7]);
    *x2col2 = atoi(argv[8]);
    
    return SUCCESS;
}

void fillPoint(double *rain, double *decin, int N, Point **field)
{
    size_t i;
    double ra,dec;
    
    *field = (Point *)malloc(N*sizeof(Point));
    
    for (i=0;i<N;i++) {      
	ra  = rain[i+1];
	dec = decin[i+1];
	
	(*field)[i].SplitDim = -1;
	(*field)[i].x[0] = cos(dec)*cos(ra);
	(*field)[i].x[1] = cos(dec)*sin(ra);
	(*field)[i].x[2] = sin(dec);
	(*field)[i].index = i;
	(*field)[i].weight = 1.0;
	
	if ((*field)[i].x[0] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[0]);
	if ((*field)[i].x[1] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[1]);
	if ((*field)[i].x[2] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[2]);
    }
}

void fillPointw(double *rain, double *decin, double *weightin, int N, Point **field)
{
    size_t i;
    double ra,dec;
    
    *field = (Point *)malloc(N*sizeof(Point));
    
    for (i=0;i<N;i++) {      
	ra  = rain[i+1];
	dec = decin[i+1];
	
	(*field)[i].SplitDim = -1;
	(*field)[i].x[0] = cos(dec)*cos(ra);
	(*field)[i].x[1] = cos(dec)*sin(ra);
	(*field)[i].x[2] = sin(dec);
	(*field)[i].index = i;
	(*field)[i].weight =  weightin[i+1];
	
	if ((*field)[i].x[0] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[0]);
	if ((*field)[i].x[1] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[1]);
	if ((*field)[i].x[2] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[2]);
    }
}

void fillPoint3d(double *rain, double *decin, double *dcin, int N, Point **field)
{
    size_t i;
    double ra,dec,dc;
    
    *field = (Point *)malloc(N*sizeof(Point));
    
    for (i=0;i<N;i++) {      
	ra  = rain[i+1];
	dec = decin[i+1];
	dc = dcin[i+1];
	
	(*field)[i].SplitDim = -1;
	(*field)[i].x[0] = dc*cos(dec)*cos(ra);
	(*field)[i].x[1] = dc*cos(dec)*sin(ra);
	(*field)[i].x[2] = dc*sin(dec);
	(*field)[i].index = i;
	(*field)[i].weight = 1.0;
	
	if ((*field)[i].x[0] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[0]);
	if ((*field)[i].x[1] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[1]);
	if ((*field)[i].x[2] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[2]);
    }
}

void fillPoint3dw(double *rain, double *decin, double *dcin, double *weightin, int N, Point **field)
{
    size_t i;
    double ra,dec,dc;
    
    *field = (Point *)malloc(N*sizeof(Point));
    
    for (i=0;i<N;i++) {      
	ra  = rain[i+1];
	dec = decin[i+1];
	dc = dcin[i+1];
	
	(*field)[i].SplitDim = -1;
	(*field)[i].x[0] = dc*cos(dec)*cos(ra);
	(*field)[i].x[1] = dc*cos(dec)*sin(ra);
	(*field)[i].x[2] = dc*sin(dec);
	(*field)[i].index = i;
	(*field)[i].weight = weightin[i+1];
	
	if ((*field)[i].x[0] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[0]);
	if ((*field)[i].x[1] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[1]);
	if ((*field)[i].x[2] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[2]);
    }
}

void fillPoint3d_i(double *rain, double *decin, double *dcin, int N, Point **field,int nsub, double lbox, double dred)
{
  size_t i,indx,indy,indz,ind;
  double ra,dec,dc;
  
  *field = (Point *)malloc(N*sizeof(Point));
  
  for (i=0;i<N;i++) {      
    ra  = rain[i+1];
    dec = decin[i+1];
    dc  = dcin[i+1];
    
    indx=1+(int)floor((ra+0.5*lbox)*nsub/lbox);
    indy=1+(int)floor((dec+0.5*lbox)*nsub/lbox);
    indz=1+(int)floor((dc+0.5*lbox-dred)*nsub/lbox);
    ind=indx+(indy-1)*nsub+(indz-1)*nsub*nsub;
    if (ind<1 || ind>nsub*nsub*nsub) printf("Error while indexing objects!\n");
	
    (*field)[i].SplitDim = -1;
    (*field)[i].x[0] = dc*cos(dec)*cos(ra);
    (*field)[i].x[1] = dc*cos(dec)*sin(ra);
    (*field)[i].x[2] = dc*sin(dec);
    (*field)[i].index = ind;
    (*field)[i].weight = 1.0;
    
    if ((*field)[i].x[0] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[0]);
    if ((*field)[i].x[1] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[1]);
    if ((*field)[i].x[2] > 1e8) printf("Warning: %lf value !\n",(*field)[i].x[2]);
  }
}

/*----------------------------------------------------------------*
 * Utils - Numeric                                                *
 *----------------------------------------------------------------*/


// ** Points are really particles!  **
int comparePoints(const void *a,const void *b)
{
    /* Compares x0 coordinates of 2 points */
 
    int SplitDim = (*(Point *)a).SplitDim;
    
    if ((*(Point *)a).x[SplitDim] > (*(Point *)b).x[SplitDim]) return 1;
    else if ((*(Point *)a).x[SplitDim] < (*(Point *)b).x[SplitDim]) return -1;
    else return 0;
}

void cpyPoints(Point *a, Point *b)
{
    /* Copies b into a */
    int j;
    
    a->SplitDim = b->SplitDim;
    for(j=0;j<NDIM;j++) a->x[j] = b->x[j];
    a->index = b->index;
}

int compareDoubles(const void *a,const void *b)
{
    /* Compares two double precision numbers */
    if (*(double *)a > *(double *)b) return 1;
    else if (*(double *)a < *(double *)b) return -1;
    else return 0;
}

double determineMachineEpsilon()
{
    double u, den;
    
    u = 1.0;
    do {
	u /= 2.0;
	den = 1.0 + u;
    } while (den>1.0);
    
    return (10.0 * u);
}

size_t determineSize_tError()
{
    /* Oui ok c'est un peu barbare, mais bon... */
    size_t count=0;
    count--;
    return count;
}

FILE *fopenAndCheck(const char *fileName,char *mode)
{
    /* Checks if fileName exists and opens it. Exits otherwise. */
    FILE *fileTmp = fopen(fileName,mode);
    if (fileTmp == NULL){
	printf("%s not found. Exiting...\n",fileName);
	exit(-1);    
    }
    return fileTmp;
}

int getStrings(char *line, char *strings, char *delimit, size_t *N)
{
    /* Extract each word/number in line separated by delimit and returns 
       the array of items in strings. */
    int i,j,begin,length;
    
    if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
    
    i = 0;
    j = 0;
    while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
	begin = i;
	while((line[i] == *delimit || line[i] == '\t') && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
	begin = i;
	while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
	length = i - begin;
	if(length > 0){
	    strncpy(strings+NCHAR*j,&line[begin],length);
	    strcpy(strings+NCHAR*j+length,"\0");
	    j++;
	}
    }
    
    (*N) = j;
    
    if(*N > 0) return 1;
    else return 0;
}

int roundToNi(double a){return (int)round(a);}
