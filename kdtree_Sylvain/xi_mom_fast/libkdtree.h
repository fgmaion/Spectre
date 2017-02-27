#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define FAILURE 0
#define SUCCESS 1

#define INF 1.0e30
#define ODD 0
#define EVEN 1

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN

#define NFIELD 200
#define NCHAR 40

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)  atoi(array+NCHAR*(col-1))
#define getLine(array,i) array+NFIELD*NCHAR*i

#define NDIM  3

/*----------------------------------------------------------------*
 *New types                                                       *
 *----------------------------------------------------------------*/

typedef struct Point
{
  int SplitDim;
  double x[NDIM];
  size_t index;
  double weight;
} Point;

typedef struct Node
{
  size_t N;
  int SplitDim;
  double SplitValue;
  double xmin[NDIM];
  double xmax[NDIM];
  Point *point;
  void *Left;
  void *Right;
} Node;

size_t IDERR;

/*----------------------------------------------------------------*
 *Tree routines                                                   *
 *----------------------------------------------------------------*/

size_t CountPt(int *C, Point o, double r1, double r2, Node *node);
size_t indexInNode(int *C,Point point,double r1,double r2,Node *node);
size_t indexInNodeAll(int *C,Point point,double r1,double r2,Node *node);
size_t check_kdt_point(double *ra,double *dec,int ngal,double pra,double pdec,double prad,int *C);

void comp_kdt_pairs(double *ra,double *dec,int ngal,double *ra_rand,double *dec_rand,int nrand,
		    double zerolog,int nlogbins,double logbinsz,double *C);
void comp_kdt_pairs_loglin(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			   double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,double linbinsz,double **C,double *norm);
void comp_kdt_pairs_loglin_smu(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			   double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,double linbinsz,double **C,double *norm);
void comp_kdt_pairs_linlin(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			   double zerolin,int nlinbins,double linbinsz,double **C,double *norm);
void comp_kdt_pairs_linlin_i(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			     double zerolin,int nlinbins,double linbinsz,double **C,double *norm,double ***CC,int nsub,double lbox,double dred);
void comp_kdt_pairs_log(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			double zerolog,int nlogbins,double logbinsz,double *C,double *norm);

void comp_kdt_pairs_loglin_pc(double *ra,double *dec,double *dc,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,int nrand,
			      double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,double linbinsz,double **C,double *norm);

void comp_kdt_pairs_w(double *ra,double *dec,double *wei,int ngal,double *ra_rand,double *dec_rand,double *wei_rand,int nrand,
		      double zerolog,int nlogbins,double logbinsz,double *C,double *norm);
void comp_kdt_pairs_loglin_w(double *ra,double *dec,double *dc,double *wei,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,
			     double *wei_rand,int nrand, double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,
			     double linbinsz,double **C,double *norm);
void comp_kdt_pairs_loglin_smu_w(double *ra,double *dec,double *dc,double *wei,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,
			     double *wei_rand,int nrand, double zerolog,int nlogbins,double logbinsz,double zerolin,int nlinbins,
			     double linbinsz,double **C,double *norm);
void comp_kdt_pairs_linlin_w(double *ra,double *dec,double *dc,double *wei,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,
			     double *wei_rand,int nrand,double zerolin,int nlinbins,double linbinsz,double **C,double *norm);
void comp_kdt_pairs_log_w(double *ra,double *dec,double *dc,double *wei,int ngal,double *ra_rand,double *dec_rand,double *dc_rand,
			  double *wei_rand,int nrand,double zerolog,int nlogbins,double logbinsz,double *C,double *norm);
void comp_kdt_density(double *ra,double *dec,double *dc,int ngal,double radius,double lbox,double dred,double *C,double *norm);
void comp_kdt_density_r(double *rarand,double *decrand,double *dcrand,int nrand,double *ra,double *dec,double *dc,int ngal,double radius,double lbox,double dred,double *C,double *norm);

size_t rangeCountPt(Point o, double r1, double r2, Node *node);
size_t countInNode(Point point,double r1,double r2,Node *node);

size_t rangeCount(double r1, double r2, Node *node1, Node *node2);
size_t countInNodes(double r1, double r2, Node *node1, Node *node2);

double rangeCountw(double r1, double r2, Node *node1, Node *node2);
double countInNodesw(double r1,double r2,Node *node1,Node *node2);
double countWeightsInNodes(Node *node1,Node *node2) ;

size_t rangeCount2(double *C, double zerolog, int nlogbins, double logbinsz, Node *node1, Node *node2);
size_t countInNodes2(double *C, double zerolog, int nlogbins, double logbinsz, Node *node1,Node *node2);

size_t rangeCount2w(double *C, double zerolog, int nlogbins, double logbinsz, Node *node1, Node *node2);
size_t countInNodes2w(double *C, double zerolog, int nlogbins, double logbinsz, Node *node1, Node *node2);

size_t rangeCountCent(double r1, double r2, Node *node1, Node *node2);
size_t countInNodesCent(double r1, double r2, Node *node1, Node *node2);

size_t rangeCount2d0(double r1, double r2, double p1, double p2, Node *node1, Node *node2);
size_t countInNodes2d0(double r1, double r2, double p1, double p2, Node *node1, Node *node2);

size_t rangeCount2d_loglin(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_loglin(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_loglin_smu(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_loglin_smu(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_loglin_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
			     double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_loglin_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
			       double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_loglin_smu_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
			     double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_loglin_smu_w(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, 
			       double linbinsz, Node *node1, Node *node2);

size_t rangeCount2dk(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2dk(double **C, double zerolog, int nlogbins, double logbinsz, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_linlin(double **C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_linlin(double **C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_linlin_w(double **C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_linlin_w(double **C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);

size_t rangeCount2d_linlin_i(double **C, double ***CC, int nsub, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_linlin_i(double **C, double ***CC, int nsub, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2d_linlin_i2(double ***CC, int nsub, int indi, int indj, Node *node1, Node *node2);
size_t normCount2d_linlin_i(double ***CC, int nsub, Node *node1, Node *node2);

size_t rangeCountR(double *C, double zerolog, int nlogbins, double logbinsz,Node *node1, Node *node2);
size_t countInNodesR(double *C, double zerolog, int nlogbins, double logbinsz,Node *node1,Node *node2);
size_t countInNodesRf(double *C, int ind, Node *node1,Node *node2);

size_t rangeCount2dR(double ***C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2dR(double ***C, double zerolin, int nlinbins, double linbinsz, Node *node1, Node *node2);
size_t countInNodes2dRf(double ***C, int indi, int indj, Node *node1, Node *node2);

double distPoint(Point a, Point b);
double distNodeMin(Point a, Node *node);
double distNodeMax(Point a, Node *node);
double distNodeNodeMin(Node *node1, Node *node2);
double distNodeNodeMax(Node *node1, Node *node2);

double distPointrp(Point a, Point b);
double distPointpi(Point a, Point b);
double distNodeNodeMinpi(Node *node1, Node *node2);
double distNodeNodeMaxpi(Node *node1, Node *node2);

double distPointmu(Point a, Point b);

double distPointpc(Point a, Point b);
int distNodeNodePc(Node *node1,Node *node2);

int insideNode(Point point,Node *node);
Point nearestPointNode(Point point,Node *node);

Node *buildTree(Point *field,size_t N);
Node *createNode(Point *field, size_t N, int SplitDim, double xmin[NDIM], double xmax[NDIM]);
int free_tree(Node *t);

/*----------------------------------------------------------------*
 *Initialization routines                                         *
 *----------------------------------------------------------------*/

int readPars(int argc, char **argv, FILE **fileIn1, int *x0col1, int *x1col1, int *x2col1, FILE **fileIn2, int *x0col2, int *x1col2, int *x2col2);
void fillPoint(double *rain, double *decin, int N, Point **field);
void fillPointw(double *rain, double *decin, double *weightin, int N, Point **field);
void fillPoint3d(double *rain, double *decin, double *dcin, int N, Point **field);
void fillPoint3dw(double *rain, double *decin, double *dcin, double *weightin, int N, Point **field);
void fillPoint3d_i(double *rain, double *decin, double *dcin, int N, Point **field,int nsub, double lbox, double dred);

/*----------------------------------------------------------------*
 *Utils - Numeric                                                 *
 *----------------------------------------------------------------*/

double round(double x);
int comparePoints(const void *a,const void *b);
void cpyPoints(Point *a, Point *b);
int compareDoubles(const void *a,const void *b);
int compareDoublesStable(const void *a,const void *b);
double determineMachineEpsilon();
size_t determineSize_tError();
FILE *fopenAndCheck(const char *filename,char *mode);
int getStrings(char *line,char *strings, char *delimit, size_t *N);
int roundToNi(double a);
