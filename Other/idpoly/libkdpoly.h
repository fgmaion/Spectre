#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#define PI    3.14159265358979323846
#define TWOPI 6.283185307179586476925287

#define FAILURE 0
#define SUCCESS 1

#define INF 1.0e30
#define ODD 0
#define EVEN 1
#define LEAF 0
#define NODE 1
#define RADEC 0
#define CART 1

#define NVERTICES 100

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN

#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)  atoi(array+NCHAR*(col-1))
#define getCharValue(array,col) array+NCHAR*(col-1)
#define getLine(array,i) array+NFIELD*NCHAR*i

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define SWAP(a,b) {swap = (a); (a) = (b); (b) = swap;}

char MYNAME[100];
size_t IDERR;

/*----------------------------------------------------------------*
 *New types                                                       *
 *----------------------------------------------------------------*/

typedef struct Complex
{
    double re;
    double im;
} Complex;

typedef struct Polygon
{
  int N, id;
  double *x; //[NVERTICES];
  double *y; //[NVERTICES];
  double *xmin; //[2];
  double *xmax; //[2];
} Polygon;


typedef struct NodeP
{
  int type, *root, id, *polysAll, SplitDim;
  double SplitValue;
  size_t Nnodes, Npolys, NpolysAll;
  int *poly_id;
  void *Left, *Right;
} NodeP;

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

NodeP *init_poly(char *file, double *x0);
int inPoly(NodeP *polyTree,double *x0,double ra,double dec);
int idPoly(NodeP *polyTree,double *x0,double ra,double dec);
void free_Polygon(Polygon *polygon, size_t N);
void free_NodeP(NodeP *node);
void cpyPolygon(Polygon *a, Polygon *b);

/*----------------------------------------------------------------*
 *Utils - geometric                                               *
 *----------------------------------------------------------------*/

int insidePolygon(Polygon *p,int Npoly, double x0,double y0,double x,double y, int *poly_id);
int insidePolygonTree(NodeP *polyTree, double x0[2], double x[2], int *poly_id);
Polygon *readPolygonFile(FILE *fileIn, int *Npolys, NodeP *polyTree);
NodeP *readPolygonFileTree(FILE *fileIn, double xmin[2], double xmax[2]);
NodeP *createNodeP(Polygon *polys, size_t Npolys, double minArea, int SplitDim, double xmin[2], double xmax[2], int firstCall);

/*----------------------------------------------------------------*
 *Utils - numeric                                                 *
 *----------------------------------------------------------------*/

int getStrings(char *line,char *strings, char *delimit, size_t *N);
void printCount(const size_t *count, const size_t *total,  const size_t step);
