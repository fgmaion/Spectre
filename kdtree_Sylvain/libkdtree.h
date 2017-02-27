// logarithmic binning in r.
double   zerolog;
double    maxlog;
int     nlogbins;
double  logbinsz;                 
// double  startlog;

// linear binning in mu. 
double   zerolin;
double    maxlin;
int     nlinbins;
double  linbinsz;

char filename[200];

int   indi; 
int   indj;

int   bruteforceCount = 0;

double pairCount = 0.0;

int   Nmin   =  200;      // Minimum number of particles which a Node can contain, should be > 1.

int   sortDim, treelevel;
int   NLeft, NRight;

int   tree_N, firstSplitDim, tree_labelCount;
  
double xmin[NDIM], xmax[NDIM];
  
double SplitValue;

double **gg, **rr, **gr, **landy_xi, *xi0, *xi2, *xi4, *logrbins, ngg, ngr, nrr;

double **gg_meanr, **gg_meanmu;
double **gr_meanr, **gr_meanmu;
double **rr_meanr, **rr_meanmu;

int*   pairsperbin;

double rmax[NDIM], rmin[NDIM], lmax[NDIM], lmin[NDIM];

typedef struct Particle{
  // int    SplitDim;
  
  double x[NDIM];
  
  int    index;
  
  double weight;
} Particle;

Particle*  point_gals;  
Particle* point_rands;

typedef struct Node{
  int              N;
  
  int          label;
  
  int      treelevel;
  
  int       SplitDim;
  
  int       Children;
  
  double  SplitValue;
  
  double  xmin[NDIM];
  
  double  xmax[NDIM];
  
  Particle* particle;
  
  void*         Left;
  void*        Right;
  
  int          NLeft;
  int         NRight;      
} Node;

Node*  galTree;
Node* randTree;

int    NodeStart_Indices[1000];

int    print_xi();

int    print_rr();

int    print_dd();

int    print_dr();

int    print_xiMultipoles();

int    CountPairs_rMu(double **C, double **r, double **mu, Node* firstTree, Node* secndTree, int sameTree);

int    assignLeafValues(Particle cat[], double xCoors[], double yCoors[], double zCoors[], int N);

Node*  buildTree(Particle cat[],  int N);

Node   *createNode(Particle cat[], int N, int SplitDim, double xmin[NDIM], double xmax[NDIM]);

int    sortby_position_alongDim_splitDim(const void *a, const void *b);

int    findSuitableNodePairs_bruteforcePairCount(double **C, double **r, double **mu, Node *node1, Node *node2, int sameTree);

int    bruteforceCountpairs_betweenNodes(double **C, double **r, double **mu, Node *node1, Node *node2, int sameTree);

double pair_zmu(Particle a, Particle b);

double log10_particleSeparation(Particle a, Particle b);

double log10_minimum_modDisplacementBetweenNodes(Node *node1, Node *node2);

double log10_maximum_modDisplacementBetweenNodes(Node *node1, Node *node2);

int    free_tree(Node *t);

double computeNorm(int Nfirst, int Nsecnd);
