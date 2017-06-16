#define NDIM 3

typedef struct Particle{
  // int    SplitDim;

  double x[NDIM];

  double disp;

  int    index;

  double weight;
} Particle;

Particle* point_gals;
Particle* point_rands;

typedef struct Node{
  long             N;
  long         label;
  long     treelevel;
  int       SplitDim;
  int       Children;

  double  SplitValue;
  double  xmin[NDIM];
  double  xmax[NDIM];

  Particle* particle;

  void*         Left;
  void*        Right;

  long         NLeft;
  long        NRight;
} Node;


Node*    galTree;
Node*    randTree;

// logarithmic binning in r.
double   zerolog;
double   maxlog;
int      nlogbins;
double   logbinsz;                 

// linear binning in mu. 
double   zerolin;
double   maxlin;
int      nlinbins;
double   linbinsz;

char     filename[200];

int      indi; 
int      indj;

long     nodeone_progresscount=0;
long     nodetwo_progresscount=0;

long     nodeone_savedcount=0;
long     nodetwo_savedcount=0;

int      bruteforceCount = 0;

double   pairCount = 0.0;

int      Nmin      = 1000;  // Minimum number of particles which a node can contain, should be > 1; previously 100.
int      sortDim, treelevel;
int      NLeft, NRight;

int      tree_N, firstSplitDim, tree_labelCount;
int      pairwisepdf_pairCount=0;
  
double   SplitValue;

double*  gg, *rr, *gr, *landy_xi, *dummy_gg, *rr_0, *rr_2, *rr_4, *rr_6, *rr_8, *rr_10, *xi0, *xi2, *xi4, *xi6, *xi8, *xi10, *logrbins, ngg, ngr, nrr;

double*  gg_meanr, *gg_meanmu;
double*  gr_meanr, *gr_meanmu;
double*  rr_meanr, *rr_meanmu, *rr_meann, *rr_meanwn;

double*  mean_xi;

int*     pairsperbin;

int      Distinct_pairCount = 0;

int      Sumof_Childrenparticles = 0;

int      NodeStart_Indices[1000];

double*  pairwisepdf_pairs;

double   max_pairSeparation = -99.0;
double   min_pairSeparation = +99.;

double   max_mu = -99.0;
double   min_mu = +99.0;

double*  meanMono;
double*  meanQuad;

double*  meanMono_error;
double*  meanQuad_error;

double*  Mono_suitable_mockCount;
double*  Quad_suitable_mockCount;


// ** Functions ** //
int    print_xi();
int    print_rr();
int    print_dd();
int    print_dr();

int    print_xiMultipoles();
double computeNorm(int Nfirst, int Nsecnd);

int    assignLeafValues(Particle cat[], double xCoors[], double yCoors[], double zCoors[], double disp[], int N);
int    bruteforce_nonodes(double* C0, double* C2, double* C4, double* C6, double* C8, double* C10, double* r, double* mu, Particle* cat, Particle* cat2, int N, int N2, int sameTree);

int                            postprocess_pairs(double* n, double* wn, double* r, double* mu, Node* node1, Node* node2);
int         bruteforceCountpairs_betweenChildren(double* n, double* wn, double* r, double* mu, Node* node1, Node* node2, int sameTree);
int    findSuitableNodePairs_bruteforcePairCount(double* n, double* wn, double* r, double* mu, Node* node1, Node* node2, int sameTree);

Node*  Create_toyChildNode();
Node*   buildTree(Particle cat[], int N);
Node*  createNode(Particle cat[], int N, int SplitDim, double xmin[NDIM], double xmax[NDIM]);

int    sortby_position_alongDim_splitDim(const void *a, const void *b);

double pair_zmu(Particle a, Particle b, double log10r);
double log10_particleSeparation(Particle a, Particle b);

double log10_minimum_modDisplacementBetweenNodes(Node *node1, Node *node2);
double log10_maximum_modDisplacementBetweenNodes(Node *node1, Node *node2);

int    free_tree(Node *t);
