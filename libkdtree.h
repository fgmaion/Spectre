#define NDIM 3

// logarithmic binning in r.
double   zerolog;
double    maxlog;
int     nlogbins;
double  logbinsz;                 

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

int   Nmin   =  100;      // Minimum number of particles which a Node can contain, should be > 1.

int   sortDim, treelevel;
int   NLeft, NRight;

int   tree_N, firstSplitDim, tree_labelCount;
int   pairwisepdf_pairCount=0;
  
double SplitValue;

double **gg, **rr, **gr, **landy_xi, **dummy_gg, **rr_0, **rr_2, **rr_4, *xi0, *xi2, *xi4, *logrbins, ngg, ngr, nrr;

double **gg_meanr, **gg_meanmu;
double **gr_meanr, **gr_meanmu;
double **rr_meanr, **rr_meanmu;

double** mean_xi;

int*   pairsperbin;

int Distinct_pairCount = 0;

int    Sumof_Childrenparticles = 0;

typedef struct Particle{
  // int    SplitDim;
  
  double x[NDIM];
  
  double disp;
  
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

int    CountPairs_rMu(double **C0, double **C2, double **C4, double **r, double **mu, Node* firstTree, Node* secndTree, int sameTree);

int    assignLeafValues(Particle cat[], double xCoors[], double yCoors[], double zCoors[], double disp[], int N);

double* pairwisepdf_pairs;

double max_pairSeparation = -99.0;
double min_pairSeparation = +99.;

double max_mu = -99.0;
double min_mu = +99.0;

Node*  buildTree(Particle cat[],  int N);

Node*  Create_toyChildNode();

Node   *createNode(Particle cat[], int N, int SplitDim, double xmin[NDIM], double xmax[NDIM]);

int    sortby_position_alongDim_splitDim(const void *a, const void *b);

int    findSuitableNodePairs_bruteforcePairCount(double **C0, double **C2,  double **C4, double **r, double **mu, Node *node1, Node *node2, int sameTree);

int    bruteforceCountpairs_betweenChildren(double **C0, double **C2, double **C4, double **r, double **mu, Node *node1, Node *node2, int sameTree);

int bruteforce_nonodes(double **C0, double **C2, double **C4, double **r, double **mu, Particle* cat, Particle* cat2, int N, int N2, int sameTree);

double pair_zmu(Particle a, Particle b);

double log10_particleSeparation(Particle a, Particle b);

double log10_minimum_modDisplacementBetweenNodes(Node *node1, Node *node2);

double log10_maximum_modDisplacementBetweenNodes(Node *node1, Node *node2);

int    free_tree(Node *t);

double computeNorm(int Nfirst, int Nsecnd);

double* meanMono;
double* meanQuad;

double* meanMono_error;
double* meanQuad_error;

double* Mono_suitable_mockCount;
double* Quad_suitable_mockCount;
