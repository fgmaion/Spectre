int          randNGP();
int          NGPCalc();
int          PkCalc();

double       (*pt2Pk)(double k)                 	= NULL;
double       (*pt2Xi)(double)       				= NULL;

double       (*pt2RSD_k)(double, double, int)       = NULL;
double       (*pt2RSD_r)(double, double, int)       = NULL;

double       (*pt2nz)(double)                       = NULL;

int          VIPERSbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);
int          Celestialbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int          projectVIPERSsystem();

int          FullCube();
int          PencilBeamSurvey(int xlow, int xhi, int ylow, int yhi);

double       interp_nz(double x);

int          FirstColumnCompare(const void* a, const void* b);
int          SecondColumnCompare(const void* a, const void* b);

double       ArrayMax(double a[], int n);
double       ArrayMin(double a[], int n);

double       SumDoubleArray(double array[], int len);

int          printInterpPk();

int          assignbinninginterval(double interval);

int          zCubeCreate();
int          assignNGPMemory();
int          assignFFTMemory();
int          assignNGPMemory();
int          assignFFTMemory();
int          assign2DPkMemory();
int          assignConvolutionMemory();

int          freeRand();
int          freeHOD();
int          freeNGP();
int          freeFFTw();
int          freePk();
int          free_sdltInterp();
int          free_linear();
int          freeClipped();
int          free2dPk();

int          setInputPk();

double       interp_comovingDistance(double z);
double       interp_inverseComovingDistance(double r);

float        (*pt2zChiIntegrand)(float);
double       (*pt2interp_comovingDistance)(double)         = &interp_comovingDistance;
double       (*pt2interp_inverseComovingDistance)(double)  = &interp_inverseComovingDistance;

double        MinChi3;                                                      // h^-1 Mpc, Approximately, fig. 14, Guzzo et al.  2013
double        MaxChi3;                                                      // Redshift limited sample, 0.7 < z < 0.9
double        IntervalChi3;

double       TotalVolume         = 0.0;
double       TotalSurveyedVolume = 0.0;

int          loopCount;

double*     sampling35   = NULL;
int*      flag_Nagoya    = NULL;
int*       flag_SSPOC    = NULL;
int*     flag_SSPOC35    = NULL;
double*       rand_sel   = NULL;
char**       pointing    = NULL;
char**       quadrant    = NULL;
double*         csr;

// derived parameters
double*     polarAngle  = NULL;

// randoms.
double       fkp_accepted_rand = 0.0;

// FFTw calc parameters.
double        kSq;
double        kmodulus;
double        mu;

double        kbinInterval;

double        k_x, k_y, k_z;

// First column is mod k, second Pk.
double**     PkArray            = NULL;
double**     twodim_pk          = NULL;

double**     muIntervalPk       = NULL;


// Binned Pk parameters.
int          kbin_no;
int*         modes_perbin        = NULL;

double       modkMax;

double*      mean_modk          = NULL;
double*      binnedPk           = NULL;
double*      logk_limits        = NULL;

// Randoms generation
int          NuRandoms    = 0; 
int          NuQuadrants  = 0;

// Jenkins scaling trick. 
double       Jenkins_foldfactor;

// Factors for theoretical prediction of redshift space P(k).
double       f;
double       beta;
double       y;

double       Om_mz;

// Clipping
double*     PkCube              = NULL;

// double*     Corrfn              = NULL;
// double*     suppressedCorrfn    = NULL;
// double*     distortedCorrfn     = NULL;
// double*     clippedPk           = NULL;

// double*     W2_veck             = NULL;
// double*     W2_vecr             = NULL;
// double*     FFTW2_vecr_re       = NULL;
// double*     FFTW2_vecr_im       = NULL;

// Binning a 2D, redshift space P(k).
double**     d2_binnedpk     = NULL;

// Comoving number density calculation.
double*      zbins    = NULL;
double*      chibins  = NULL;
double*      Nchi     = NULL;
double*      nbar     = NULL;
double*      nbar_2d  = NULL;

double*      comovVol = NULL;

double chi_interval;
int    chibin_no;

// FKP weights. 
double       fkpPk;

// VIPERS ra and dec of cell co-ordinates.

int CatalogueInput(char filepath[]);

double sdltNz(double z);

int   CatalogNumber;

double steradians2sqdegs(double inSteradians);

int    len;


double*  weights = NULL;

int   maxiter    = 20000000;


int lineNo;

int dof;

double sqdegs2steradians(double inSqdegs);

int     muBinNumb;

double* muBinLimits;

double  CentreRA;
double  CentreDec;

double  W1area;
double  W4area;


double  TotalObservedGalaxies    =    0.0;
double  dimmestAcceptedMagnitude =  -99.0;

double splintLinearPk(double k);

char      theoryPk_flag[200];
char      theoryRSD_flag[200];

int polarPk_modeCount = 0;


// Calculation fof the age of the universe, leading to linear growth rate.  Change to double precision.
double*  AgeOftheUniverse;  // Units of H_0
double*  HubbleCnstWithTime;

// Calculation of the linear growth rate. 
void  (*pt2derivs)(float, float[], float[]) = NULL;
void  (*pt2rkqs)(float [], float [], int, float *, float, float, float [], float *, float *, void (*)(float, float [], float[])) = NULL;

int    nVar    = 2;

float* yStartArray;
float* yFinalArray;
float* y2derivsStartArray;

int   linearGrowth_nPoints;

float eps                = 0.00001;
float defaultStepSize    = 0.001;
float MinAllowedStepSize = 0.000001;

float InitialStartTime;
float FinalTime;
float Initial_lnScalefactor;
float Final_lnScalefactor;
    
// Defining declarations for these variables, with memory allocation xp[1..kmax] and yp[1..nvar][1..kmax]
// for the arrays, should be in the calling program.
double* Age2derivatives;
double* Om_mOfa;
double* f_Om_mOfa545;
double* f_Om_mOfa545_2derivs;
double* approx2linear_growthfactor;

double* lnAarray;
double* lnA2derivatives;
double* linear_growthfactor;
double* SplineParams_ofgrowthfactor;
double* SplineParams_ofgdot;

float* xVals;
float* length_scales;
float* gdot;
float* derivs_error;

float* Second_derivs;
float* Second_derivs_error;
float* Second_deriv_lengthscales;

float* SplineParams_ofgdotdot;

int    nDerivs;
double growthfactor_today;
double approx_growthfactor_today;

double* HubbleConstant2derivatives;

// Likelihood calculation.
double  app_sigma8; 

double  ChiSq_kmax;
double  ChiSq_kmin;

gsl_matrix* Covariance;
gsl_matrix* sigma_norm;

gsl_vector* eval;
gsl_matrix* evec;

gsl_eigen_symmv_workspace* w;

gsl_vector* col;

double*   MeanMultipoles;

double    ChiSqEval();

double    clippedVolume;

double*  eigenVals;
double** eigenVecs;
   
double  smallestEigenvalue;

int     chiSq_kmaxIndex;
int     chiSq_kminIndex;

// Index at which 
int     jenkins_foldIndex_foldedfile;
int     jenkins_foldIndex_unfoldedfile;
double  jenkins_fold_kjoin;

int     lineNo;

float   (*pt2_pkIntegrand)(float);

double Pk_powerlaw(double k, double r0, double gamma);

double* mono;
double* quad;
double*  hex;

double Pk_powerlaw_truncated_xi(double r);

double HODPk_Gaussian(double k);

double LegendrePolynomials(double x, int n);

double u0, variance;

double* monoCorr;
double* monoCorr2d;
double* rmonoCorr;

double  dra;
double  ddec;

int rabin;
int decbin;

int radecbinNumb = 300;

double rabinInterval;
double decbinInterval;

double*  ra_bins;
double* dec_bins;

double** binned_pairs;

double* fiber_flag;
double remaining;

double* rvals;
double*    xi;
double*  xi2D;

int splinexi_N =0;

int*    shuffle_rows;
double* rand_redshift;

double* xDisplacement;
double* yDisplacement;
double* zDisplacement;

int startint = -99;
int endint =-99;

double dummy;
double* kVals;

double** invA;

double zel_r;
double zel_xi;
double zel_err;
double zel_scale;

double splint_sigma2perp(double r);
double splint_sigma2para(double r);

double ukm_nfw_profile(double k);

double haloModel_pk(double k, double nbar, double beta, int Order);
double ukm2(double k);
double linearPk_Gaussian(double k);

double nfw_rvir;
double nfw_conc;

double* r_nfwinversion;
double* q_nfwinversion;
double* r2D_nfwinversion;

double rhobox, Mbox, Mhalo, Mpart;
double Delta_crit = 200.;
double Dplus      = 0.5;
double fGR        = 0.7;

double* hod_disp;

fftw_complex *outx, *outy, *outz;

fftw_plan    iplan_x, iplan_y, iplan_z;
fftw_plan    plan;

int*     shuffle_rows;

double anisoGauss_delta;
double anisoGauss_asigma;
double anisoGauss_bsigma;

int*   fftlog_indices;


// mass functions etc.
double sig2_R;

int   sigres;

// for splint of sigma^2 of R.
double*     sig2;
double*    radii;
double*  sig2_2D;

double*  lnm_jenkins;
double*  ln_invSig_jenkins;
double*  ln_invSig_jenkins_2D;

int      pk_lineNo;

double dummy;

double pk_hin;
double pk_lon;
double pk_hiA;
double pk_loA;

double modeCorrections(int k, int j, int i);

double* Rr;
double* Ir;
double* PkCube;

double chi_invertedStefanoBasis(double* xval, double* yval, double* zval);

double calc_volavg_fkpweights();

double spherical_tophat(double k, double R);

double underlyingGaussian_sigma;

double GaussianFilter_radius;
double depletion_factor;

int    Nz_fit(double output[]);
double model_NzGaussian(double z, double lnA, double z0, double sigma);

double nbar_dV(double chi);
double ChiSqEval_ap();






