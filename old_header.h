// Available functions. 
int          ApplyFKPweights();

int          CalcCorrections();

int          randNGP();
int          CoordinateCalc();
int          NGPCalc();
int          PkCalc();

double       (*pt2Pk)(double k)                 	= NULL;
double       (*pt2Xi)(double)       				= NULL;

double       (*pt2RSD_k)(double, double, int)       = NULL;
double       (*pt2RSD_r)(double, double, int)       = NULL;

double       (*pt2nz)(double)                       = NULL;

double       volavg_invnbar;
double       volavg_redshift;

int          VIPERSbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int          Celestialbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int          projectVIPERSsystem();


// Artificial window fn's for cubic run
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

// padded window fn. calculation.

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

// Pointers to interpolation functions. 
float        (*pt2zChiIntegrand)(float);
double       (*pt2interp_comovingDistance)(double)         = &interp_comovingDistance;
double       (*pt2interp_inverseComovingDistance)(double)  = &interp_inverseComovingDistance;

double        lo_MBlim;
double        hi_MBlim;  
         
double        MinChi3;                                                      // h^-1 Mpc, Approximately, fig. 14, Guzzo et al.  2013
double        MaxChi3;                                                      // Redshift limited sample, 0.7 < z < 0.9
double        IntervalChi3;

double*      meanCellRedshift    = NULL;

double       TotalVolume         = 0.0;
double       TotalSurveyedVolume = 0.0;

double               fkpWeightedVolume   = 0.0;
double               fkpSqWeightsVolume  = 0.0;
double               fkpShotNoiseCorr    = 0.0;
double       analyticfkpWeightedVolume   = 0.0;

int          loopCount;

// Spline and splint holders. 
double       z_ComovingDistance_2derivatives[1001];
double       ComovingDistance_z_2derivatives[1001];

// Memory assignment for ZADE catalogue input. 
int         fieldFlag;

double*             ra  = NULL;
double*            dec  = NULL;
double*           zobs  = NULL;
double*           zcos  = NULL;
double*            M_B  = NULL;
double*          zflag  = NULL;
int*              type  = NULL;
int*          photoMask = NULL;


// Value added catalogue parameters. 
double*           zpec   = NULL;
double*          zphot   = NULL;

double*	         gal_z   = NULL;
double*            csr   = NULL;
double*       sampling   = NULL;  
double*     sampling35   = NULL;
char**       pointing    = NULL;
char**       quadrant    = NULL;
int*      flag_Nagoya    = NULL;
int*       flag_SSPOC    = NULL;
int*     flag_SSPOC35    = NULL;
double*       rand_sel   = NULL;
double*   fkp_galweight  = NULL;
double*   clip_galweight = NULL;


// derived parameters
bool*    Acceptanceflag = NULL;
double*     polarAngle  = NULL;
double*          rDist  = NULL;
double*          xCoor  = NULL;
double*          yCoor  = NULL;
double*          zCoor  = NULL;
double*          xVel   = NULL;
double*          yVel   = NULL;
double*          zVel   = NULL; 

// randoms.
double       alpha;

int          accepted_gals;
int          rand_number   = 0;
int          accepted_rand = 0;
double       fkp_accepted_rand = 0.0;

// TotalWeight is the sum of ZADE weight for the ZADE catalogue = Number of spec z + Number of used Photometric galaxies used (including compensation for sampling).
double       TotalZADEWeight   = 0.0;
double       MeanNumberDensity = 0.0;

int          boxlabel;
int          xlabel, ylabel, zlabel;

// FFTw calc parameters.
int          Index;
int          local_Index;
int          xNyquistIndex; 
int          yNyquistIndex;
int          zNyquistIndex;

double        kSq;
double        kIntervalx;
double        kIntervaly;
double        kIntervalz;

double        kmodulus;
double        mu;

double        kbinInterval;

double        k_x, k_y, k_z;

double        H_kReal;          
double        H_kImag;           

double        xNyquistWaveNumber;
double        yNyquistWaveNumber;
double        zNyquistWaveNumber;

// First column is mod k, second Pk.
double**     PkArray            = NULL;
double**     twodim_pk          = NULL;

double**     muIntervalPk       = NULL;

double*      legendre2weights   = NULL;
// Binned Pk parameters.

int          kbin_no;
int          Num_ModesInMuInterval = 0;
int*         modes_perbin        = NULL;

double       modkMax;

double*      del2               = NULL;
double*      mean_modk          = NULL;
double*      binnedPk           = NULL;
double*      logk_limits        = NULL;

// Randoms generation
int          lineNum      = 0;
int          NuRandoms    = 0; 
int          NuQuadrants  = 0;

// Jenkins scaling trick. 
double       Jenkins_foldfactor;

double*      sdltk               = NULL;
double*      sdltPk              = NULL;
double*      sdlt2d              = NULL;

// Camb linear P(k)
double*      lineark             = NULL;
double*      linearPk            = NULL;
double*      linear2d            = NULL;

// Factors for theoretical prediction of redshift space P(k).
double       kaiserFactor;
double       velDispersion;
double       f;
double       beta;
double       y;

double       fsigma8;
double       bsigma8;
double        alpha_pad;
double      epsilon_pad;

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
int**        zSpacemodesPerBin  = NULL;
double**     mean_perpk         = NULL;
double**     mean_losk          = NULL;
double**     d2_binnedpk     = NULL;

// Comoving number density calculation.
// double       zBinWidth;

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
double       TotalFKPweight;
double       Interim;

char         surveyType[200];

// clipping threshold

// rolling periodic cube. 
double       xcentre;
double       ycentre;
double       zcentre;

double       xroll;
double       yroll;
double       zroll;

double*      inputPk;

// VIPERS ra and dec of cell co-ordinates.

int CatalogueInput(char filepath[]);

int   CalcCellraDec();

double sdltNz(double z);

int   CatalogNumber;

double steradians2sqdegs(double inSteradians);

int    len;


double*  weights = NULL;

int   maxiter    = 20000000;


int lineNo;

int dof;

double sqdegs2steradians(double inSqdegs);

double  lightconeShot();
double  CubeShot();

int     muBinNumb;
int     modkBinNumb;

double* muBinLimits;

double  CentreRA;
double  CentreDec;

double  W1area;
double  W4area;


double  TotalObservedGalaxies    =    0.0;
double  dimmestAcceptedMagnitude =  -99.0;

double xtranslateDist;
double ytranslateDist;

double kaiserGauss_Monofactor(double ks, double beta);
double kaiserGauss_Quadfactor(double ks, double beta);
double kaiserGauss_Hexfactor(double ks, double beta);

double splintLinearPk(double k);

char      theoryPk_flag[200];
char      theoryRSD_flag[200];

double    shotNoiseCorrection_clipfactor;

// translation parameters for Stefano's co-ordinates. 
double    stefano_trans_x;
double    stefano_trans_y;
double    stefano_trans_z;

// mean sampling rate, e.g for FKP calc. 
double    meanSampling;

int min_zShift = 99, max_zShift = 0;
int min_yShift = 99, max_yShift = 0;
int min_xShift = 99, max_xShift = 0;

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

// float* HubbleConstantArray;
double* HubbleConstant2derivatives;

// Likelihood calculation.
double  app_sigma8; 

double  ChiSq_kmax;
double  ChiSq_kmin;

double    detCovariance;

// int*      ModeNumber;
gsl_matrix* Covariance;
gsl_matrix* sigma_norm;

gsl_vector* eval;
gsl_matrix* evec;

gsl_eigen_symmv_workspace* w;

gsl_vector* col;

double*   MeanMultipoles;

double    A11Sq;
double    linearBias;

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

int     besseltransform_order;
float   q0;

double Pk_powerlaw(double k, double r0, double gamma);

double* mono;
double* quad;
double*  hex;

double* monop;
double* quadp;
double*  hexp;

double Pk_powerlaw_truncated_xi(double r);

double HODPk_Gaussian(double k);

double LegendrePolynomials(double x, int n);

double u0, variance;

double splintConvQuad(double k);
double splintConvMono(double k);

double app_mean;

double* monoCorr;
double* monoCorr2d;
double* rmonoCorr;

double  AbelIntegrand(double);

double (*pt2funcforAbel)(double) = NULL;

double volAverage_Corrfn(double rmax);
double volNormedCorrfn(double r);

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

double AcceptedMax(double a[], bool b[], int n);
double AcceptedMin(double a[], bool b[], int n);

double chi_invertedStefanoBasis(double* xval, double* yval, double* zval);

double* expected_foldedRandCounts;
double calc_volavg_fkpweights();

double clipping_fSq;

double spherical_tophat(double k, double R);

double underlyingGaussian_sigma;

double GaussianFilter_radius;
double depletion_factor;

double (*pt2maskMultipoles)(double r, int transformOrder) = NULL;

int    Nz_fit(double output[]);
double model_NzGaussian(double z, double lnA, double z0, double sigma);

double nbar_dV(double chi);
double ChiSqEval_ap();






