// Available functions. 
int          comovDistReshiftCalc();

int          DisplacementCalc();

int          JenkinsCoordinates();
int          JenkinsFold(double original[], int lenArray, int axis);
int          ApplyJenkins();

int          EvaluateGridParameters();
double       SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval);

int          ApplyFKPweights();

int          CalcCorrections();

int          randNGP();

// double        zChi_Integrand(double x);
int          CoordinateCalc();
int          NGPCalc();
int          PkCalc();

double       splintMatterPk(double EvalPoint);
double       splintWindowfunc(double EvalPoint);
double       AnalyticSpherical(double k);
int          EvaluateConvolution(char type[]);

double 		 (*pt2AnisoWf)(double, double, double)  = NULL;
double       (*pt2Windowfn)(double)                 = NULL;
double       (*pt2Pk)(double)       				= NULL;
double       (*pt2Xi)(double)       				= NULL;

double       (*pt2RSD_k)(double, double, int)       = NULL;
double       (*pt2RSD_r)(double, double, int)       = NULL;

double       (*pt2nz)(double)                       = NULL;
double       (*pt2shot)(double)                     = NULL;

double       CubeShot();

// Artificial window fn's for cubic run
int          FullCube();
int          PencilBeamSurvey(int xlow, int xhi, int ylow, int yhi);

int          printWindowfuncSlices();
int          WindowfuncSlice(double kintervali, int ni, int x0, int y0, int z0, char filepath[]);

int          inputWindowfnBinNumb;

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
int          padwfPkCalc(int sidepad);
int          assignPaddedWindowfn(int padfactor);

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
int          freeConvolutionMemory();


// Convolution with anisotropic window fn.
int          AnisoGauss(double x, double y, double z);
int          AnisoConvolution();
int          setInputPk();
int          SetWfKernel();
double       ConvolveCell(int x, int y, int z);
int 		 convolve3DInputPk();

double       anisokGauss(double x, double y, double z);

int          MeasureAnisoWfKernel();
int          setMeasuredWfKernel();

double       interp_comovingDistance(double z);
double       interp_inverseComovingDistance(double r);

// Pointers to interpolation functions. 
float        (*pt2zChiIntegrand)(float);
double       (*pt2interp_comovingDistance)(double)         = &interp_comovingDistance;
double       (*pt2interp_inverseComovingDistance)(double)  = &interp_inverseComovingDistance;

char         vipersHOD_dir[200];

// directories
char         WindowfuncSlices_dir[200];


// filepaths
char         Windowfunc_xSlice[200];
char         Windowfunc_ySlice[200];
char         Windowfunc_zSlice[200];


// VIPERS HOD mock parameters
const double  Om_v      =   0.73;
const double  Om_r      =   0.0;
const double  Om_m      =   0.27;
const double  Om_b      =   0.0469;
const double  Om_tot    =   1.0;
const double  h         =   0.7;
const double  sigma_8   =   0.82;
const double  ns        =   0.95;

double        hz        =    0.0;

double        H_0;
double        H_0inPerSec;
double        HubbleTime;

// Selection parameters. Volume limited sample between redshift 0.7 and 0.9
double        redshiftHiLimit;
double        redshiftLowLimit;
double        absMagCut;

// Array to hold the coordinate limits of the VIPERS survey. 
double        AxisLimsArray[2][3];
         
double        CellSize;                                                     // Cell size, comoving distance, h^-1 Mpc

double        MinChi3;                                                      // h^-1 Mpc, Approximately, fig. 14, Guzzo et al.  2013
double        MaxChi3;                                                      // Redshift limited sample, 0.7 < z < 0.9
double        IntervalChi3;

double        CellVolume;

double*      densityArray        = NULL;
double*      FKPweights          = NULL;
double*      booldensity         = NULL;
double*      meanCellRedshift    = NULL;

double       TotalVolume         = 0.0;
double       TotalSurveyedVolume = 0.0;

double       fkpWeightedVolume   = 0.0;
double       fkpSqWeightsVolume  = 0.0;
double       SumOverStates_W2k   = 0.0;
double       fkpShotNoiseCorr    = 0.0;

// Dimensions of the padded volume (TotalVolume) in units of CellSize. 
int          n0, n1, n2;      
int          loopCount;

// z - Comoving distance data points for interpolation. 
int          nPoints = 1000;
double       z_Array[1001];
double       ComovingDistance_z[1001];
double       redshiftInterval  = 0.015;

// Spline and splint holders. 
double       z_ComovingDistance_2derivatives[1001];
double       ComovingDistance_z_2derivatives[1001];

// Memory assignment for ZADE catalogue input. 

int                ch;
int        Vipers_Num;

int*               id  = NULL;
double*             ra  = NULL;
double*            dec  = NULL;
double*           zobs  = NULL;
double*           zcos  = NULL;
double*            M_B  = NULL;
int*             type  = NULL;


// Value added catalogue parameters. 
double*           zpec  = NULL;
double*          zphot  = NULL;

double*	    zUtilized  = NULL;
double*            csr  = NULL;
double*       sampling  = NULL;  
double*     sampling35  = NULL;
char**       pointing  = NULL;
char**       quadrant  = NULL;
int*      flag_Nagoya  = NULL;
int*       flag_SSPOC  = NULL;
int*     flag_SSPOC35  = NULL;
double*       rand_sel  = NULL;


// derived parameters
bool*   Acceptanceflag = NULL;
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

int          rand_number   = 0;
int          accepted_rand = 0;

double*      rand_ra     = NULL;
double*      rand_dec    = NULL;
double*      rand_chi    = NULL; 
double*      rand_x      = NULL;
double*      rand_y      = NULL;
double*      rand_z      = NULL;
bool*        rand_accept = NULL;

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

float         fkmodulus;

double        kmodulus;
double        mu;

double        kbinInterval;

double        k_x, k_y, k_z;

double        WindowFunc;     
double        GaussianFilter;    

double        H_kReal;          
double        H_kImag;           

double        NyquistWaveNumber;

// First column is mod k, second Pk.
double**     PkArray            = NULL;
double**     TwoDpkArray        = NULL;

double**     muIntervalPk       = NULL;

double*      legendre2weights   = NULL;
// Binned Pk parameters.

int          kBinNumb;
int          LowerBinIndex;
int          UpperBinIndex;
int          Num_ModesInMuInterval = 0;
int*         modesPerBin        = NULL;

double       modkMax;

double*      del2               = NULL;
double*      midKBin            = NULL;
double*      meanKBin           = NULL;        
double*      binnedPk           = NULL;
double*      kBinLimits         = NULL;
double*      linearErrors       = NULL;

// Randoms generation
int          lineNum      = 0;
int          NuRandoms    = 0; 
int          NuQuadrants  = 0;

// Jenkins scaling trick. 
double       JenkinsScalefactor;

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

const double gamma_GR           =    0.545;
const double gamma_DGP          =  11./16.;

//const double gamma_fR         = 0.41-0.21*z;          

double       Om_mz;

// Clipping
double*     PkCube              = NULL;

double*     Corrfn              = NULL;
double*     suppressedCorrfn    = NULL;
double*     distortedCorrfn     = NULL;
double*     clippedPk           = NULL;

double*     W2_veck             = NULL;
double*     W2_vecr             = NULL;
double*     FFTW2_vecr_re       = NULL;
double*     FFTW2_vecr_im       = NULL;

// Binning a 2D, redshift space P(k).
int**        zSpacemodesPerBin  = NULL;
double**     mean_perpk         = NULL;
double**     mean_losk          = NULL;
double**     zSpaceBinnedPk     = NULL;

// Comoving number density calculation.
// double       zBinWidth;
double       VIPERS_SolidAngle;

double*      zbins   = NULL;
double*      chibins = NULL;
double*      Nchi      = NULL;
double*      nbar    = NULL;
double*       nbar_2d = NULL;

double*      comovVol = NULL;

double chi_interval;
int    chibin_no;

// double*        filteredComovingNumberDensity = NULL;
// double*       MeanSliceRedshift             = NULL;
// double*       LuminosityDistance            = NULL;
// double*       Schechter_fn                  = NULL;

// spline and splint comoving number density, n(z).

// FKP weights. 
double        fkpPk;
double        TotalFKPweight;
double        Chi;                       //  Comoving distance at redshift z for weight calculation.
double        Interim;

// Window fn. convolution, spherical symmetry.
int          InterpK_binNumber;
double*       kVals;                            // changed from float to double. 
float*       interpolatedPk;

int*         midKmodesperbin;
float*       midKBinNR;
float*       WindowFuncNR;
float*       WindowFunc2d;

double*      ConvolvedPk;
double*      ConvolvedPk2d;

float        inputPkEval         = 0.0;
float        WindowFuncEval      = 0.0;
float        PkEvalPoint         = 0.0;
float        muVal;
float        muInterval;

int          MuIntegralPrecision;

/*              */

double       Analytic2powerlaw(double k);
double       AnalyticGaussian(double k);

// Integral constraint correction.
float        WindowfnZeroPointEval;

// Nomenclature for applied survey window.
char         surveyType[200];

// padded window fn. calculation
int          sidepad;
int          padfactor;
int          padcellsNumber;
int          padIndex;

// clipping threshold
double       appliedClippingThreshold;

// rolling periodic cube. 
double       xcentre;
double       ycentre;
double       zcentre;

double       xroll;
double       yroll;
double       zroll;

// Anisotropic convolution
// int 		 xwfKernelsize;
// int          ywfKernelsize;
// int          zwfKernelsize;

float        fPkCubeEntry;

double 		 PkCubeEntry;
double 		 ConvNorm;

double*      inputPk;
double*		 windowFunc3D; 
double*		 convolvedPk3d;

double**     flattenedConvolvedPk3D;

double*      AnisoWfKernel;
int*         AnisoWfKernel_ModeNumb;


// VIPERS ra and dec of cell co-ordinates.
double*   Cell_rotatedXvals;
double*   Cell_rotatedYvals;
double*   Cell_rotatedZvals;

double*   Cell_raVIPERSsystem;
double*   Cell_decVIPERSsystem;
double*   Cell_chiVIPERSsystem;

double*   Cell_VIPERSweights;
double*   Cell_VIPERSbools;
double*   Cell_SurveyLimitsMask;
double*   Cell_AppliedWindowFn;

int VIPERSbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int Celestialbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len);

int projectVIPERSsystem();

double UpperChiLimit;
double LowerChiLimit;

double UpperRAlimit;
double LowerRAlimit;

double UpperDecLimit;
double LowerDecLimit;

int CatalogueInput(char filepath[]);

int xminKernelshift;
int xmaxKernelshift;

int yminKernelshift;
int ymaxKernelshift;

int zminKernelshift;
int zmaxKernelshift;

int   iishift, kkshift, jjshift;

double KernelMaxk;

int   filterIndex;
int   convIndex;
int   PkIndex;

int   CalcCellraDec();

double sdltNz(double z);

int   CatalogNumber;

double steradians2sqdegs(double inSteradians);

double ZeroPointNorm();


// Convolution calc. spherical symmetry. 

int          splineInputWindowfn(char filepath[]);

float        kConvScope;
float       muConvScope;
float        qConvScope;

double  NormKernel(double q);
float  fNormKernel(float);

double SphericalW2NormEval();

int EvaluateSphericalConv(char Convtype[]);

double SphConvPk(double kval);
float fqConvKernel(float q);
double qConvKernel(double q);
float fmuConvKernel(float q);
double muConvKernel(double mu);
double SphConvZeroPoint();
double ZeroPoint;
double ConvPkZeroPoint;

int    IntegralConstraintCorrection();

double ConvolutionNorm(double array[]);


// fitting n(z) using Powell's method.
int    len;
 
double (*pt2Theoryfn)(float array[], double z);

double*    xdata = NULL;
double*    ydata = NULL;
double*    edata = NULL;
double*  weights = NULL;

double*    zdata = NULL;
int*      Nzdata = NULL;

int   maxiter    = 20000000;

float returnval;


float **unitMatrix;

double HubbleCnst(double z);

int lineNo;

float* stepArray;
float* nz_startParams;
float* NewParams;

int mc_jump(float paramsArray[], float* minChi2, int paramNumber, float fracPrecision);

const gsl_rng_type* gsl_rng_T;
gsl_rng*            gsl_rng_r;

int dof;

int mc_loopCount;

double sqdegs2steradians(double inSqdegs);

// float* filteredNumberAtRedshift;
// float* filtered_divfuncln_Atz;

// float nzSigma;

//  Apodise the window fn. to supress the Gibb's phenomenon. 
double  GibbsSkinDepth;


double* Cell_SurveyEdge;
double* Cell_ApodiseWeights; 
double* Cell_ShortDist2edge;

double  apodisedVolume;

float   nbarChi2(float Chi);
float   Chi2(float Chi);
double  lightconeShot();
double  CubeShot();

double* loskBinLimits;
double* perpkBinLimits;

int     loskBinNumb;
int     perpkBinNumb;

double  perpkInterval;

int     forCount;

int     muBinNumb;
int     modkBinNumb;

double* muBinLimits;

int**    polar_modesPerBin;

double** mean_mu;
double** mean_modk;
double** polar2Dpk;
double** polar2DBinnedPk;

double*  kHexadecapole; 
double*  kQuadrupole;
double*  kMonopole;

double* kMonopole_expError;
double* kQuadrupole_expError;

double*  wfMonopole;
double*  wfQuadrupole;

float*  fkBinLimits;

float*  convMono;
float*  convQuad;

float*  convMono2d;
float*  convQuad2d;

float   TotalW1W4area;                                                                          

double  CentreRA;
double  CentreDec;

double  W1area;
double  W4area;

const gsl_rng_type* gsl_ran_T;
gsl_rng*            gsl_ran_r;

float* tophat;

double  TotalObservedGalaxies    =    0.0;
double  dimmestAcceptedMagnitude =  -99.0;

double xtranslateDist;
double ytranslateDist;

double kaiserGauss_Monofactor(double ks, double beta);
double kaiserGauss_Quadfactor(double ks, double beta);
double kaiserGauss_Hexfactor(double ks, double beta);

double splintLinearPk(double k);

int**  wfKernel_minAmpIndices;
int    largeAmpIndices = 0;
double minAmp_ConvolutionNorm(double array[]);

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

double wfKernel_minAmp;

double Wfzeropoint;
double convolution_modkmax;
int    subVol1, subVol2, subVol3, subVolCount;

int polarPk_modeCount = 0;


// Calculation fof the age of the universe, leading to linear growth rate. 
float*  AgeOftheUniverse;  // Units of H_0
float*  HubbleCnstWithTime;

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

float (*pt2f_Om_545)(float);
float (*pt2AgeIntegrand)(float);

float (*pt2linearGrowthRate)(float)         = NULL;
float (*pt2lnlinearGrowthRate)(float)       = NULL;
float (*pt2linearGrowthRate_deriv)(float)   = NULL;
float (*pt2linearGrowthRate_2deriv)(float)  = NULL;
    
    // Defining declarations for these variables, with memory allocation xp[1..kmax] and yp[1..nvar][1..kmax]
    // for the arrays, should be in the calling program.

// float* TimeArray;
float* Age2derivatives;
float* Om_mOfa;
float* f_Om_mOfa545;
float* f_Om_mOfa545_2derivs;
float* approx2linear_growthfactor;

float* lnAarray;
float* lnA2derivatives;
float* linear_growthfactor;
float* SplineParams_ofgrowthfactor;
float* SplineParams_ofgdot;
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
float* HubbleConstant2derivatives;

// Likelihood calculation.
double**     Multipoles;
double**    dMultipoles;

double*       sigmaNorm;

// flattened covariance matrix. 
double* flatCov;

double  ChiSq_kmax;
double  ChiSq_kmin;

double* betaPosterior;
double* sigmaPosterior;
double** betaSigmaPosterior;

double  PosteriorNorm = 0.0;


double    detCovariance;

// int*      ModeNumber;
gsl_matrix* Covariance;
gsl_matrix* sigma_norm;

gsl_vector* eval;
gsl_matrix* evec;

gsl_eigen_symmv_workspace* w;

gsl_vector* col;

// double**  Covariance;
// double**  Covariance_CofactorMatrix;

double*   MeanMultipoles;

double*   mvGauss;
// double**  invCov;

double    A11Sq;
double    linearBias;
double    ChiSqEval();

int       Res;
double    dRes;

double*** ChiSqGrid;
double*** lnLikelihoodGrid;
double    minChiSq;

double    minChiSq_beta;
double    minChiSq_A11Sq;
double    minChiSq_sigma;

double    paramNumber;

double    min_beta;
double    max_beta;

double    min_velDisperse;
double    max_velDisperse;

double    min_A11Sq;
double    max_A11Sq;

double    min_linearbias;
double    max_linearbias;

double    clippedVolume;

// Boot strap catalogue generation.
int*      BootStrap_flag;
double*   BootStrap_Wght;

double*   xBoot;
double*   yBoot;
double*   zBoot;
double*   deltaBoot;

double*  eigenVals;
double** eigenVecs;
    
// double** sigmaNorm;
int order, nrotations;

double* xdata;
double* xtheory;
    
double* ydata;
double* ytheory;

double  smallestEigenvalue;
int     chiSq_kmaxIndex;

int     hiMultipoleOrder;
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

double*  cmono;
double* cmonop;

float* rVals;
float* fcmono;
float* fcmono2D;

int    FFTlogRes;

// double* fftlogr;
// double* fftlogk;

double Pk_powerlaw_truncated_xi(double r);

double HODPk_Gaussian(double k);

double*  theta;
double** mixingmatrix;

double**  rmodulus_vec;
double*   kmodulus_vec;
double*         mu_vec;

double LegendrePolynomials(double x, int n);

double unitTheory(double r, double k, double d, int order);

double u0, variance;

double splintConvQuad(double k);
double splintConvMono(double k);

double apparent_mean;

double* windfn_rvals;
double* windfn_rMonopole;
double* windfn_rQuadrupole;
double* windfn_rHexadecapole;

double* windfn_rMonopole2d;
double* windfn_rQuadrupole2d;
double* windfn_rHexadecapole2d;

int     rbinNumb;

double* monoCorr;
double* monoCorr2d;
double* rmonoCorr;

double  AbelIntegrand(double);

double (*pt2funcforAbel)(double) = NULL;

double Abel_rp;

double volAvgOfXi;

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

double* sigma2rvals;

double* sigma2perp;
double* sigma2para;

double* sigma2perp_2d;
double* sigma2para_2d;

double zel_vecr[3];
double zel_vecx[3];

double zel_qmin[3];
double zel_qmax[3];

double dummy;

double** invA;

double zel_r;
double zel_xi;
double zel_err;
double zel_scale;

double splint_sigma2perp(double r);
double splint_sigma2para(double r);

double NFW_densityprofile_NormedFourier(double k);

double haloModel_pk(double k, double beta, int Order);
double ukm2(double k);
double linearPk_Gaussian(double k);

double NFW_rvir;
double NFW_conc;
int    haloNumber;

double* r_Nfwinversion;
double* q_Nfwinversion;
double* rq2D_Nfwinversion;

double rhobox, Mbox, Mhalo, Mpart;
double Delta_crit = 200.;
double Dplus = 0.5;
double fGR   = 0.7;

double* hod_disp;

fftw_complex *outx, *outy, *outz;

fftw_plan    iplan_x, iplan_y, iplan_z;

int*     shuffle_rows;

double anisoGauss_delta;
double anisoGauss_asigma;
double anisoGauss_bsigma;


int     VIPERS_mask_lineNo_lo;

double* VIPERS_maskr_lo;

double* VIPERS_maskMono_lo;
double* VIPERS_maskQuad_lo;
double* VIPERS_maskHex_lo;

double* VIPERS_maskMono2D_lo;
double* VIPERS_maskQuad2D_lo;
double* VIPERS_maskHex2D_lo;

double mask_monopolenorm_lo;

int     VIPERS_mask_lineNo_hi;
double* VIPERS_maskr_hi;
double* VIPERS_maskMono_hi;
double* VIPERS_maskQuad_hi;
double* VIPERS_maskHex_hi;

double* VIPERS_maskMono2D_hi;
double* VIPERS_maskQuad2D_hi;
double* VIPERS_maskHex2D_hi;

double mask_monopolenorm_hi;

double loRes_highRes_join;

double splint_VIPERS_maskMono(double r);
double splint_VIPERS_maskQuad(double r);
double splint_VIPERS_maskHex(double r);

double splint_VIPERS_maskMultipoles(double r, int transformOrder);

int VIPERS_kSpace_multipoles_lineNo;

double* VIPERS_k;
double* VIPERS_kMono;
double* VIPERS_kMono2D;
double* VIPERS_kQuad;
double* VIPERS_kQuad2D;

double splint_VIPERS_kSpaceMono(double k);
double splint_VIPERS_kSpaceQuad(double k);

int*   fftlog_indices;
