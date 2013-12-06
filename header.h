// Available functions. 
int          comovDistReshiftCalc();

int          JenkinsCoordinates();
int          JenkinsFold(float original[], int lenArray, int axis);
int          ApplyJenkins();

int          EvaluateGridParameters();
float        SolidAngleCalc(float decLowerBound, float decUpperBound, float raInterval);

int          randGenerate();
int          pointRandGen(int lineNum, float Min_SinDec, float Max_SinDec, float Min_rA, float Max_rA);
int          RandCoorCalc();
int          randNGP();

int          ApplyFKPweights();

int          CalcCorrections();

int          randNGP();

float        Integrand(float x);
int          CoordinateCalc();
int          NGPCalc();
int          PkCalc();

float        splintMatterPk(float EvalPoint);
float        splintWindowfunc(float EvalPoint);
float        AnalyticSpherical(float k);
int          EvaluateConvolution(char type[]);

float        EvaluatefilterNormalisation();
float        NormalisationKernel(float q);

float 		 (*pt2AnisoWf)(float, float, float) = NULL;
float        (*pt2Windowfn)(float) 				= NULL;
float        (*pt2Pk)(float)       				= NULL;

float        SphericalWindowfuncEval(float k, float R);

// Artificial window fn's for cubic run
int          FullCube();
int          PencilBeamSurvey(int xlow, int xhi, int ylow, int yhi);

int          printWindowfuncSlices();
int          WindowfuncSlice(float kintervali, int ni, int x0, int y0, int z0, char filepath[]);

int          inputWindowfnBinNumb;

float        interp_nz(float x);
float        testWindowfn(float k2);

int          FirstColumnCompare(const void* a, const void* b);
int          SecondColumnCompare(const void* a, const void* b);

float        arrayMax(float a[], int n);
float        arrayMin(float a[], int n);

double       DoubleArrayMax(double a[], int n);
double       DoubleArrayMin(double a[], int n);
double       SumDoubleArray(double array[], int len);
float        SumFloatArray(float array[], int len);

int          ConvolveTheory();
int          splineInputWindowfn(char filepath[]);

float         kConvScope;
float        muConvScope;
float         qConvScope;

float        ConvolvedPkQrombCalc(float kval);
float        muConvKernel(float mu);

float        qConvKernel(float q);

float        ConvolvedPkCalculation(float kval);

int          IntegralConstraintCorrection();

int          printInterpPk();

int          assignbinninginterval(float interval);

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
float        anisoGauss(float x, float y, float z);
int          AnisoConvolution();
int          setInputPk();
int          SetWfKernel();
float        ConvolveCell(int x, int y, int z);
float 		 filter3Dnorm();
int 		 convolve3DInputPk(float convolvedPk[], float input[]);

int          EstimateAnisoWfKernel();
int          setMeasuredWfKernel();

// Pointers to interpolation functions. 
float        (*pt2Integrand)(float)                         = &Integrand;
float        (*pt2interp_comovingDistance)(float)           = &interp_comovingDistance;
float        (*pt2interp_inverseComovingDistance)(float)    = &interp_inverseComovingDistance;

char         vipersHOD_dir[200];

// directories
char         WindowfuncSlices_dir[200];


// filepaths
char         Windowfunc_xSlice[200];
char         Windowfunc_ySlice[200];
char         Windowfunc_zSlice[200];


// VIPERS HOD mock parameters
const float  Om_v      =   0.73;
const float  Om_r      =   0.0;
const float  Om_m      =   0.27;
const float  Om_b      =   0.0469;
const float  Om_tot    =   1.0;
const float  h         =   0.7;
const float  sigma_8   =   0.82;
const float  ns        =   0.95;

float        hz        =    0.0;
// Selection parameters. Volume limited sample between redshift 0.7 and 0.9
float        redshiftHiLimit;
float        redshiftLowLimit;
float        absMagCut;

// Array to hold the coordinate limits of the VIPERS survey. 
float        AxisLimsArray[2][3];
         
float        CellSize;                                                     // Cell size, comoving distance, h^-1 Mpc

float        MinChi3;                                                      // h^-1 Mpc, Approximately, fig. 14, Guzzo et al.  2013
float        MaxChi3;                                                      // Redshift limited sample, 0.7 < z < 0.9
float        IntervalChi3;

float        CellVolume;

double*      densityArray        = NULL;
double*      FKPweights          = NULL;
double*      booldensity         = NULL;

float        SumOfVIPERSbools    = 0.0;

float        TotalVolume         = 0.0;
float        TotalSurveyedVolume = 0.0;

double       fkpWeightedVolume   = 0.0;
double       fkpSqWeightsVolume  = 0.0;
double       fkpShotNoiseCorr    = 0.0;

// Dimensions of the padded volume (TotalVolume) in units of CellSize. 
int          n0, n1, n2;      
int          loopCount;

// Declaration  for Random Generation. 

float        Min_rA;  
float        Max_rA;
float        Min_SinDec;
float        Max_SinDec;

long         randCall_chi3    = 45;
long         randCall_sinDec  = 95;
long         randCall_rA      = 153;

float*       rand_chi         = NULL;
float*       rand_dec         = NULL;
float*       rand_rA          = NULL;
float*       rand_rshift      = NULL;
float*       rand_polar       = NULL;
float*       rand_x           = NULL;
float*       rand_y           = NULL;
float*       rand_z           = NULL;

// z - Comoving distance data points for interpolation. 
int          nPoints = 1000;
float        z_Array[1001];
float        ComovingDistance_z[1001];
float        redshiftInterval  = 0.015;

// Spline and splint holders. 
float        z_ComovingDistance_2derivatives[1001];
float        ComovingDistance_z_2derivatives[1001];

// Memory assignment for ZADE catalogue input. 

int                ch;
int        Vipers_Num;

int*               id  = NULL;
float*             ra  = NULL;
float*            dec  = NULL;
float*           zobs  = NULL;
float*           zcos  = NULL;
float*            M_B  = NULL;
int*             type  = NULL;


// Value added catalogue parameters. 
float*           zpec  = NULL;
float*          zphot  = NULL;

float*	    zUtilized  = NULL;
float*            csr  = NULL;
float*       sampling  = NULL;  
float*     sampling35  = NULL;
char**       pointing  = NULL;
char**       quadrant  = NULL;
int*      flag_Nagoya  = NULL;
int*       flag_SSPOC  = NULL;
int*     flag_SSPOC35  = NULL;
float*       rand_sel  = NULL;


// derived parameters
float*     polarAngle  = NULL;
float*          rDist  = NULL;
float*          xCoor  = NULL;
float*          yCoor  = NULL;
float*          zCoor  = NULL;
float*          xVel   = NULL;
float*          yVel   = NULL;
float*          zVel   = NULL; 

// TotalWeight is the sum of ZADE weight for the ZADE catalogue = Number of spec z + Number of used Photometric galaxies used (including compensation for sampling).
float        TotalZADEWeight   = 0.0;
float        MeanNumberDensity = 0.0;

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

double        WindowFunc;     
double        GaussianFilter;    

double        H_kReal;          
double        H_kImag;           

double        NyquistWaveNumber;

// First column is mod k, second Pk.
double**     PkArray            = NULL;
double**     TwoDpkArray        = NULL;

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

float*      sdltk               = NULL;
float*      sdltPk              = NULL;
float*      sdlt2d              = NULL;

// Camb linear P(k)

float*      lineark             = NULL;
float*      linearPk            = NULL;
float*      linear2d            = NULL;

// Factors for theoretical prediction of redshift space P(k).
double       kaiserFactor;
double       velocityDispersion;
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

// Binning a 2D, redshift space P(k).
int**        zSpacemodesPerBin  = NULL;
double**     mean_perpk         = NULL;
double**     mean_losk          = NULL;
double**     zSpaceBinnedPk     = NULL;

// Comoving number density calculation
int          zBinNumber;

float        zBinWidth;
float        VIPERS_SolidAngle;

float*       redshiftSlices     = NULL;
float*       ChiSlices          = NULL;
float*       NumberAtReshift    = NULL;
float*       ComovingVolumeAtZ  = NULL;
float*       ComovingNumberDensity  = NULL;


// spline and splint comoving number density, n(z).
float*       ComovingNumberDensity2d    = NULL;

// 
float        linearBias;
float        A11;

// FKP weights. 
float        fkpPk;
float        TotalFKPweight;
float        Chi;                       //  Comoving distance at redshift z for weight calculation.
float        Interim;


// Convolving theory with window fn.
int          InterpK_binNumber;
float*       kVals;
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
float        FilterNormalisation = 0.0;
float        muVal;
float        muInterval;

int          MuIntegralPrecision;

float        Analytic2powerlaw(float k);
float        AnalyticGaussian(float k);

// Integral constraint correction.
double       ConvolvedPkZeroPoint;
float        WindowfnZeroPointEval;

// Nomenclature for applied survey window.
char         surveyType[200];


// padded window fn. calculation
int          sidepad;
int          padfactor;
int          padcellsNumber;
int          padIndex;

// clipping threshold
float        appliedClippingThreshold;

// rolling periodic cube. 
double       xcentre;
double       ycentre;
double       zcentre;

double       xroll;
double       yroll;
double       zroll;

// Anisotropic convolution
int 		 wfKernelsize;

float 		 PkCubeEntry;
float 		 wf3Dnorm;

float*       inputPk;
float*		 windowFunc3D; 
float*		 convolvedPk3d;

double**     flattenedConvolvedPk3D;

float*       AnisoWfKernel;
int*         AnisoWfKernel_ModeNumb;


// VIPERS ra and dec of cell co-ordinates.
float* 	 Cell_rotatedXvals;
float* 	 Cell_rotatedYvals;
float* 	 Cell_rotatedZvals;

float*   Cell_raVIPERSsystem;
float*   Cell_decVIPERSsystem;
float*   Cell_chiVIPERSsystem;

float*   Cell_VIPERSweights;
float*   Cell_VIPERSbools;
float*   Cell_SurveyLimitsMask;
float*   Cell_AppliedWindowFn;

int VIPERSbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len);

int Celestialbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len);

int projectVIPERSsystem();

float UpperChiLimit;
float LowerChiLimit;

float UpperRAlimit;
float LowerRAlimit;

float UpperDecLimit;
float LowerDecLimit;

int CatalogueInput(char filepath[]);