//-- Misc. --//
double       Interim, dummy;
int          Index, local_Index, loopCount, lineNo;

//-- Cosmology -- //
double   H_0;
double   H_0inPerSec;
double   HubbleTime;

const double gamma_GR           =    0.545;
const double gamma_DGP          =  11./16.;

// z - comoving distance data points for interpolation.
int                     nPoints = 1000;
double                   z_Array[1001];
double        ComovingDistance_z[1001];
double       redshiftInterval  = 0.015;

double       z_ComovingDistance_2derivatives[1001]; // z_chi_2d.
double       ComovingDistance_z_2derivatives[1001];

// Calculation of the age of the universe, leading to linear growth rate.  Change to double precision.
double*  AgeOftheUniverse;  // Units of H_0
double*  HubbleCnstWithTime;

// Calculation of the linear growth rate.
void  (*pt2derivs)(float, float[], float[]) = NULL; // FLOAT?!
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


//-- spline/splint real-space P(k) model --//
int      pk_lineNo;

double*      sdltk               = NULL; // hod p(k).
double*      sdltPk              = NULL;
double*      sdlt2d              = NULL;

double*      lineark             = NULL; // Camb linear p(k)
double*      linearPk            = NULL;
double*      linear2d            = NULL;

double pk_hin; // power law index and amplitude for k<<1 and k>>1 extrapolations of real-space p(k).
double pk_lon;
double pk_hiA;
double pk_loA;

//-- VIPERS --//
int          fieldFlag;
double   TotalW1W4area;

double  W1area;
double  W4area;

double  CentreRA;
double  CentreDec;

int   CatalogNumber;


// -- Clipping -- //
double u0, variance;


//-- Functions --//
double LegendrePolynomials(double x, int n);

double HubbleCnst(double z);

float  (*pt2zChiIntegrand)(float);  // FLOAT?!

double interp_nz(double);

double interp_comovingDistance(double z);
double interp_inverseComovingDistance(double r);

double AcceptedMax(double a[], bool b[], int n);
double AcceptedMin(double a[], bool b[], int n);
