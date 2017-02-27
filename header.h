//-- Cosmology -- //
double   H_0;
double   H_0inPerSec;
double   HubbleTime;

const double gamma_GR           =    0.545;
const double gamma_DGP          =  11./16.;

//-- Misc. --//
double       Interim, dummy;
int          Index, local_Index;


// z - comoving distance data points for interpolation.
int                     nPoints = 1000;
double                   z_Array[1001];
double        ComovingDistance_z[1001];
double       redshiftInterval  = 0.015;

double       z_ComovingDistance_2derivatives[1001]; // z_chi_2d.
double       ComovingDistance_z_2derivatives[1001];

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


//-- Functions --//
double   HubbleCnst(double z);


double AcceptedMax(double a[], bool b[], int n);
double AcceptedMin(double a[], bool b[], int n);
