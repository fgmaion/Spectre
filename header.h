//-- Cosmology --//
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

//-- VIPERS --//
double   TotalW1W4area;


double   HubbleCnst(double z);
