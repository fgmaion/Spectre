//-- Cosmology -- //
double   H_0;
double   H_0inPerSec;
double   HubbleTime;

const double gamma_GR           =    0.545;
const double gamma_DGP          =  11./16.;

//-- Misc. --//
int          boxlabel, Index, local_Index;
int          xlabel, ylabel, zlabel;


// z - comoving distance data points for interpolation.
int                     nPoints = 1000;
double                   z_Array[1001];
double        ComovingDistance_z[1001];
double       redshiftInterval  = 0.015;

double       z_ComovingDistance_2derivatives[1001]; // z_chi_2d.
double       ComovingDistance_z_2derivatives[1001];

//-- VIPERS --//
int          fieldFlag;

double   TotalW1W4area;


//-- Functions --//
double   HubbleCnst(double z);
