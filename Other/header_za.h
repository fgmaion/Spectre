double zel_r;
double zel_xi;
double zel_err;
double zel_scale;

double* sigma2rvals;

double* sigma2perp;
double* sigma2para;

double* sigma2perp_2d;
double* sigma2para_2d;

double zel_vecr[3];
double zel_vecx[3];

double zel_qmin[3];
double zel_qmax[3];


// -- Functions -- //
int   DisplacementCalc();

double splint_sigma2perp(double r);
double splint_sigma2para(double r);
