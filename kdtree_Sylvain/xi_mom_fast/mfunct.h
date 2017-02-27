#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

int nint(double x);
double ran1(long *idum);
double gasdev(long *idum);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
double trapzd(double (*func)(double), double a, double b, int n);
double qtrap(double (*func)(double), double a, double b);
double trapzd2(double (*func)(double,double), double a, double b, int n, double c);
double qtrap2(double (*func)(double,double), double a, double b, double c);
double midpnt(double (*func)(double), double a, double b, int n);
double qromo(double (*func)(double), double a, double b,
	     double (*choose)(double(*)(double), double, double, int));
double midpnt2(double (*func)(double,double), double a, double b, int n, double amc);
double qromo2(double (*func)(double,double), double a, double b,
	      double (*choose)(double(*)(double,double), double, double, int, double),double amc);
double midpnt3(double (*func)(double,double,double), double a, double b, int n, 
	       double amc, double rlog);
double qromo3(double (*func)(double,double,double), double a, double b,
	      double (*choose)(double(*)(double,double,double), double, double, int, 
	      double, double),double amc,double rlog);

double gammln(double xx);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
double gammq(double a, double x);
void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
	 double *b, double *siga, double *sigb, double *chi2, double *q);

int ct_lines(char *file);
int ct_cols(char *file);
