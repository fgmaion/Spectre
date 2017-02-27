#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

//#define free(aa) {printf("[%s][ligne %d] Liberation bloc %s à %p\n",__FILE__,__LINE__,#aa,aa);free(aa);}

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define NR_END 1
#define FREE_ARG char*

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5
#define FUNC(x) ((*func)(x))

#define TEPS 1.0e-5
#define TJMAX 20
#define TFUNC(x,c) ((*func)(x,c))

#define ITMAXB 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20

#define ITMAX 200
#define TOL 2.0e-4

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define GITMAX 100
#define GEPS 3.0e-7
#define FPMIN 1.0e-30

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *) calloc((size_t)(nh-nl+1+NR_END),sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *) calloc((size_t)(nh-nl+1+NR_END),sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **)malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) calloc((size_t)(nrow+NR_END),sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) calloc((size_t)(nrow*ncol+NR_END),sizeof(double));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **)malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **)malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **)malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***)malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***)malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **)malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *)malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double f3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}

/***! Calculate the nearest integer of a double !**/

int nint(double x)
{
  if ((floor(x)-x)*(floor(x)-x) <= (ceil(x)-x)*(ceil(x)-x)) return (int)floor(x);
  else return (int)ceil(x);
}

/***! Random number generator (from numerical recipies) !***/

double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

/***! Normal deviate (from numerical recipies) !***/

double gasdev(long *idum)
{
  double ran1(long *idum);
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

/***! Polynomial interpolation !***/

void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
        int i,m,ns=1;
        double den,dif,dift,ho,hp,w;
        double *c,*d;

        dif=fabs(x-xa[1]);
        c=dvector(1,n);
        d=dvector(1,n);
        for (i=1;i<=n;i++) {
                if ( (dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        *y=ya[ns--];
        for (m=1;m<n;m++) {
                for (i=1;i<=n-m;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+m]-x;
                        w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
        }
        free_dvector(d,1,n);
        free_dvector(c,1,n);
}

/***! Spline interpolation !***/

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
    int i,k;
    double p,qn,sig,un,*u;

    u=dvector(1,n-1);
    if (yp1 > 0.99e30)
        y2[1]=u[1]=0.0;
    else {
        y2[1] = -0.5;
        u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
    }
    for (i=2;i<=n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
    }
    y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
    for (k=n-1;k>=1;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    free_dvector(u,1,n-1);
}

void splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
    int klo,khi,k;
    double h,b,a;

    klo=1;
    khi=n;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) printf("Bad xa input to routine splint\n");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

/***! Integration !***/

double trapzd(double (*func)(double), double a, double b, int n)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double qtrap(double (*func)(double), double a, double b)
{
        double trapzd(double (*func)(double), double a, double b, int n);
        void nrerror(char error_text[]);
        int j;
        double s,olds;

        olds = -1.0e30;
        for (j=1;j<=TJMAX;j++) {
                s=trapzd(func,a,b,j);
                if (fabs(s-olds) < TEPS*fabs(olds)) return s;
                olds=s;
        }
        nrerror("Too many steps in routine qtrap");
        return 0.0;
}

double trapzd2(double (*func)(double,double), double a, double b, int n, double c)
{
        double x,tnm,sum,del;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*(TFUNC(a,c)+TFUNC(b,c)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += TFUNC(x,c);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

double qtrap2(double (*func)(double,double), double a, double b, double c)
{
        double trapzd2(double (*func)(double,double), double a, double b, int n, double c);
        void nrerror(char error_text[]);
        int j;
        double s,olds;

        olds = -1.0e30;
        for (j=1;j<=TJMAX;j++) {
                s=trapzd2(func,a,b,j,c);
                if (fabs(s-olds) < TEPS*fabs(olds)) return s;
                olds=s;
        }
        nrerror("Too many steps in routine qtrap");
        return 0.0;
}

double midpnt(double (*func)(double), double a, double b, int n)
{
        double x,tnm,sum,del,ddel;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=(b-a)*FUNC(0.5*(a+b)));
        } else {
                for(it=1,j=1;j<n-1;j++) it *= 3;
                tnm=it;
                del=(b-a)/(3.0*tnm);
                ddel=del+del;
                x=a+0.5*del;
                sum=0.0;
                for (j=1;j<=it;j++) {
                        sum += FUNC(x);
                        x += ddel;
                        sum += FUNC(x);
                        x += del;
                }
                s=(s+(b-a)*sum/tnm)/3.0;
                return s;
        }
}

double qromo(double (*func)(double), double a, double b,
        double (*choose)(double(*)(double), double, double, int))
{
        int j;
        double ss,dss,h[JMAXP+1],s[JMAXP+1];

        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=(*choose)(func,a,b,j);
                if (j >= K) {
                        polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (fabs(dss) < 1.0e-5*fabs(ss)) return ss;
                }
                s[j+1]=s[j];
                h[j+1]=h[j]/9.0;
        }
        nrerror("Too many steps in routing qromo");
        return 0.0;
}

double midpnt2(double (*func)(double,double), double a, double b, int n, double amc)
{
        double x,tnm,sum,del,ddel;
        static double s;
        int it,j;

        if (n == 1) {
                return (s=(b-a)*((*func)(0.5*(a+b),amc)));
        } else {
                for(it=1,j=1;j<n-1;j++) it *= 3;
                tnm=it;
                del=(b-a)/(3.0*tnm);
                ddel=del+del;
                x=a+0.5*del;
                sum=0.0;
                for (j=1;j<=it;j++) {
                        sum += ((*func)(x,amc));
                        x += ddel;
                        sum += ((*func)(x,amc));
                        x += del;
                }
                s=(s+(b-a)*sum/tnm)/3.0;
                return s;
        }
}

double qromo2(double (*func)(double,double), double a, double b,
        double (*choose)(double(*)(double,double), double, double, int, double),double amc)
{
        int j;
        double ss,dss,h[JMAXP+1],s[JMAXP+1];

        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=(*choose)(func,a,b,j,amc);
                if (j >= K) {
                        polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (fabs(dss) < 1.0e-4*fabs(ss)) return ss;
                }
                s[j+1]=s[j];
                h[j+1]=h[j]/9.0;
        }
        nrerror("Too many steps in routing qromo2");
        return 0.0;
}

double midpnt3(double (*func)(double,double,double), double a, double b, int n, double amc, 
	       double rlog)
{
    double x,tnm,sum,del,ddel;
    static double s;
    int it,j;
    
    if (n == 1) {
	return (s=(b-a)*((*func)(0.5*(a+b),amc,rlog)));
    } else {
	for(it=1,j=1;j<n-1;j++) it *= 3;
	tnm=it;
	del=(b-a)/(3.0*tnm);
	ddel=del+del;
	x=a+0.5*del;
	sum=0.0;
	for (j=1;j<=it;j++) {
	    sum += ((*func)(x,amc,rlog));
	    x += ddel;
	    sum += ((*func)(x,amc,rlog));
	    x += del;
	}
	s=(s+(b-a)*sum/tnm)/3.0;
	return s;
    }
}

double qromo3(double (*func)(double,double,double), double a, double b,
	      double (*choose)(double(*)(double,double,double), double, double, int, 
	      double, double),double amc,double rlog)
{
    int j;
    double ss,dss,h[JMAXP+1],s[JMAXP+1];
    
    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
	s[j]=(*choose)(func,a,b,j,amc,rlog);
	if (j >= K) {
	    polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
	    if (fabs(dss) < 1.0e-4*fabs(ss)) return ss;
	}
	s[j+1]=s[j];
	h[j+1]=h[j]/9.0;
    }
    printf("Too many steps in routing qromo3\n");
    return 0.0;
}

/***! Linear fit !***/

/***! log of Gamma function !***/

double gammln(double xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}

void gcf(double *gammcf, double a, double x, double *gln)
{
        double gammln(double xx);
        int i;
        double an,b,c,d,del,h;

        *gln=gammln(a);
        b=x+1.0-a;
        c=1.0/FPMIN;
        d=1.0/b;
        h=d;
        for (i=1;i<=GITMAX;i++) {
                an = -i*(i-a);
                b += 2.0;
                d=an*d+b;
                if (fabs(d) < FPMIN) d=FPMIN;
                c=b+an/c;
                if (fabs(c) < FPMIN) c=FPMIN;
                d=1.0/d;
                del=d*c;
                h *= del;
                if (fabs(del-1.0) < GEPS) break;
        }
        if (i > GITMAX) printf("a too large, GITMAX too small in gcf\n");
        *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

void gser(double *gamser, double a, double x, double *gln)
{
        double gammln(double xx);
        int n;
        double sum,del,ap;

        *gln=gammln(a);
        if (x <= 0.0) {
                if (x < 0.0) printf("x less than 0 in routine gser\n");
                *gamser=0.0;
                return;
        } else {
                ap=a;
                del=sum=1.0/a;
                for (n=1;n<=GITMAX;n++) {
                        ++ap;
                        del *= x/ap;
                        sum += del;
                        if (fabs(del) < fabs(sum)*GEPS) {
                                *gamser=sum*exp(-x+a*log(x)-(*gln));
                                return;
                        }
                }
                //printf("a too large, GITMAX too small in routine gser\n");
                return;
        }
}

double gammq(double a, double x)
{
        void gcf(double *gammcf, double a, double x, double *gln);
        void gser(double *gamser, double a, double x, double *gln);
        double gamser,gammcf,gln;

        if (x < 0.0 || a <= 0.0) printf("Invalid arguments in routine gammq\n");
        if (x < (a+1.0)) {
                gser(&gamser,a,x,&gln);
                return 1.0-gamser;
        } else {
                gcf(&gammcf,a,x,&gln);
                return gammcf;
        }
}

void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
        double *b, double *siga, double *sigb, double *chi2, double *q)
{
        double gammq(double a, double x);
        int i;
        double wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;

        *b=0.0;
        if (mwt) {
                ss=0.0;
                for (i=1;i<=ndata;i++) {
                        wt=1.0/SQR(sig[i]);
                        ss += wt;
                        sx += x[i]*wt;
                        sy += y[i]*wt;
                }
        } else {
                for (i=1;i<=ndata;i++) {
                        sx += x[i];
                        sy += y[i];
                }
                ss=ndata;
        }
        sxoss=sx/ss;
        if (mwt) {
                for (i=1;i<=ndata;i++) {
                        t=(x[i]-sxoss)/sig[i];
                        st2 += t*t;
                        *b += t*y[i]/sig[i];
                }
        } else {
                for (i=1;i<=ndata;i++) {
                        t=x[i]-sxoss;
                        st2 += t*t;
                        *b += t*y[i];
                }
        }
        *b /= st2;
        *a=(sy-sx*(*b))/ss;
        *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
        *sigb=sqrt(1.0/st2);
        *chi2=0.0;
        if (mwt == 0) {
                for (i=1;i<=ndata;i++)
                        *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
                *q=1.0;
                sigdat=sqrt((*chi2)/(ndata-2));
                *siga *= sigdat;
                *sigb *= sigdat;
        } else {
                for (i=1;i<=ndata;i++)
                        *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
                *q=gammq(0.5*(ndata-2),0.5*(*chi2));
        }
}

/***! Count the number of lines in a file !***/

int ct_lines(char *file)
{
    FILE *fic;
    int n;
    char com[BUFSIZ];
    
    fic=fopen(file,"r");
    if (fic) {
	n=0;
	while(!feof(fic)) {
	    fscanf(fic,"%[^\n]\n",com);
	    n++;
	}
	fclose(fic);
	return n;
    } else {
        /*fprintf(stderr,"Can't read %s file !\n",file);*/
	return 0;
    }
}

/***! Count the number of columns in a file !***/

int ct_cols(char *file)
{
    FILE *fic;
    int n;
    char com[BUFSIZ];
    
    fic=fopen(file,"r");
    if (fic) {
	n=0;
	while(!feof(fic)) {
	    fscanf(fic,"%s",com);
	    n++;
	}
	fclose(fic);
	return n/ct_lines(file);
    } else {
        /*fprintf(stderr,"Can't read %s file !\n",file);*/
	return 0;
    }
}
