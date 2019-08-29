/*
        Routine ggb_coeff_jacobi () 

        This routine computes N Gegenbauer coefficients of order l for the function f() by using the Gauss-Jacobi quadrature 
	for computing the Abel integral and one FFT, according to the algorithm by De Micheli & Viano.

        Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
                   by Enrico De Micheli and Giovanni Alberto Viano.
                   Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008

        Input: 
           f(x): the function whose Gegenbauer coefficients are to be computed;
           N: the number of Gegenbauer coefficients to compute;
           lambda: the order of the Gegenbauer polynomials;
	   order: the order of the Gauss-Jacobi quadrature.

        Output: Array with N Gegenbauer coefficients of the function f(x).

	The Fast Fourier Transform routine comes from the gsl library.

	Part of this code, precisely, the computation of abscissae and weights for the Gauss-Jacobi quadrature,
        is written by other authors as acknowledged in each single subroutine of this file.

        COMPILE with: g++ -O2 -c ggb_coeff_jacobi.cpp

	The calling main program must be compiled with: g++ -O2 calling.cpp -o calling ggb_coeff_jacobi.o -lgsl -lgslcblas -lm

	Author: Enrico De Micheli

        Licensing: This code is distributed under the GNU GPL license 
                   published by the Free Software Foundation, Inc., 51 Franklin Street, Boston, MA 02110-1301, USA.
        
        This program is distributed WITHOUT ANY WARRANTY. 

	Copyright 2012, 2013 Enrico De Micheli, Giovanni Alberto Viano
*/
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_fft_complex.h>
#include <iostream>

#define REAL(z,i)       ((z)[2*(i)])
#define IMAG(z,i)       ((z)[2*(i)+1])

using namespace std;

int 		ggb_coeff_jacobi (double (*f)(double x), double *ggb, int N, double lambda, int order);
static void 	Jacobi_quadrature (int order, double alpha, double beta, double a, double b, double abscissae[], double weights[]);
static double 	Jacobi_integral (double (*f)(double x), double t, double l, int n);
static void 	abel_transform (double (*f)(double x), double t, double l, int order, double hf[]);
//
// The following routines have been made by various authors for computing abscissas and weights for the Jacobi quadrature. 
// See the author at the top of each routine.
//
static void 	cdgqf (int nt, int kind, double alpha, double beta, double t[], double wts[]);
static void 	cgqf (int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[]);
static double 	class_matrix (int kind, int m, double alpha, double beta, double aj[], double bj[]);
static void 	imtqlx (int n, double d[], double e[], double z[]);
static void 	parchk (int kind, int m, double alpha, double beta);
static double 	r8_abs (double x);
static double 	r8_epsilon ();
static double 	r8_sign (double x);
static void 	scqf (int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[], double st[], 
                      int kind, double alpha, double beta, double a, double b);
static void 	sgqf (int nt, double aj[], double bj[], double zemu, double t[], double wts[]);
//-----------------------------------------------------------------------------------------

int
ggb_coeff_jacobi (double (*f)(double x), double *ggb, int N, double lambda, int order)
{
register int    k;
int		N2;
double		delta,		// Time sampling rate 
		factor,
		t; 
double		*data, hf[2];

  if (ggb == NULL) {
      fprintf (stderr, "\nggb_coeff_jacobi: Missing output Gegenbauer coefficients memory array\n\n");
      return (1);
  }
  if (N <= 0) {
      fprintf (stderr, "\nggb_coeff_jacobi: The number of Gegenbauer coefficients must be positive and power of 2\n\n");
      return (1);
  }
  if (N & (N - 1) != 0) {
      fprintf (stderr, "\nggb_coeff_jacobi: The number of Gegenbauer coefficients is not a power of 2 \n\n");
      return (1);
  }

  if (lambda <= -0.5) {
      fprintf (stderr, "\nggb_coeff_jacobi: The parameter LAMBDA must be larger than -0.5\n\n");
      return (1);
  }

  if (order <= 0) {
      fprintf (stderr, "\nggb_coeff_jacobi: The order of the Jacobi quadrature must be positive\n\n");
      return (1);
  }

  N2 = 2*N;			// This is doubled because fft works on doubled data, they are symmetric
  data = new double[2*N2];	// This is again doubled because data is complex

  for (k = 0, t = -M_PI, delta = M_PI / N; k < N2; k++, t += delta) {
       abel_transform (f, t, lambda, order, hf);	 	// The result is returned in the variable hf[]
       REAL (data, k) = hf[0];
       IMAG (data, k) = hf[1];
  }

  gsl_fft_complex_radix2_backward (data, 1, N2);

// oldfactor = pow (2.0, lambda-1.0) / M_PI; newfactor = pow (2.0, lambda-1.0) / M_PI * delta = pow (2.0, lambda-1.0) / N;
  factor = pow (2.0, lambda-1.0) / N;
  for (k = 0; k < N; k++)
       *(ggb+k) = ((k % 2) ? -factor : factor) * REAL (data, k) * (k+lambda);  		// Gegenbauer coefficients

  return (0);
}
/* ---------------------------------------------------------------------- */

static double
Jacobi_integral (double (*f)(double x), double t, double lambda, int n)
{
double		sum, s2; 
double		*xi;
double		*wi;
register int	j;

  xi = new double[n];
  wi = new double[n];

  Jacobi_quadrature (n, 0.0, lambda-1.0, -1.0, 1.0, xi, wi);

  s2  = sin (t/2.0);
  s2 *= s2;		// Squared sin(t/2)

  for (j = 0, sum = 0.0; j < n; j++)
       sum += (wi[j] * (*f) (s2*(xi[j]-1)+1)); 

  delete [] wi;
  delete [] xi;

  return pow (s2, lambda) * sum;
}
/*******************************************************************************/

static void
abel_transform (double (*f)(double x), double t, double lambda, int order, double hf[])
{
double          val, re, im;

// Compute exp[i*lambda*(t-Epsilon(t)*Pi)]
  if (t >=0) {
       re = cos (lambda*(t-M_PI));
       im = sin (lambda*(t-M_PI));
  }
  else {
       re = cos (lambda*(t+M_PI));
       im = sin (lambda*(t+M_PI));
  }

//val = pow (2.0, lambda-1.0) * Jacobi_integral (f, t, lambda, order) / M_PI;
  val = Jacobi_integral (f, t, lambda, order);

  hf[0] = re*val;
  hf[1] = im*val;
}
/*******************************************************************************/

static void 
Jacobi_quadrature (int order, double alpha, double beta, double a, double b, double abscissae[], double weights[])
{
int 	kind, j;
double 	*w, *x;

  w = new double[order];
  x = new double[order];
  
  kind = 4;
  cgqf (order, kind, alpha, beta, a, b, x, w);

  for (j = 0; j < order; j++) {
       weights[j] = w[j];
       abscissae[j] = x[j];
  }

  delete [] w;
  delete [] x;
}
//****************************************************************************

static void 
cdgqf (int nt, int kind, double alpha, double beta, double t[], double wts[])

//****************************************************************************
//
//  Purpose:
//
//    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with a classical weight function with default values for A and B,
//    and only simple knots.
//
//    There are no moments checks and no printing is done.
//
//    Use routine EIQFS to evaluate a quadrature computed by CGQFS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
double *aj;
double *bj;
double zemu;

  parchk (kind, 2 * nt, alpha, beta);
//
//  Get the Jacobi matrix and zero-th moment.
//
  aj = new double[nt];
  bj = new double[nt];

  zemu = class_matrix (kind, nt, alpha, beta, aj, bj);
//
//  Compute the knots and weights.
//
  sgqf (nt, aj, bj, zemu, t, wts);

  delete [] aj;
  delete [] bj;

  return;
}
//****************************************************************************

static void 
cgqf (int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[])

//****************************************************************************
//
//  Purpose:
//
//    CGQF computes knots and weights of a Gauss quadrature formula.
//
//  Discussion:
//
//    The user may specify the interval (A,B).
//
//    Only simple knots are produced.
//
//    Use routine EIQFS to evaluate this quadrature formula.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints, or
//    other parameters.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
int 	i;
int 	*mlt;
int 	*ndx;
//
//  Compute the Gauss quadrature formula for default values of A and B.
//
  cdgqf (nt, kind, alpha, beta, t, wts);
//
//  Prepare to scale the quadrature formula to other weight function with 
//  valid A and B.
//
  mlt = new int[nt];
  ndx = new int[nt];

  for (i = 0; i < nt; i++){
       mlt[i] = 1;
       ndx[i] = i + 1;
  }

  scqf (nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

  delete [] mlt;
  delete [] ndx;

  return;
}
//****************************************************************************

static double 
class_matrix (int kind, int m, double alpha, double beta, double aj[], double bj[])

//****************************************************************************
//
//  Purpose:
//
//    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
//
//  Discussion:
//
//    This routine computes the diagonal AJ and sub-diagonal BJ
//    elements of the order M tridiagonal symmetric Jacobi matrix
//    associated with the polynomials orthogonal with respect to
//    the weight function specified by KIND.
//
//    For weight functions 1-7, M elements are defined in BJ even
//    though only M-1 are needed.  For weight function 8, BJ(M) is
//    set to zero.
//
//    The zero-th moment of the weight function is returned in ZEMU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the Jacobi matrix.
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
//    of the Jacobi matrix.
//
//    Output, double CLASS_MATRIX, the zero-th moment.
//
{
double 	a2b2;
double 	ab;
double 	aba;
double 	abi;
double 	abj;
double 	abti;
double 	apone;
int 	i;
double 	pi = 3.14159265358979323846264338327950;
double 	temp;
double 	temp2;
double 	zemu;

  temp = r8_epsilon ();

  parchk (kind, 2 * m - 1, alpha, beta);

  temp2 = 0.5;

  if (500.0 * temp < r8_abs (pow (tgamma (temp2), 2) - pi)){
      cout << "\n";
      cout << "CLASS_MATRIX - Fatal error!\n";
      cout << "Gamma function does not match machine parameters.\n";
      exit (1);
  }

  if ( kind == 1 )
  {
    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 2 )
  {
    zemu = pi;

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] =  sqrt ( 0.5 );
    for ( i = 1; i < m; i++ )
    {
      bj[i] = 0.5;
    }
  }
  else if ( kind == 3 )
  {
    ab = alpha * 2.0;
    zemu = pow ( 2.0, ab + 1.0 ) * pow ( tgamma ( alpha + 1.0 ), 2 )
      / tgamma ( ab + 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
    for ( i = 2; i <= m; i++ )
    {
      bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * pow ( i + alpha, 2 ) - 1.0 ) );
    }
  }
  else if ( kind == 4 )
  {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow ( 2.0, ab + 1.0 ) * tgamma ( alpha + 1.0 ) 
      * tgamma ( beta + 1.0 ) / tgamma ( abi );
    aj[0] = ( beta - alpha ) / abi;
    bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
      / ( ( abi + 1.0 ) * abi * abi ) );
    a2b2 = beta * beta - alpha * alpha;

    for ( i = 2; i <= m; i++ )
    {
      abi = 2.0 * i + ab;
      aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi * abi;
      bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
        / ( ( abi - 1.0 ) * abi ) );
    }
  }
  else if ( kind == 5 )
  {
    zemu = tgamma ( alpha + 1.0 );

    for ( i = 1; i <= m; i++ )
    {
      aj[i-1] = 2.0 * i - 1.0 + alpha;
      bj[i-1] = sqrt ( i * ( i + alpha ) );
    }
  }
  else if ( kind == 6 )
  {
    zemu = tgamma ( ( alpha + 1.0 ) / 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
    }
  }
  else if ( kind == 7 )
  {
    ab = alpha;
    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 8 )
  {
    ab = alpha + beta;
    zemu = tgamma ( alpha + 1.0 ) * tgamma ( - ( ab + 1.0 ) ) 
      / tgamma ( - beta );
    apone = alpha + 1.0;
    aba = ab * apone;
    aj[0] = - apone / ( ab + 2.0 );
    bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for ( i = 2; i <= m; i++ )
    {
      abti = ab + 2.0 * i;
      aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
    }

    for ( i = 2; i <= m - 1; i++ )
    {
      abti = ab + 2.0 * i;
      bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
        / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
    }
    bj[m-1] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      bj[i] =  sqrt ( bj[i] );
    }
  }

  return zemu;
}
//****************************************************************************

static void 
imtqlx (int n, double d[], double e[], double z[])

//****************************************************************************
//
//  Purpose:
//
//    IMTQLX diagonalizes a symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This routine is a slightly modified version of the EISPACK routine to 
//    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 
//
//    The authors thank the authors of EISPACK for permission to use this
//    routine. 
//
//    It has been modified to produce the product Q' * Z, where Z is an input 
//    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
//    The changes consist (essentialy) of applying the orthogonal transformations
//    directly to Z as they are generated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//    Roger Martin, James Wilkinson,
//    The Implicit QL Algorithm,
//    Numerische Mathematik,
//    Volume 12, Number 5, December 1968, pages 377-383.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D(N), the diagonal entries of the matrix.
//    On output, the information in D has been overwritten.
//
//    Input/output, double E(N), the subdiagonal entries of the 
//    matrix, in entries E(1) through E(N-1).  On output, the information in
//    E has been overwritten.
//
//    Input/output, double Z(N).  On input, a vector.  On output,
//    the value of Q' * Z, where Q is the matrix that diagonalizes the
//    input symmetric tridiagonal matrix.
//
{
double 	b;
double 	c;
double 	f;
double 	g;
int 	i;
int 	ii;
int 	itn = 30;
int 	j;
int 	k;
int     l;
int     m;
int 	mml;
double 	p;
double 	prec;
double 	r;
double 	s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( r8_abs ( e[m-1] ) <= prec * ( r8_abs ( d[m-1] ) + r8_abs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        cout << "\n";
        cout << "IMTQLX - Fatal error!\n";
        cout << "  Iteration limit exceeded\n";
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + r8_abs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( r8_abs ( g ) <= r8_abs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
//
//  Sorting.
//
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
//****************************************************************************

static void 
parchk (int kind, int m, double alpha, double beta)

//****************************************************************************
//
//  Purpose:
//
//    PARCHK checks parameters ALPHA and BETA for classical weight functions. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
//
//    Input, int M, the order of the highest moment to
//    be calculated.  This value is only needed when KIND = 8.
//
//    Input, double ALPHA, BETA, the parameters, if required
//    by the value of KIND.
//
{
double tmp;

  if ( kind <= 0 )
  {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  KIND <= 0.\n";
    exit ( 1 );
  }
//
//  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
//
  if ( 3 <= kind && alpha <= -1.0 )
  {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  3 <= KIND and ALPHA <= -1.\n";
    exit ( 1 );
  }
//
//  Check BETA for Jacobi.
//
  if ( kind == 4 && beta <= -1.0 )
  {
    cout << "\n";
    cout << "PARCHK - Fatal error!\n";
    cout << "  KIND == 4 and BETA <= -1.0.\n";
    exit ( 1 );
  }
//
//  Check ALPHA and BETA for rational.
//
  if ( kind == 8 )
  {
    tmp = alpha + beta + m + 1.0;
    if ( 0.0 <= tmp || tmp <= beta )
    {
      cout << "\n";
      cout << "PARCHK - Fatal error!\n";
      cout << "  KIND == 8 but condition on ALPHA and BETA fails.\n";
      exit ( 1 );
    }
  }
  return;
}
//****************************************************************************

static double 
r8_abs (double x)

//****************************************************************************
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
double value;

  if (0.0 <= x)
      value = x;
  else
      value = -x;

  return value;
}
//****************************************************************************

static double 
r8_epsilon ()

//****************************************************************************
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the 
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
double value;

  value = 1.0;

  while (1.0 < (double) (1.0 + value))
         value = value / 2.0;

  value = 2.0 * value;

  return value;
}
//****************************************************************************

static double 
r8_sign (double x)

//****************************************************************************
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
double value;

  if (x < 0.0)
      value = -1.0;
  else
      value = 1.0;

  return value;
}
//****************************************************************************

static void 
scqf (int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
      double swts[], double st[], int kind, double alpha, double beta, double a, double b) 
//****************************************************************************
//
//  Purpose:
//
//    SCQF scales a quadrature formula to a nonstandard interval.
//
//  Discussion:
//
//    The arrays WTS and SWTS may coincide.
//
//    The arrays T and ST may coincide.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double T[NT], the original knots.
//
//    Input, int MLT[NT], the multiplicity of the knots.
//
//    Input, double WTS[NWTS], the weights.
//
//    Input, int NWTS, the number of weights.
//
//    Input, int NDX[NT], used to index the array WTS.  
//    For more details see the comments in CAWIQ.
//
//    Output, double SWTS[NWTS], the scaled weights.
//
//    Output, double ST[NT], the scaled knots.
//
//    Input, int KIND, the rule.
//    1, Legendre,             (a,b)       1.0
//    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
//    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
//    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
//    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
//    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
//    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
//    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
//    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
//
//    Input, double ALPHA, the value of Alpha, if needed.
//
//    Input, double BETA, the value of Beta, if needed.
//
//    Input, double A, B, the interval endpoints.
//
{
double al;
double be;
int i;
int k;
int l;
double p;
double shft;
double slp;
double temp;
double tmp;

  temp = r8_epsilon ( );

  parchk ( kind, 1, alpha, beta );

  if ( kind == 1 )
  {
    al = 0.0;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 2 )
  {
    al = -0.5;
    be = -0.5;
    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 3 )
  {
    al = alpha;
    be = alpha;
    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 4 )
  {
    al = alpha;
    be = beta;

    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 5 )
  {
    if ( b <= 0.0 )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  B <= 0\n";
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 6 )
  {
    if ( b <= 0.0 )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  B <= 0.\n";
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / sqrt ( b );
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 7 )
  {
    al = alpha;
    be = 0.0;
    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 8 )
  {
    if ( a + b <= 0.0 )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  A + B <= 0.\n";
      exit ( 1 );
    }
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  }
  else if ( kind == 9 )
  {
    al = 0.5;
    be = 0.5;
    if ( r8_abs ( b - a ) <= temp )
    {
      cout << "\n";
      cout << "SCQF - Fatal error!\n";
      cout << "  |B - A| too small.\n";
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }

  p = pow ( slp, al + be + 1.0 );

  for ( k = 0; k < nt; k++ )
  {
    st[k] = shft + slp * t[k];
    l = abs ( ndx[k] );

    if ( l != 0 )
    {
      tmp = p;
      for ( i = l - 1; i <= l - 1 + mlt[k] - 1; i++ )
      {
        swts[i] = wts[i] * tmp;
        tmp = tmp * slp;
      }
    }
  }
  return;
}
//****************************************************************************

static void 
sgqf (int nt, double aj[], double bj[], double zemu, double t[], double wts[])

//****************************************************************************
//
//  Purpose:
//
//    SGQF computes knots and weights of a Gauss Quadrature formula.
//
//  Discussion:
//
//    This routine computes all the knots and weights of a Gauss quadrature
//    formula with simple knots from the Jacobi matrix and the zero-th
//    moment of the weight function, using the Golub-Welsch technique.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 January 2010
//
//  Author:
//
//    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Sylvan Elhay, Jaroslav Kautsky,
//    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
//    Interpolatory Quadrature,
//    ACM Transactions on Mathematical Software,
//    Volume 13, Number 4, December 1987, pages 399-415.
//
//  Parameters:
//
//    Input, int NT, the number of knots.
//
//    Input, double AJ[NT], the diagonal of the Jacobi matrix.
//
//    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
//    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.
//
//    Input, double ZEMU, the zero-th moment of the weight function.
//
//    Output, double T[NT], the knots.
//
//    Output, double WTS[NT], the weights.
//
{
int i;
//
//  Exit if the zero-th moment is not positive.
//
  if ( zemu <= 0.0 )
  {
    cout << "\n";
    cout << "SGQF - Fatal error!\n";
    cout << "  ZEMU <= 0.\n";
    exit ( 1 );
  }
//
//  Set up vectors for IMTQLX.
//
  for ( i = 0; i < nt; i++ )
  {
    t[i] = aj[i];
  }
  wts[0] = sqrt ( zemu );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
//
//  Diagonalize the Jacobi matrix.
//
  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i];
  }

  return;
}
//****************************************************************************
