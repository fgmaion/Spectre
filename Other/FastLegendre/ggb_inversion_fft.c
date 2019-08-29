/*
        Routine: ggb_inversion_fft ()

        This routine computes the sum f(cos u) (or f(x)) of a Gegenbauer expansion from the Gegenbauer coefficients f_lambda
	by computing the Abel-type integral of a suitable (inverse) Fourier transformation of the Gegenbauer coefficients, 
	according to the algorithm by E. De Micheli and G.A. Viano.

        Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
		   by Enrico De Micheli and Giovanni Alberto Viano.
	           Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008

        Input:
           f_lambda: the set of Gegenbauer coefficients;
	   nmb_ggb_coeff: number of input Gegenbauer coefficients;
           lambda: the order lambda of the Gegenbauer polynomials;
	   quadrature: number of knots to use for the Newton-Cotes quadrature;
	   plottype: for choosing to have the output as f(x_n) or f(cos u_n);
	   zeropadding: factor for the extension of the length of the input Gegenbauer coefficients array to increase accuracy;
	   nmb_out_samples: number of output samples.

        Output: array with the values of the function f at nmb_out_samples abscissae.
                abscissae: array with the abscissae where the function has been computed.

        Use is made of a GSL routine for computing the FFT.

        COMPILE WITH: cc -O2 -c ggb_inversion_fft.c

        The calling program must be compiled with: cc -O2 calling.c ggb_inversion_fft.o -o calling -lgsl -lgslcblas -lm

	Author: Enrico De Micheli

        Licensing: Free to use giving credit to the author
        
        This program is distributed WITHOUT ANY WARRANTY. 

*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_fft_complex.h>

#define TRAPEZOIDAL	2
#define SIMPSON		3
#define SIMPSON38	4
#define BOOLE		5

#define REAL(z,i)       ((z)[2*(i)])
#define IMAG(z,i)       ((z)[2*(i)+1])

static double trapezoidal ();
static double simpson ();
static double simpson38 ();
static double boole ();
static double Abeltype ();
static void   Chi_function ();





// ------------------------------------------------------------------------------

int
ggb_inversion_fft (f_lambda, nmb_ggb_coeff, lambda, quadrature, plottype, zeropadding, abscissae, f, nmb_out_samples)
double  *f_lambda,      // Input. Array with the Gegenbauer coefficients.
	*f,		// Output. Array with the samples f(cos u_n) or f(x_n). This array of length 'nmb_out_samples' must be provided by the caller.
	*abscissae,	// Output. Array with the abscissae u_n or x_n. This array of length 'nmb_out_samples' must be provided by the caller.
        lambda;         // Input. Order lambda of the Gegenbauer polynomials. It must be l > 0.0.
int     nmb_ggb_coeff,  // Input. Number of input Gegenbauer coefficients. It must be a power of 2.
	quadrature, 	// Quadrature rule, number of points of the Newton-Cotes formula. [2-5].
	plottype,	// plottype = 0 to have f(cos u_n) with u_n \in [-pi,pi]; plottype = 1 to have f (x_n) with x_n \in[-1,1]. Equispaced grid always.
        zeropadding, 	// Input vector is extended to be (nmb_ggb_coeff * 2^zeropadding) long. The entries from nmb_ggb_coeff
		        // through (nmb_ggb_coeff * 2^zeropadding-1) are filled with zeros. This is to increase FFT accuracy under the
			// hypotheses that the Gegenbauer coefficients are nearly null for n sufficiently large, e.g., in the case of analytic functions.
        nmb_out_samples;// Input. Number of output samples of f(cos u) to be computed.
{
register int	k;
double		*c_lambda;						// Contains the transformed Gegenbauer coefficients c_n^lambda
double		s, delta, begin, factor, *chi;
double		(*quadrature_rule) ();

  if (zeropadding < 0) {
      fprintf (stderr, "\nggb_inversion_fft: Zeropadding factor must be positive\n\n");
      return (1);
  }
  zeropadding = 1 << zeropadding;	// 2^zeropadding
  if (nmb_out_samples < 2) {
      fprintf (stderr, "\nggb_inversion_fft: Number of output samples must be >= 2\n\n");
      return (1);
  }
  if (plottype > 1){
      fprintf (stderr, "\nggb_inversion_fft: Plot type > 1: Output is set to f(x_n) with x_n in [-1,1]\n\n");
      plottype = 1;
  }
  if (plottype < 0){
      fprintf (stderr, "\nggb_inversion_fft: Plot type < 0: Output is set to f(cos u_n) with u_n in [-pi, pi]\n\n");
      plottype = 0;
  }
  if (lambda <= 0.0){
      fprintf (stderr, "\nggb_inversion_fft: The parameter Lambda must be positive\n\n"); 
      return (1);
  }
// Sampling of the u-interval MUST NOT be the same as the sampling of the integrating t-interval: u_k != t_j for all suitable k and j
  if (nmb_out_samples == zeropadding*nmb_ggb_coeff+1){
      fprintf (stderr, "\nggb_inversion_fft: Choose a different number of sampling points: A sampling point coincides with the pole of the integrand\n\n");
      return (1);
  }
  if (nmb_ggb_coeff <= 0) {
      fprintf (stderr, "\nggb_inversion_fft: The number of Gegenbauer coefficients must be positive and power of 2\n\n"); // For FFT constraint
      return (1);
  }
  if (nmb_ggb_coeff & (nmb_ggb_coeff - 1) != 0) {
      fprintf (stderr, "\nggb_inversion_fft: The number of Gegenbauer coefficients is not a power of 2 \n\n"); // For FFT constraint
      return (1);
  }
  if (quadrature < 2 || quadrature > 5){
      fprintf (stderr, "\nggb_inversion_fft: Unknown quadrature rule. Number of integration knots in Newton-Cotes rule must be in [2, 5]\n\n");
      return (1);
  }
  switch (quadrature) {	// Switch according to number of points in the chosen QUADRATURE RULE
	  case TRAPEZOIDAL: quadrature_rule = trapezoidal;
		  break;
          case SIMPSON:     quadrature_rule = simpson;
		  break;
	  case SIMPSON38:   quadrature_rule = simpson38;
		  break;
	  case BOOLE:       quadrature_rule = boole;
		  break;
  }
// Compute the transformed Gegenbauer coefficients c_lambda
  c_lambda = (double *) malloc ((unsigned) nmb_ggb_coeff*sizeof(double));
  if (!c_lambda){
      fprintf (stderr, "\nggb_inversion_fft: Allocation failure in c_lambda\n\n");
      return (1);
  }
  for (k = 0; k < nmb_ggb_coeff; k++)
       c_lambda[k] = exp (lgamma ((double) k+2.0*lambda)-lgamma ((double) k+1.0)) * f_lambda[k];

// The array chi contains the values of the function Chi^lambda(t) at nmb_ggb_coeff points of the interval [-Pi,Pi]
// Update number of actual Gegenbauer coefficients: chi is zero padded, i.e., the input of the FFT is filled of zeros to increase 
// accuracy, assuming the GGB coefficients are almost null for large values of k. 
  chi = (double *) malloc ((unsigned) zeropadding*nmb_ggb_coeff*sizeof(double));
  if (!chi){
      fprintf (stderr, "\nggb_inversion_fft: Allocation failure in chi\n\n");
      return (1);
  }
  Chi_function (nmb_ggb_coeff, zeropadding, c_lambda, lambda, chi);

  nmb_out_samples += 2; 			// Avoid computaton at the endpoints, there is a pole.

  if (plottype == 0) { 			// f(cos u), then the array f contains f(cos u_n)
      delta = M_PI/(nmb_out_samples-1);	// This is delta_u
      begin = delta;			// Avoid u = 0, pi
  }
  else {				// f(x), then the array f contains f(x_n)
      delta = 2.0/(nmb_out_samples-1);	// This is delta_x
      begin = -1.0 + delta;		// Avoid x = +-1
  }

  factor = pow (2.0, (double) (lambda-1.0));
  for (s = begin, k = 0; k < nmb_out_samples-2; k++, s+= delta) {	// The abscissa s can be either u or x
       abscissae[k] = s;						// Save the abscissa
       f[k] = factor * Abeltype ((plottype) ? acos (s) : s, lambda, chi, zeropadding*nmb_ggb_coeff, quadrature_rule);
  }

  return 0;
}
// --------------------------------------------------------------------------------

// Compute the integral from u to Pi of Chi^lambda(t)*(cos(u)-cos(t))^(lambda-1)
// by using the Newton-Cotes equispaced quadrature.

static double
Abeltype (u, lambda, chi, N, quadrature)
double  u, 		// Running coordinate
	lambda,		// Order of the Gegenbauer coefficients
	*chi;		// Array containing the samples of the function Chi^lambda (t)
int	N;		// Number of Gegenbauer coefficients used for the inversion
double	(*quadrature)();
{
register int	k;
register double	*run;
double		*integrand, delta_t, t, t0, result, factor;
int		ki, kf, nu;

  if (fabs(u-M_PI) < 1e-14)  /* Is u = Pi? Avoid computing the integral, it is zero */
      return (0.0);

//  oldfactor = pow (2.0, (double) (lambda-1.0)) / pow (sin (u), (double) (2.0*lambda-1.0)); 
  factor = 1.0 / pow (sin (u), (double) (2.0*lambda-1.0));

// t-variable sampling period
  delta_t = 2.0 * M_PI / N;

// Set the two endpoints of the integration interval in terms of values of k. t_{N/2} = Pi !!!
  ki = (int) (u/delta_t) + 1;
  t0 = delta_t * ki;
  kf = N/2;	// t_{N/2} = Pi
  nu = kf-ki+1;
  run = chi + ki;

// Compute the value of the integrand at the samples t_j of the interval [0,2Pi]. Actually it'd suffice the interval [u, Pi]
  integrand = (double *) malloc ((unsigned) nu*sizeof(double));
  if (!integrand){
      fprintf (stderr, "\nAbeltype: Allocation failure in integrand\n\n");
      return (1);
  }
  for (k = 0, t = t0; k < nu; k++, t += delta_t) //integrand[k] = chi[k+ki] * pow (cos(u)-cos(t), (double) (lambda-1.0)); 
       integrand[k] = *run++ * pow (cos(u)-cos(t), (double) (lambda-1.0)); 

  result = (* quadrature) (integrand, nu) * delta_t; // Compute the integral

// Left endpoint correction.
  result += ((t0-u) * integrand[0]);

  return factor * result;
}
// --------------------------------------------------------------------------------

static void
Chi_function (nmb_ggb, zp, c_lambda, lambda, chi)
double *c_lambda, *chi, lambda;
int   	nmb_ggb,			// Number of Gegenbauer coefficients 
	zp;				// Zero padding factor
{
register int	k;
double		g, factor, t, delta;
int		N;
double		*data;

  g = tgamma (lambda);
  factor = 1.0 / (g*g);
  if (lambda != 1.0)
      factor *= pow ((double) 4.0, (double) 1.0-lambda);

  N = zp * nmb_ggb;
  data = (double *) calloc ((unsigned int) (2*N), sizeof (double)); 	// It must be initialized to zero

// Fill in the actually read Gegenbauer coefficients, the rest of the data are zeros.
  for (k = 0; k < nmb_ggb; k++)
       REAL (data, k) = c_lambda[k];

  gsl_fft_complex_radix2_backward (data, 1, N);		// FFT

  delta = 2.0 * M_PI / (double) N;
  for (k = 0, t = 0.0; k < N; k++, t += delta)
       chi[k] = factor * (cos (lambda*(t-M_PI)) * REAL (data, k) - sin (lambda*(t-M_PI)) * IMAG (data, k));
}
// ----------------------------------------------------------------------

static double
boole (integrand, nu)
double	*integrand;
int	nu;
{
int	k;
double	sum;

  for (k = 0, sum = 0.0; k <= nu-5; k += 4)
       sum += (7.0*integrand[k]+32.0*integrand[k+1]+12.0*integrand[k+2]+32.0*integrand[k+3]+7.0*integrand[k+4]);
  sum *= (2.0 / 45.0);

// Correction at the right end
  if (k == nu-4)
      sum += 3.0 * ((integrand[k]+3.0*integrand[k+1]+3.0*integrand[k+2]+integrand[k+3]) / 8.0);
  else
      if (k == nu-3) 
          sum += ((integrand[k]+4.0*integrand[k+1]+integrand[k+2]) / 3.0);
      else
          if (k == nu-2) 
              sum += ((integrand[k]+integrand[k+1]) / 2.0);
//  if k == nu-1 do nothing, it has already reached the right endpoint

  return sum;
}
// ----------------------------------------------------------------------

static double
simpson38 (integrand, nu)
double	*integrand;
int	nu;
{
int	k;
double	sum;

  for (k = 0, sum = 0.0; k <= nu-4; k += 3)
       sum += (integrand[k]+3.0*integrand[k+1]+3.0*integrand[k+2]+integrand[k+3]);
  sum *= 3.0 / 8.0;

// Correction at the right end
  if (k == nu-3) 
      sum += ((integrand[k]+4.0*integrand[k+1]+integrand[k+2]) / 3.0);
  else
     if (k == nu-2) 
         sum += ((integrand[k]+integrand[k+1]) / 2.0);
//  if k == nu-1 do nothing, it has already reached the right endpoint

  return sum;
}
// ---------------------------------------------------------------------- 

static double
simpson (integrand, nu)
double	*integrand;
int	nu;
{
int	k;
double	sum;

  for (k = 0, sum = 0.0; k <= nu-3; k += 2)
       sum += (integrand[k]+4.0*integrand[k+1]+integrand[k+2]);
  sum /= 3.0;

// Correction at the right end
  if (k == nu-2) 
      sum += ((integrand[k]+integrand[k+1]) / 2.0);
//  if k == nu-1 do nothing, it already reached the right endpoint

  return sum;
}
// ----------------------------------------------------------------------

static double
trapezoidal (integrand, nu)
double	*integrand;
int	nu;
{
int	k;
double	sum;

  for (k = 0, sum = 0.0; k <= nu-2; k += 1)
       sum += (integrand[k]+integrand[k+1]);
  sum /= 2.0;

  return sum;
}
// ----------------------------------------------------------------------
