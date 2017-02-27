/*
	Routine: ggb_coeff_fft ()

  	This routine computes the Gegenbauer coefficients f_lambda by using the algorithm by E. De Micheli and G.A. Viano.
	Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
        	   by Enrico De Micheli and Giovanni Alberto Viano.
                   Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008
	Input: 
           f(x): the function whose Gegenbauer coefficients are to be computed;
           N: the number of Gegenbauer coefficients to compute;
	   lambda: the order of the Gegenbauer polynomials;
	   worksize: gsl parameter for allocating memory in computing the Abel integral;
	   key: gsl parameter, key representing the number of points to use in the Gauss-Kronrod quadrature.
           eps: gsl parameter, tolerance admitted in the computation of the Abel integral.

        Output: Array with N Gegenbauer coefficients of the function f(x).

        Use is made of GSL routines for computing the Abel-type integral and the FFT.

        COMPILE WITH: cc -O2 -c ggb_coeff_fft.c

	The calling program must be compiled with: cc -O2 calling.c ggb_coeff_fft.o -o calling -lgsl -lgslcblas -lm

	Author: Enrico De Micheli

        Licensing: Free to be used giving credit to the author
        
        This program is distributed WITHOUT ANY WARRANTY. 

*/

#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fft_complex.h>

struct integrand_par {		// Parameters passed to the integrand:
	double (*f) (); 	// function to be expanded;
	double t; 		// for the upper limit of integration;
	double lambda; 		// Gegenbauer order;
	int n;			// degree of the Gegenbauer polynomial;
};

static void 	abel_transform ();
static double 	integral_abel ();
static double 	integrand ();

int	ggb_coeff_fft ();

//--------------------------------------------------------------------------------------------------------

int
ggb_coeff_fft (f, ggb, N, lambda, worksize, key, eps)
double	(*f) (), 	// Input. Function whose Gegenbauer coefficients are to be computed.
	*ggb, 		// Output. Array provided by the calling routine. It will contain the Gegenbauer coefficients.
	lambda, 	// Input. Order lambda of the Gegenbauer polynomials.It must be lambda > -0.5.
	eps;		// Input. Tolerance for the computation of the Abel-type integral. Typically eps = 1e-11.
int	N, 		// Input. Number of Gegenbauer coefficients to be computed. It must be a power of 2.
	worksize, 	// Input. Worksize for computing the Abel integral. Typically worksize = 5000.
	key;		// Input. Number of Gauss-Kronrod quadrature points. Typically key = 2. It must be in [1, 6].
{
register unsigned int	k;
double			*data = NULL, 
			delta,        	// Sampling period of the Abel transform which feeds the FFT.
			factor,
                	t; 

  if (ggb == NULL) {
      fprintf (stderr, "\nggb_coeff_fft: Missing output Gegenbauer coefficients memory array\n\n");
      return (1);
  }
  if (N <= 0) {
      fprintf (stderr, "\nggb_coeff_fft: The number of Gegenbauer coefficients must be positive and power of 2\n\n");
      return (1);
  }
  if (N & (N - 1) != 0) {
      fprintf (stderr, "\nggb_coeff_fft: The number of Gegenbauer coefficients is not a power of 2 \n\n");
      return (1);
  }
  N *= 2;

  if (lambda <= -0.5) {
      fprintf (stderr, "\nggb_coeff_fft: The parameter LAMBDA must be larger than -0.5\n\n");
      return (1);
  }

  if (key < 1 || key > 6) {
      fprintf (stderr, "\nggb_coeff_fft: The parameter KEY must be 1 <= KEY <= 6\n\n");
      return (1);
  }

  if (worksize < 1000) {
      fprintf (stderr, "\nggb_coeff_fft: The parameter WORKSIZE must be > 1000\n\n");
      return (1);
  }

  if (eps <= 0) {
      fprintf (stderr, "\nggb_coeff_fft: Integration relative error must be positive\n\n");
      return (1);
  }

  data = (double *) malloc ((unsigned) 2*N*sizeof(double));
  if (!data) {
      fprintf (stderr, "\nggb_coeff_fft: Allocation failure for data array\n\n");
      return (1);
  }

// Compute the Abel Transform at sampled points t_k on the regular grid of the interval [-Pi, Pi),
// and put the (complex) result in the array data, according to the following rule (which comes from the gsl library): 
// data[2k] = Real part -- data[2k+1] = Imaginary part.
  for (k = 0, t = -M_PI, delta = 2.0 * M_PI / N; k < N; k++, t += delta)
       abel_transform (f, t, lambda, worksize, key, eps, data+2*k);

// Compute the FFT of the Abel transform contained in the array data, where the final result in then stored.
  gsl_fft_complex_radix2_backward (data, 1, N);

// Compute the Gegenbauer coefficients from the array data.
// Remind: (k & 1 == 1) ? k is ODD : k is EVEN). Useful to generate (-1)^k.
  factor = pow (2.0, lambda) / (N*lambda); // It is the result of (pow (2.0, lambda-1.0) / (M_PI*lambda)) * delta
  for (k = 0; k < N/2; k++)
//     *ggb++ = oldfactor * (k+lambda) * delta * *(data+2*k) * ((k & 1 == 1) ? /* k is odd */ -1 : /* k is even */ 1); 
       *ggb++ = factor * (k+lambda) * *(data+2*k) * ((k & 1 == 1) ? /* k is odd */ -1 : /* k is even */ 1); 

  return (0);
}
/* -------------------------------------------------------------------------------- */

/*
    Keys for the Gauss-Kronrad quadrature rules in the routine gsl_integration_qag:

          GSL_INTEG_GAUSS15  (key = 1)
          GSL_INTEG_GAUSS21  (key = 2)
          GSL_INTEG_GAUSS31  (key = 3)
          GSL_INTEG_GAUSS41  (key = 4)
          GSL_INTEG_GAUSS51  (key = 5)
          GSL_INTEG_GAUSS61  (key = 6)
*/

static double
integral_abel (f, t, lambda, worksize, key, eps)
double 	(*f) (), t, lambda, eps;
int	worksize, key;
{
double				result, abs_error;
gsl_function                    F;                      		/* Struct funzione da integrare */
struct integrand_par            par;                         	/* Parametri da passare alla funzione da integrare */
gsl_integration_workspace      *ws;

  par.f = f;
  par.t = t;
  par.lambda = lambda;
  F.params = &par;
  F.function = &integrand;

  ws = gsl_integration_workspace_alloc (worksize);

  gsl_integration_qag (&F, 0.0, pow (1.0-cos(t), lambda), eps, eps, worksize, key, ws, &result, &abs_error);

  gsl_integration_workspace_free (ws);

  return (result);
}
/* -------------------------------------------------------------------------------- */

static double
integrand (x, par)
double  x;
void   *par;
{
//return f (pow (x, 1.0/lambda)+cos(t), par);
  return ((struct integrand_par *) par)->f (pow (x, 1.0/((struct integrand_par *) par)->lambda)+cos(((struct integrand_par *) par)->t), par);
}
/* -------------------------------------------------------------------------------- */

static void
abel_transform (f, t, lambda, worksize, key, eps, run)
double  (*f) (), t, lambda, eps, *run;
int	worksize, key;
{
double	tmp, re, im;

  if (fabs(t) < 1e-14){  /* Is t == 0? Avoid computing the integral, it's null */
      *run = 0.0;
      *(run+1) = 0.0;
      return;
  }
  else
//**  tmp = integral_abel (f, t, lambda, worksize, key, eps) * pow (2.0, lambda-1.0) / (M_PI*lambda); Numerical factor moved.
      tmp = integral_abel (f, t, lambda, worksize, key, eps);

  if (t >=0) {
       re = cos (lambda*(t-M_PI));
       im = sin (lambda*(t-M_PI));
  }
  else {
       re = cos (lambda*(t+M_PI));
       im = sin (lambda*(t+M_PI));
  }

  *run = re*tmp;
  *(run+1) = im*tmp;
}
/* -------------------------------------------------------------------------------- */
