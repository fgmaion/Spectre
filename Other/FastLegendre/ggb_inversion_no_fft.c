/*
        Routine: ggb_inversion_no_fft ()

        This routine computes the sum f(cos u) (or f(x)) of a Gegenbauer expansion from the Gegenbauer coefficients f_lambda
        by computing the Abel-type integral of a suitable (inverse) cosine-transform of the Gegenbauer coefficients, 
        according to the algorithm by E. De Micheli and G.A. Viano,

        Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
                   by Enrico De Micheli and Giovanni Alberto Viano.
                   Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008

        Input:
           f_lambda: array with the set of Gegenbauer coefficients;
           nmb_ggb_coeff: number of input Gegenbauer coefficients;
           lambda: the order lambda of the Gegenbauer polynomials;
           plottype: choose the output as f(x_n) or f(cos u_n);
	   worksize: gsl parameter, memory dimension to compute the Abel-type integral;
	   eps: gsl parameter, tolerance in computation of the Abel-type integral;
           nmb_out_samples: number of output samples.

        Output: 
	   f: array with the values of the function f(.) computed in nmb_out_samples abscissae.
           abscissae: array with the abscissae where the function has been computed.

        Use is made of a GSL routine for computing the Abel-type integral.

        COMPILE WITH: cc -O2 -c ggb_inversion_no_fft.c

        The calling program must be compiled with: cc -O2 calling.c ggb_inversion_no_fft.o -o calling -lgsl -lgslcblas -lm

 	Author: Enrico De Micheli

        Licensing: Free to use giving credit to the author
        
        This program is distributed WITHOUT ANY WARRANTY. 

*/

#include <stdio.h>
#include <gsl/gsl_integration.h>

struct my_f_params {double u; int N; double *cnl; double lambda;};    /* Parameters passed to the integrand: 
									u = lower limit of integration;
                                                                        N = number of Gegenbauer coeffiecients used in summation;
                                                                        cnl = pointer to the transformed Gegenbauer coefficients.  
								       */
static double	Abel();
static double 	integrand ();
static double   kernel ();
static double   Chi ();
int		ggb_inversion_no_fft ();
//------------------------------------------------------------------------------------------------------------------------

int
ggb_inversion_no_fft (f_lambda, nmb_ggb_coeff, nmb_ggb_coeff_to_use, lambda, plottype, worksize, eps, abscissae, f, nmb_out_samples)
double  *f_lambda,      // Input. Array with the Gegenbauer coefficients.
        *f,             // Output. Array with the samples f(cos u_n) or f(x_n). This array of length 'nmb_out_samples' must be provided by the caller.
        *abscissae,     // Output. Array with the abscissae u_n or x_n. This array of length 'nmb_out_samples' must be provided by the caller.
	eps,
        lambda;            // Input. Order lambda of the Gegenbauer polynomials. It must be l > 0.0.
int     nmb_ggb_coeff,  // Input. Number of input Gegenbauer coefficients. It must be a power of 2.
	nmb_ggb_coeff_to_use, // The actual number of Gegenbauer coefficients to use in the sum.
        plottype,       // plottype = 0 to have f(cos u_n) with u_n \in [-pi,pi]; plottype = 1 to have f (x_n) with x_n \in[-1,1]. Equispaced grid always.
	worksize,
        nmb_out_samples;// Input. Number of output samples of f(cos u) to be computed.
{
register int	k;
double          *c_lambda;          // Contains the transformed Gegenbauer coefficients c_n^lambda
double		s, delta, begin;

  if (plottype > 1){
      fprintf (stderr, "\nggb_inversion_no_fft: Plot type > 1: Output is set to f(x_n) with x_n in [-1,1]\n\n");
      plottype = 1;
  }
  if (plottype < 0){
      fprintf (stderr, "\nggb_inversion_no_fft: Plot type < 0: Output is set to f(cos u_n) with u_n in [-pi, pi]\n\n");
      plottype = 0;
  }
  if (lambda <= 0.0){
      fprintf (stderr, "\nggb_inversion_no_fft: The parameter Lambda must be positive\n\n");
      return (1);
  }
  if (nmb_ggb_coeff <= 0) {
      fprintf (stderr, "\nggb_inversion_no_fft: The number of Gegenbauer coefficients must be positive\n\n"); 
      return (1);
  }
  if (nmb_ggb_coeff_to_use < 1) {
      fprintf (stderr, "\nggb_inversion_no_fft: The number of Gegenbauer coefficients to use must be >= 1\n\n"); 
      return (1);
  }
  if (worksize < 1000) {
      fprintf (stderr, "\nggb_inversion_no_fft: The parameter WORKSIZE must be > 1000\n\n");
      return (1);
  }
  if (eps <= 0){
      fprintf (stderr, "\nggb_inversion_no_fft: Integration error must be positive\n\n");
      return (1);
  }
  if (nmb_out_samples < 2) {
      fprintf (stderr, "\nggb_inversion_no_fft: nmb_out_samples must be >= 2\n\n");
      return (1);
  }
  if (nmb_ggb_coeff_to_use > nmb_ggb_coeff) {
      nmb_ggb_coeff_to_use = nmb_ggb_coeff;
      fprintf (stderr, "\nggb_inversion_no_fft: nmb_ggb_coeff_to_use set to: %d\n", nmb_ggb_coeff);
  }

// Compute transformed Gegenbauer coefficients c_lambda
  c_lambda = (double *) malloc ((unsigned) nmb_ggb_coeff*sizeof(double));
  if (!c_lambda) {
      fprintf (stderr, "\nggb_inversion_no_fft: Allocation failure\n\n");
      return (1);
  }
  for (k = 0; k < nmb_ggb_coeff; k++)
       c_lambda[k] = exp (lgamma ((double) k+2.0*lambda)-lgamma ((double) k+1.0)) * f_lambda[k];

  nmb_out_samples += 2;       // Avoid computing the extrema

  if (plottype == 0) { // Compute f(x)
      delta = M_PI/(nmb_out_samples-1);
      begin = delta;
  }
  else { // Compute f(cos u)
      delta = 2.0/(nmb_out_samples-1);
      begin = -1.0 + delta;     // Avoid x = -1
  }
  for (s = begin, k = 0; k < nmb_out_samples-2; k++, s += delta) {
       abscissae[k] = s;
       f[k] = Abel ((plottype) ? acos (s) : s, c_lambda, lambda, worksize, eps, nmb_ggb_coeff_to_use);
  }

  return 0;
}
/* -------------------------------------------------------------------------------- */

// Compute the integral from u to Pi of Chi^lambda(t)*[2(cos(u)-cos(t)]^(lambda-1)

static double
Abel (u, cnl, lambda, worksize, eps, N)
double  u, 		// Running coordinate
	lambda,		// Order lambda of the Gegenbauer coefficients
	eps,
	*cnl;		// Array with the transformed Gegenbauer coefficients
int	N,		// Number of Gegenbauer coefficients used for the inversion
	worksize;
{
double                          result, error;
gsl_function                    F;                                     
struct my_f_params              params;         // Prameters to feed the integrand
gsl_integration_workspace      *ws;

  if (fabs(u-M_PI) < 1e-14)  // Is it Pi? Avoid computing the integral, it is zero
      return (0.0);

  params.u = u;
  params.N = N;
  params.cnl = cnl;
  params.lambda = lambda;
  F.params = &params;
  F.function = &integrand;

// Integral from u to Pi of Integrand 
  ws = gsl_integration_workspace_alloc (worksize);
  gsl_integration_qags (&F, u, M_PI, eps, eps, worksize, ws, &result, &error);
  gsl_integration_workspace_free (ws);

  return result / pow (sin (u), (double) (2.0*lambda-1.0));
}
/* -------------------------------------------------------------------------------- */

static double
integrand (x, par)
double  x;
void   *par;
{

  return kernel (x, ((struct my_f_params *) par)->u, ((struct my_f_params *) par)->lambda) * Chi(x, par);
}
/* ---------------------------------------------------------------------- */

static double
kernel (x, u, lambda)
double  x, u, lambda;
{

  return pow ((2.0*(cos(u)-cos(x))), (double) lambda-1.0);
}
/* ---------------------------------------------------------------------- */

static double
Chi (x, par)
double	x;
void   *par;
{
double		g, factor;
register int	k;
register double sum;

//  for (k = 0, val = 0.0; k < N; k++)
//       val += cnl[k] * cos((k+lambda)*x-lambda*M_PI);
//  run = ((struct my_f_params *) par)->cnl;

  for (k = 0, sum = 0.0; k < ((struct my_f_params *) par)->N; k++)
       sum += *(((struct my_f_params *) par)->cnl+k) * cos((k+((struct my_f_params *) par)->lambda)*x-((struct my_f_params *) par)->lambda*M_PI);

  g = tgamma (((struct my_f_params *) par)->lambda);	
  factor = pow ((double) 4.0, (double) 1.0-((struct my_f_params *) par)->lambda) / (g*g);

  return sum*factor;
}
/* ---------------------------------------------------------------------- */
