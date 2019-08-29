#include <math.h>
#include "fLegendre_header.h"

double fLegendre_func(double x){
  // Generating function for the Gegenbauer polynomials
  // parameters[0] = lambda
  // parameters[1] = alpha
  //  fprintf (stderr, "P[0] = %lf -- P[1] = %lf\n", parameters[0], parameters[1]);
  //  fprintf (stderr, "P[2] = %lf -- P[3] = %lf -- P[4] = %lf \n", parameters[2], parameters[3], parameters[4]);

  // return pow (1.0-2.0*x*parameters[1]+parameters[1]*parameters[1], -parameters[0]);
  // return pow(fabs(x), 3./2.);
  
  // k=0.1, sigma=2. -> ks = 0.2
  return  (1. + pow(parameters[2]*x*x, 2.) + 2.*parameters[2]*x*x)/(1. + 0.5*pow(parameters[3]*parameters[4], 2.)*x*x);
}

