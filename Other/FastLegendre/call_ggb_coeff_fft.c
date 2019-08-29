#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "fLegendre_header.h"


int fastLegendre_init(){  
  parameters = malloc(5*sizeof(*parameters));
  ggb        = malloc ((unsigned) 256*sizeof(double));
  
  // Default values of parameters to feed the input function
  parameters[0] = 0.5;	    // parameters[0] is set to lambda.
  parameters[1] = 0.5;		// Default value of parameters[1]

  return 0;
}


int fastLegendre(){  
  FILE		                *fpout;
  char		     filename_out[200];
  
  fastLegendre_init();
  
  sprintf(filename_out, "/disk1/mjw/HOD_MockRun/Scripts/FastLegendre/multipole_moments.dat");

  fpout = fopen(filename_out, "w");
  
  int ii, k;
  
  double beta  = 0.5;
  double sigma = 2.0;
  double  kmod;

  parameters[2] = beta;
  parameters[3] = sigma;
  
  for(ii=0; ii<200; ii++){
    kmod = pow(10., -2. + ii/50.);
    
    parameters[4] = kmod;

    if(ggb_coeff_fft(fLegendre_func, ggb, 256, 0.5, 5000, 2, 1e-12) == 1)  error (" Error computing the Gegenbauer coefficients");

    fprintf (fpout, "%.6e \t", kmod);

    for(k=0; k<256; k+=2)  fprintf (fpout, "%.6e \t", ggb[k]);

    fprintf (fpout, "\n");
  }
  
  fclose (fpout); 

  printf("\n\nFast Legendre called.");

  return 0;
}
