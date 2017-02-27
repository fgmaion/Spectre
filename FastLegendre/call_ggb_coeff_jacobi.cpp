/*
 
    	This program computes N Gegenbauer coefficients by using Gauss-Jacobi quadrature and FFT. 

    	This program calls the routine ggb_coeff_jacobi () placed in the file ggb_coeff_jacobi.cpp
 
    	COMPILE WITH: g++ call_ggb_coeff_jacobi.cpp -O2 ggb_coeff_jacobi.o function.o -o call_ggb_coeff_jacobi -lgsl -lgslcblas -lm

        Licensing: This code is distributed under the GNU GPL license 
                   published by the Free Software Foundation, Inc., 51 Franklin Street, Boston, MA 02110-1301, USA.
        
        This program is distributed WITHOUT ANY WARRANTY. 

        Copyright 2012, 2013 Enrico De Micheli, Giovanni Alberto Viano
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

char	*progname;
double	parameters[5];  // Five parameters if needed, for the input function could be sufficient

int 	    main (int argc, char **argv);
int	    ggb_coeff_jacobi (double (*f)(double x), double *ggb, int N, double lambda, int order);
static void error (char *message);

extern "C" double f (double x);
// ------------------------------------------------------------------------------

main (int argc, char **argv)
{
register int    k;
int 		order = 20, 	// Order of the Gauss-Jacobi quadrature 
		N = 128;   	// Number of Gegenbauer coefficients to be computed
double		*ggb = NULL, lambda = 0.5;
FILE		*fpout;
char		filename_out [64], err[64], common[64];
int		time = 1, verbose = 1, option = 0, write_output = 0, flag_lambda = 0;
struct timeval  start, end;
long            mtime, seconds, useconds;

  // Default values of parameters to feed the input function
  parameters[0] = lambda;       // parameters[0] is set to lambda.
  parameters[1] = 0.5;          // Default value of parameters[1]

  progname = argv[0];

  opterr = 0;
  while ((option = getopt(argc, argv,":vtl:n:1:2:3:4:j:")) != -1) {
        switch (option) {
             case 'v' : verbose = 0;
                        break;
             case 't' : time = 0;
                        break;
             case 'l' : lambda = atof(optarg);
                        if (lambda <= -0.5)
                            error ((char *) "GTJACOBI: The parameter LAMBDA must be larger than -0.5");
  			parameters[0] = lambda;
                        break;
             case '1' : parameters[1] = atof(optarg);
                        break;
             case '2' : parameters[2] = atof(optarg);
                        break;
             case '3' : parameters[3] = atof(optarg);
                        break;
             case '4' : parameters[4] = atof(optarg);
                        break;
             case 'n' : N = atoi(optarg);
                        if (N <= 0)
                            error ((char *) "GTJACOBI: The number of Gegenbauer coefficients must be positive");
                        break;
             case 'j' : order = atoi(optarg);
                        if (order <= 0)
                            error ((char *) "GTJACOBI: The order of the Jacobi quadrature must be positive");
                        break;
             case ':' : fprintf (stderr, "\n %s -GTJACOBI Error: Option: `-%c' needs a value \n\n", progname, optopt);
                        exit (1);
             case '?' : sprintf (err, "GTJACOBI: No such option: `-%c' ", optopt);
                        error (err);
        }
  }

  if (argc-optind > 1)
      error ((char *) "GTJACOBI: Cannot parse correctly the COMMON output filename");
  else
      if (argc-optind == 1) {
          strcpy (common, argv[optind]);
          sprintf (filename_out, "%s_ggb_coefficients_jacobi.dat", common);
          if ((fpout = fopen (filename_out, "w")) == NULL)
              error ((char *) "GTJACOBI: Can't open output file for the Gegenbauer coefficients");
          write_output = 1;
      }

  if (write_output) {
      fprintf (fpout, "# Gegenbauer coefficients for the function f(x) computed with Gauss-Jacobi quadrature by the De Micheli & Viano algorithm \n");
      fprintf (fpout, "# Lambda = %lf \n", lambda);
      fprintf (fpout, "#@ %lf \n", lambda);
      fprintf (fpout, "# n \t f_n^Lambda \n");
  }

  if (verbose == 0 && write_output == 0)
      verbose = 1;

  if (time)
      gettimeofday(&start, NULL);

// THE CORE OF THE COMPUTATION
  ggb = new double[N];
  if (ggb_coeff_jacobi (f, ggb, N, lambda, order))
      exit (1);

  if (time) {
      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
  }

  if (verbose) {
      fprintf (stderr, "GTJACOBI: Gegenbauer coefficients computed by Gauss-Jacobi quadrature and FFT \n");
      fprintf (stderr, "GTJACOBI: Lambda = %lf -- Jacobi quadrature order = %d \n", lambda, order);
      fprintf (stderr, "# n \t f_n^Lambda \n");
  }

  if (verbose)
      for (k = 0; k < N; k++)
           fprintf (stderr, "%d \t %.16e \n", k, ggb[k]);

  if (write_output)
      for (k = 0; k < N; k++)
           fprintf (fpout, "%d \t %.16e \n", k, ggb[k]);

  if (time)
      fprintf (stderr, "\nGTJACOBI: Elapsed Time = %ld millisec. [%.3lf secs] \n\n", mtime, (double) mtime / 1000.0);

  if (write_output)
      fprintf (stderr, "GTJACOBI: Gegenbauer coefficients with Lambda = %.2lf saved in file --> %s \n\n", lambda, filename_out);
}
/* ---------------------------------------------------------------------- */

static void
error (char *message)
{
  fprintf (stderr, "\n Error: %s\n", message);
  fprintf (stderr, " Usage: %s [common output filename] \n", progname);
  fprintf (stderr, "        [-j order (Def.: 20)] [-l lambda (Def.: 0.5)]\n");
  fprintf (stderr, "        [-n number of computed Gegenbauer coefficients (Def.: 128)] \n");
  fprintf (stderr, "        [-1 parameter#1 to feed the inputfunction] ... [-4 parameter#4 to feed inpufunction] \n");
  fprintf (stderr, "        [-v for no verbose][-t for not showing computing time values] \n\n");
  exit (-1);
}
/* -------------------------------------------------------------------------------- */
