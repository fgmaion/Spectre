/*
  	This program calls the routine ggb_inversion_fft() to compute the Gegenbauer sum at selected points. 

        Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
                   by Enrico De Micheli and Giovanni Alberto Viano.
                   Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008

	Input: File with the Gegenbauer coefficients f_n^lambda, in the form: n1 f_n1^lambda, n2 f_n2^lambda, ...
        output: File with the values of the Gegenbauer sum at: x1 f(x1), x2 f(x2), ...
   
        COMPILE WITH: cc call_ggb_inversion_fft.c -O2 ggb_inversion_fft.o -o call_ggb_inversion_fft -lgsl -lgslcblas -lm

	Author: Enrico De Micheli

        Licensing: Free to use giving credit to the author
        
        This program is distributed WITHOUT ANY WARRANTY. 

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

static double	*double_vector ();
static double   *read_Gegenbauer_coefficients ();
static void 	error ();

char	progname[64];
// ------------------------------------------------------------------------------

main(argc, argv)
int	  argc;
char	**argv;
{
register int	k;
FILE		*fpin = NULL, *fpout = NULL;
char		output_filename[64], input_ggb_filename[64], err[64], common[64];
int		flag_input = 0, flag_common = 0, flag_lambda = 0, zeropadding = 0, time = 1, verbose = 1, plottype = 0;
double		lambda = 0.5, lambda_from_file = -1.0, delta;
double		*f_lambda; 	// Contains the input Gegenbauer coefficients
double		*abscissa, 	// Abscissae of the computed samples of Gegenbauer sum
		*f;		// Contains the output value of the Gegenbauer sum
int		n_samples = 32, // Function computed at n_samples points
		quad = 5, 	// quad = 2 is trapezoidal, quad = 3 is simpson, quad = 4 is simpson38, quad = 5 is boole 
		nmb_ggb_coeff;
struct timeval  start, end;
long            mtime, seconds, useconds;
int		option = 0; 
double		*read_Gegenbauer_coefficients ();

  strcpy (progname, argv[0]);

  opterr = 0;
  while ((option = getopt(argc, argv,":vtpq:z:n:l:o:")) != -1) {
        switch (option) {
             case 'v' : verbose = 0;
                        break;
             case 't' : time = 0;
                        break;
             case 'p' : plottype = 1;	// plottype = 0 -> f(cos u); plottype = 1 -> f(x)
                        break;
             case 'q' : quad = atoi(optarg);
			if (quad < 2 || quad > 5)
			    error ("Unknown quadrature rule");
			break;
             case 'z' : zeropadding = atoi(optarg); 
			if (zeropadding < 0)
			    error ("Zeropadding must be positive");
                        break;
             case 'n' : n_samples = atoi(optarg);
			if (n_samples < 2)
			    error ("n_samples must be >= 2");
                        break;
             case 'l' : lambda = atof(optarg);
			if (lambda <=0)
			    error ("Lambda must be positive");
			flag_lambda = 1;
                        break;
             case 'o' : strcpy (common, optarg);
			flag_common++;
			break;
             case ':' : sprintf (err, "%s - Error: Option: `-%c' needs a value \n", progname, optopt);
                        error (err);
     	     case '?' : sprintf (err, "%s - Error: No such option: `-%c' \n", progname, optopt);
			error (err);
        }
  }

  if (optind == argc)
      error ("Missing input Gegenbauer coefficients filename");
  else
      if (argc-optind != 1)
          error ("Cannot parse correctly the input Gegenbauer coefficients filename");
      else
          strcpy (input_ggb_filename, argv[optind]);

  if (plottype){
      strcpy (output_filename, "IGTFFT_x.dat");
      delta = 2.0/(n_samples-1);
  }
  else {
      strcpy (output_filename, "IGTFFT_cos_u.dat");
      delta = M_PI/(n_samples-1);
  }

  if (flag_common) {
      if (plottype)
          sprintf (output_filename, "%s_IGTFFT_x.dat", common);
      else
          sprintf (output_filename, "%s_IGTFFT_cos_u.dat", common);
  }

  fpin = fopen(input_ggb_filename, "r");
  if (fpin == NULL)
      error ("Can't open input Gegenbauer coefficients file");
  f_lambda = read_Gegenbauer_coefficients (fpin, &nmb_ggb_coeff, &lambda_from_file);
  if (nmb_ggb_coeff & (nmb_ggb_coeff - 1) != 0)
      error ("The number of input Gegenbauer coefficients is not a power of 2");
  fclose (fpin);

  if (flag_lambda == 0)
      if(lambda_from_file == -1.0) // No value of lambda read from file and no lambda value provided in input
         fprintf (stderr, "\nIGTFFT Warning: lambda set to the default value: lambda = 0.5\n");
      else {  // lambda read from file as: #@lambdavalue at the beginning of the Gegenbauer coefficients file
         lambda = lambda_from_file;
         fprintf (stderr, "\nIGTFFT Warning: Value of lambda read from input file: lambda = %lf \n\n", lambda);
      }

  if (verbose)
      fprintf (stderr, " IGTFFT: Read %d Gegenbauer Coefficients from file %s, zero padded up to %d \n", 
               nmb_ggb_coeff, input_ggb_filename, (1<<zeropadding)*nmb_ggb_coeff);

// Sampling of the u-interval MUST NOT be the same as the sampling of the integrating t-interval: u_k != t_j for all (suitable) k and j
  if (n_samples == (1<<zeropadding)*nmb_ggb_coeff+1)
      error ("Choose a different number of sampling points: A sampling point coincides with the pole of the integrand");

  if (time)
      gettimeofday(&start, NULL);

  f = double_vector (n_samples);	// Can contain either f(cos u) or f(x), according to the value of plottype
  abscissa = double_vector (n_samples);	// Can contain either u or x = acos(u), according to the value of plottype

  if(ggb_inversion_fft (f_lambda, nmb_ggb_coeff, lambda, quad, plottype, zeropadding, abscissa, f, n_samples))
     exit (1);

  if (time) {
      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
  }

  if (verbose) {
      fprintf (stderr, " IGTFFT: Output function computed by using %d Gegenbauer coefficients with lambda = %lf \n", nmb_ggb_coeff, lambda);
      fprintf (stderr, " IGTFFT: Reconstructed output function sampled in %d points \n", n_samples-2);
  }

// Save sample values
  fpout = fopen(output_filename, "w");
  if (fpout == NULL)
      error ("Can't open output file");

  fprintf (fpout, "# Inverse Gegenbauer Transform \n");
  fprintf (fpout, "# Data file generated by the program: %s \n", progname);
  fprintf (fpout, "# Reconstruction made with %d Gegenbauer coefficients with lambda = %lf \n", nmb_ggb_coeff, lambda);
  fprintf (fpout, "# Input Gegenbauer coefficients zero padded up to %d \n", nmb_ggb_coeff*(1<<zeropadding));
  fprintf (fpout, "# Gegenbauer coefficients read from file: %s \n#\n", input_ggb_filename);
  if (plottype) {
      fprintf (fpout, "# Plot of the function f(x) versus x in [-1+D,1-D] with D = %lf \n", delta);
      fprintf (fpout, "#x \t\t f(x) \n");
  }
  else {
      fprintf (fpout, "# Plot of the function f(cos(u)) versus u in [D,Pi-D] with D = %lf \n", delta);
      fprintf (fpout, "#u \t\t f(cos(u)) \n");
  }

  for (k = 0; k < n_samples-2; k++)
       fprintf (fpout, "%lf \t %lf \n", abscissa[k], f[k]);
  fclose (fpout);

  if (verbose) 
      if (plottype)
          fprintf (stderr, " IGTFFT: Plot of f(x) for x in [-1+D,1-D] with D = %lf saved in file --> %s \n", delta, output_filename);
      else
          fprintf (stderr, " IGTFFT: Plot of f(cos(u)) for u in [D,Pi-D] with D = %lf saved in file --> %s \n", delta, output_filename);

  if (time)
      fprintf (stderr, "\n Elapsed Time = %ld millisec. [%lf secs] \n\n", mtime, (double) mtime / 1000.0);
}
// --------------------------------------------------------------------------------

static void 
error (message)
char	*message;
{
  fprintf (stderr, "\n Error: %s\n", message);
  fprintf (stderr, " Usage: %s gegenbauer_coefficients_filename \n", progname);
  fprintf (stderr, "        [-l lambda (Def.: 0.5)][-z zero_padding factor (Def.: 0 - no padding)] \n");
  fprintf (stderr, "        [-n number of output samples for the plot (Def: 32)] \n");
  fprintf (stderr, "        [-o common output filename (Def.: output_filename: IGTFFT_cos_u.dat)] \n");
  fprintf (stderr, "        [-q quadrature_rule (quadrature_rule = 2,3,4,5 points)(Def.: 5)] \n");
  fprintf (stderr, "        [-p to switch to plot f(x) (instead of f(cos u))] \n");
  fprintf (stderr, "        [-v for no verbose][-t for not showing computing time values] \n\n");
  exit (-1);
}
// --------------------------------------------------------------------------------

#define	SAVEBUF	256

// Input data are expected to be in format:
//	n1	f_n1
//	n2	f_n2
//	etc

static double *
read_Gegenbauer_coefficients (fp, n_elements, lambda)
FILE	*fp;
int	*n_elements;
double	*lambda;
{
register int	nmb_entry = 0;
double		*c1, *c2; 
char		buf[256];
int		n, size = SAVEBUF;
double		valr;

  c1 = double_vector (size); 

  for (;;) {                            // Skip initial comments #, if any 
       if (fgets (buf, 256, fp) == NULL) error ("Error reading input Gegebauer coefficients");
       if (buf[0] == '#')
           if (buf[1] == '@') // It is the string #@, read lambda
               sscanf (buf+2, "%lf", lambda);
           else
               continue;
       else {
          sscanf (buf, "%d %lf", &n, &valr);	// Read first datum
          c1[n] = valr;	
          nmb_entry++;
          break;
       }
  }

  while (fgets (buf, 256, fp) != NULL) {
         sscanf (buf, "%d %lf", &n, &valr);
         c1[n] = valr;
         nmb_entry++;
         if (nmb_entry >= size) {
             size += SAVEBUF;
	     c2 = double_vector (size);
             memcpy (c2, c1, nmb_entry * sizeof(double));
             free ((char *) c1);
             c1 = c2;
             c2 = NULL;
         }
  }
  *n_elements = nmb_entry;

  return (c1);
}
// ---------------------------------------------------------------------- 

static double 
*double_vector (n)
int     n;
{
double  *v;

  v = (double *) malloc ((unsigned) n*sizeof(double));
  if (!v) 
      error ("Allocation failure in double_vector ()");

  return v;
}
// ---------------------------------------------------------------------- 
