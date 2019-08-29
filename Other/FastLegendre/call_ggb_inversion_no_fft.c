/*
	This program calls the routine ggb_inversion_no_fft() to compute the Gegenbauer sum at selected points.

        Reference: The expansion in Gegenbauer polynomials: a simple method for the computation of the Gegenbauer coefficients
                   by Enrico De Micheli and Giovanni Alberto Viano.
                   Journal of Computational Physics. http://dx.doi.org/10.1016/j.jcp.2013.01.008

        Input: File with the Gegenbauer coefficients f_n^lambda in the form: n1 f_n1^lambda, n2 f_n2^lambda, ...
        output: File with the values of the Gegenbauer sum at: x1 f(x1), x2 f(x2), ...

        COMPILE WITH: cc call_ggb_inversion_no_fft.c -O2 ggb_inversion_no_fft.o -o call_ggb_inversion_no_fft -lgsl -lgslcblas -lm

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
static void	error ();

char    progname[64];

/* ------------------------------------------------------------------------------ */

main(argc, argv)
int	  argc;
char	**argv;
{
register int	k;
FILE		*fpin = NULL, *fpout = NULL;
char		tmp_filename[64], output_filename[64], input_ggb_coeff_filename[64], common[64], err[64];
int		nmb_ggb_coeff, nmb_ggb_coeff_to_use = 32, nmb_out_samples = 32;
int		flag_output = 0, flag_lambda = 0, plottype = 0;
int		option = 0, worksize = 2000, verbose = 1, time = 1;
double		*f_lambda, 		/* Contains the Gegenbauer Coefficients */
		*c_lambda,		/* Contains the transformed Gegenbauer coefficients c_n^lambda */
		lambda = 0.5, tmp_lambda = -1.0, delta;
double		*abscissae, *f;
double		eps = 1.0e-12;    /* Tolerances for integral computation */
struct timeval  start, end;
long            mtime, seconds, useconds;

  strcpy (progname, argv[0]);

  opterr = 0;
  while ((option = getopt(argc, argv,":vtpw:E:e:n:l:o:")) != -1) {
        switch (option) {
             case 'v' : verbose = 0;
                        break;
             case 't' : time = 0;
                        break;
             case 'p' : plottype = 1;   // plottype = 0 -> f(cos u); plottype = 1 -> f(x)
                        break;
             case 'w' : worksize = atoi(optarg);
                        if (worksize < 1000)
                            error ("The parameter WORKSIZE must be > 1000");
                        break;
             case 'e' : eps = atof(optarg);
                        if (eps <= 0)
                            error ("Integration absolute error must be positive");
                        break;
             case 'n' : nmb_out_samples = atoi(optarg);
                        if (nmb_out_samples < 2)
                            error ("nmb_out_samples must be >= 2");
                        break;
             case 'N' : nmb_ggb_coeff_to_use = atoi(optarg);
                        if (nmb_ggb_coeff_to_use < 1)
                            error ("nmb_ggb_coeff_to_use must be >= 1");
                        break;
             case 'l' : lambda = atof(optarg);
                        if (lambda <=0)
                            error ("Lambda must be positive");
			flag_lambda++;
                        break;
             case 'o' : strcpy (common, optarg);
                        flag_output++;
                        break;
             case ':' : fprintf (stderr, "\n %s - IGTNOFFT Error: Option: `-%c' needs a value \n", progname, optopt);
                        error ("");
             case '?' : fprintf (stderr, "\n %s - IGTNOFFT Error: No such option: `-%c' \n", progname, optopt);
                        error ("");
        }
  }

  if (optind == argc)
      error ("IGTNOFFT: Missing input Gegenbauer coefficients filename");
  else
      if (argc-optind != 1)
          error ("Cannot parse correctly the input Gegenbauer coefficients filename");
      else
          strcpy (input_ggb_coeff_filename, argv[optind]);

  if (plottype){
      strcpy (tmp_filename, "IGTNOFFT_x.dat");
      delta = 2.0/(nmb_out_samples-1);
  }
  else{
      strcpy (tmp_filename, "IGTNOFFT_cos_u.dat");
      delta = M_PI/(nmb_out_samples-1);
  }
  if (flag_output)
      sprintf (output_filename, "%s_%s", common, tmp_filename);

  fpin = fopen(input_ggb_coeff_filename, "r");
  if (fpin == NULL)
      error ("Can't open input Gegenbauer coefficients file");
  f_lambda = read_Gegenbauer_coefficients (fpin, &nmb_ggb_coeff, &tmp_lambda);
  fclose (fpin);

  if (flag_lambda == 0)
      if(tmp_lambda == -1.0) // Didn't read any lambda from file
         fprintf (stderr, "\n IGTNOFFT Warning: lambda set to the default value: lambda = 0.5\n\n");
      else {  // lambda read from file
         lambda = tmp_lambda;
         fprintf (stderr, "\n IGTNOFFT Warning: Value of lambda read from input file: lambda = %lf \n\n", lambda);
      }

  if (verbose)
      fprintf (stderr, "\n IGTNOFFT: Read %d Gegenbauer Coefficients from file %s \n", nmb_ggb_coeff, input_ggb_coeff_filename);

  if (nmb_ggb_coeff_to_use > nmb_ggb_coeff)
      nmb_ggb_coeff_to_use = nmb_ggb_coeff;

  if (time)
      gettimeofday(&start, NULL);

  f = double_vector (nmb_out_samples);		// Contains the values of f(cos u) or f(x) at the sampled points
  abscissae = double_vector (nmb_out_samples);	// Contains the u-(or x-)coordinates of the output sampled points
  if (ggb_inversion_no_fft (f_lambda, nmb_ggb_coeff, nmb_ggb_coeff_to_use, lambda, plottype, worksize, eps, abscissae, f, nmb_out_samples))
      exit (1);

  if (time) {
      gettimeofday(&end, NULL);
      seconds  = end.tv_sec  - start.tv_sec;
      useconds = end.tv_usec - start.tv_usec;
      mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
  }

  if (verbose) {
      fprintf (stderr, " IGTNOFFT: Output function computed by using %d Gegenbauer coefficients with lambda = %lf \n", nmb_ggb_coeff_to_use, lambda);
      fprintf (stderr, " IGTNOFFT: Reconstructed output function sampled in %d points \n", nmb_out_samples-2);
  }

// Save output plot file
  fpout = fopen(output_filename, "w");
  if (fpout == NULL)
      error ("Can't open output file");

  fprintf (fpout, "# Inverse Gegenbauer Transform \n");
  fprintf (fpout, "# Data file generated by the program: %s \n", progname);
  fprintf (fpout, "# Reconstruction made with %d Gegenbauer coefficients with lambda = %lf \n", nmb_ggb_coeff_to_use, lambda);
  fprintf (fpout, "# Gegenbauer coefficients read from file: %s \n#\n", input_ggb_coeff_filename);
  if (plottype) {
      fprintf (fpout, "# Plot of the function f(x) versus x in [-1+D,1-D] with D = %lf \n", delta);
      fprintf (fpout, "#x \t\t f(x) \n");
  }
  else {
      fprintf (fpout, "# Plot of the function f(cos(u)) versus u in [D,Pi-D] with D = %lf \n", delta);
      fprintf (fpout, "#u \t\t f(cos(u)) \n");
  }

  for (k = 0; k < nmb_out_samples-2; k++)
       fprintf (fpout, "%lf \t %lf \n", abscissae[k], f[k]);
  fclose (fpout);

  if (verbose)
      if (plottype)
          fprintf (stderr, " IGTNOFFT: Plot of f(x) for x in [-1+D,1-D] with D = %lf saved in file --> %s \n", delta, output_filename);
      else
          fprintf (stderr, " IGTNOFFT: Plot of f(cos(u)) for u in [D,Pi-D] with D = %lf saved in file --> %s \n", delta, output_filename);

  if (time)
      fprintf (stderr, "\n IGTNOFFT: Elapsed Time = %ld millisec. [%lf secs] \n\n", mtime, (double) mtime / 1000.0);
}
/* -------------------------------------------------------------------------------- */

static void 
error (message)
char	*message;
{
  fprintf (stderr, "\n Error: %s\n", message);
  fprintf (stderr, " Usage: %s input_ggb_coeff_filename \n", progname);
  fprintf (stderr, "        [-l lambda (Def: 0.5)][-e abs_err (Def.: 1e-12)][-E rel_err (Def.: 1.0e-12)] \n");
  fprintf (stderr, "        [-N nmb_of_used_gegenbauer_coefficients (Def: 32)] \n");
  fprintf (stderr, "        [-n number of output samples for the plot (Def: 32)][-w integration worksize (Def.: 2000)] \n");
  fprintf (stderr, "        [-o common output filename (Def.: output_filename: IGTNOFFT_{cos_u,x}.dat)] \n");
  fprintf (stderr, "        [-v for no verbose][-t for not showing computing time values] \n");
  fprintf (stderr, "        [-p to switch to plot f(x) (instead of f(cos u))] \n\n");
  exit (-1);
}
/* -------------------------------------------------------------------------------- */

#define SAVEBUF	256

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
double		l, valr;

  c1 = double_vector (size); 

  for (;;) {                            // Skip initial comments #, if any
       if (fgets (buf, 256, fp) == NULL) 
	   error ("Error reading Gegenbauer coefficients");
       if (buf[0] == '#')
           if (buf[1] == '@') // It is lambda!
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
/* ---------------------------------------------------------------------- */

double *double_vector (n)
int     n;
{
double  *v;

  v = (double *) malloc ((unsigned) n*sizeof(double));
  if (!v)
      error ("double_vector: Allocation failure in double_vector ()");
  return v;
}
/* ---------------------------------------------------------------------- */
