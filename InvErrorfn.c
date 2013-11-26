/* inverse cerf from Fred Schwab 

  The following approximations to the inverse of the error function are
  taken from J. M. Blair, C. A. Edwards, and J. H. Johnson, "Rational
  Chebyshev Approximations for the Inverse of the Error Function",
  Mathematics of Computation, 30 (1976) 827-830 + microfiche appendix. 

*/


#define MAXDUB (1.0e10000) /* made this up , is there a macro? */

double sign(double x)
{
  if( x >= 0.0 )
    return(1.0);
  else 
    return(-1.0);
}


double inverse_erf( double x )
{
  double ax, t, ret;

  ax = fabs(x);

/* This approximation, taken from Table 10 of Blair et al., is valid
   for |x|<=0.75 and has a maximum relative error of 4.47 x 10^-8. */

  if( ax <= 0.75 ) {

    double p[] = {-13.0959967422,26.785225760,-9.289057635};
    double q[] = {-12.0749426297,30.960614529,-17.149977991,1.00000000};

    t = x*x-0.75*0.75;
    ret = x*(p[0]+t*(p[1]+t*p[2]))/(q[0]+t*(q[1]+t*(q[2]+t*q[3])));

  } else if( ax >= 0.75 && ax <= 0.9375 ) {
    double p[] = {-.12402565221,1.0688059574,-1.9594556078,.4230581357};
    double q[] = {-.08827697997,.8900743359,-2.1757031196,1.0000000000};

/* This approximation, taken from Table 29 of Blair et al., is valid
   for .75<=|x|<=.9375 and has a maximum relative error of 4.17 x 10^-8. */

    t = x*x - 0.9375*0.9375;
    ret = x*(p[0]+t*(p[1]+t*(p[2]+t*p[3])))/
         (q[0]+t*(q[1]+t*(q[2]+t*q[3])));

  } else if( ax >= 0.9375 && ax <= (1.0-1.0e-100) ) {
    double p[] = {.1550470003116,1.382719649631,.690969348887,
       -1.128081391617, .680544246825,-.16444156791};
    double q[] = {.155024849822,1.385228141995,1.000000000000};

/* This approximation, taken from Table 50 of Blair et al., is valid
   for .9375<=|x|<=1-10^-100 and has a maximum relative error of 2.45 x 10^-8. */

    t=1.0/sqrt(-log(1.0-ax));
    ret = sign(x)*(p[0]/t+p[1]+t*(p[2]+t*(p[3]+t*(p[4]+t*p[5]))))/
          (q[0]+t*(q[1]+t*(q[2])));
    } else {
      ret = MAXDUB;
    }

   return(ret);
}


int InvErrorfnTest(){
    sprintf(filepath, "%s/Data/ClippedPk/InvErrorfnTest.dat", root_dir);
    output = fopen(filepath, "w");
   
    float x;
    float y;
    
    for(j=0; j<1000; j++){
         x = -1. + pow(1000., -1.)*2.*j;
         y = inverse_erf(x);
    
        fprintf(output, "%f \t %f \n", x, y);
    }

    fclose(output);

    return 0;
}
