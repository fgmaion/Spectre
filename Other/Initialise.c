int init_gsl_randgen(){
  gsl_rng_env_setup();  // random generation setup.

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);  // returns a pointer to a newly-created instance of a random number generator of type T        

  // seed buffers needed for random number generation for each thread. 
  randBuffers = malloc(omp_get_max_threads()*sizeof(*randBuffers));
                       
  for(j=0; j<omp_get_max_threads(); j++)  srand48_r((long) j, &randBuffers[j]);

  return 0;
}

int init_padding(){
  n0 = n1 = n2 = fft_size;
    
  AxisLimsArray[0][0]       =        0.0;      // Embedding volume for P(k) measurement. Stefano basis.
  AxisLimsArray[1][0]       =     1000.0;      // required for 0.7 < z < 1.2; previously 800. 

  AxisLimsArray[0][1]       =        0.0;
  AxisLimsArray[1][1]       =     1000.0;      // h^-1 Mpc

  AxisLimsArray[0][2]       =        0.0;      // h^-1 Mpc
  AxisLimsArray[1][2]       =     1000.0;
  
  return 0;
}

int init_cell_info(){
  // assumes cubic cells.
  xCellSize             = (AxisLimsArray[1][2] - AxisLimsArray[0][2])/n2;
  yCellSize             = (AxisLimsArray[1][1] - AxisLimsArray[0][1])/n1;
  zCellSize             = (AxisLimsArray[1][0] - AxisLimsArray[0][0])/n0;

  // Clipping weights. 
  CellVolume            = xCellSize*yCellSize*zCellSize;                            // h^-3 Mpc^3
  
  return 0;
}
/*
int init_dist_z(){
  chi_zcalc(); // Cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h
  
  printf("\n\nRedshift limits, lower bound: %.6lf\t%.10lf h^-1 Mpc (%.10lf Mpc), \n\t\t upper bound: %.6lf\t%.10lf h^-1 Mpc (%.10lf Mpc)", lo_zlim, loChi, loChi/h, hi_zlim, hiChi, hiChi/h);

  printf("\n\nVolume surveyed by field W%d: %.10lf [10^-3 (h^-1 Gpc)^3]", fieldFlag, pow(10., 3.)*calc_vol());
  
  return 0;
}

double check_radialextent(double lochi, double hichi, double lopad){
  double radial_extent = hichi - lochi;

  printf("\n\nRadial extent: %.3lf [h^-1 Mpc] \n", radial_extent);
  
  if(radial_extent > (lopad + AxisLimsArray[1][0] - AxisLimsArray[0][0])){
    exit(EXIT_FAILURE);
  }
  
  return 0;
}
*/
int init_fftgrid(){
  xNy = pi/xCellSize;  // k = 2*pi x Nyquist frequency
  yNy = pi/yCellSize;  // k = 2*pi x Nyquist frequency
  zNy = pi/zCellSize;  // k = 2*pi x Nyquist frequency

  fund_kx = 2.*pi*pow(n2, -1.)*pow(xCellSize, -1.);
  fund_ky = 2.*pi*pow(n1, -1.)*pow(yCellSize, -1.);
  fund_kz = 2.*pi*pow(n0, -1.)*pow(zCellSize, -1.);
  
  return 0;
}

int Initialise(){
  init_gsl_randgen();

  //  init_dist_z();

  init_padding();
    
  init_cell_info();

  init_fftgrid();
    
  return 0;
}
/*
double SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval){
    double SolidAngle    = 0.0;
    
    double thetaLowerRad = pi/2. - decUpperBound*pi/180.;
    double thetaUpperRad = pi/2. - decLowerBound*pi/180.;
    
    SolidAngle = (raInterval*pi/180.)*(-cos(thetaUpperRad) + cos(thetaLowerRad));

    printf("\nSolid angle observed in a given mock: %e steradians.", SolidAngle);
    printf("\n                                    : %e sq degrees.", steradians2sqdegs(SolidAngle));
 
    return SolidAngle;
}
*/
