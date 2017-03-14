int init_gsl_randgen(){
  gsl_rng_env_setup();  // random generation setup.

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);  // returns a pointer to a newly-created instance of a random number generator of type T        

  return 0;
}


int init_padding(){
  n0 = n1 = n2 = fft_size;
    
  AxisLimsArray[0][0]       =        0.0;      // Embedding volume for P(k) measurement. Stefano basis.
  AxisLimsArray[1][0]       =      800.0;

  AxisLimsArray[0][1]       =        0.0;
  AxisLimsArray[1][1]       =      800.0;      // h^-1 Mpc

  AxisLimsArray[0][2]       =        0.0;      // h^-1 Mpc
  AxisLimsArray[1][2]       =      800.0;
  
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


int initi_dist_z(){
  comovDistReshiftCalc(); // Cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h

  loChi                 = interp_comovingDistance(lo_zlim);
  hiChi                 = interp_comovingDistance(hi_zlim);

  printf("\n\nRedshift limits, lower bound: %.4lf \t %.4lf h^-1 Mpc, \n\t\t upper bound: %.4lf \t %.4lf h^-1 Mpc", lo_zlim, loChi, hi_zlim, hiChi);
  
  return 0;
}


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

  initi_dist_z();

  init_padding();
    
  init_cell_info();

  init_fftgrid();
    
  return 0;
}


double SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval){
    double SolidAngle    = 0.0;
    
    double thetaLowerRad = pi/2. - decUpperBound*pi/180.;
    double thetaUpperRad = pi/2. - decLowerBound*pi/180.;
    
    SolidAngle = (raInterval*pi/180.)*(-cos(thetaUpperRad) + cos(thetaLowerRad));

    printf("\nSolid angle observed in a given mock: %e steradians.", SolidAngle);
    printf("\n                                    : %e sq degrees.", steradians2sqdegs(SolidAngle));
 
    return SolidAngle;
}
