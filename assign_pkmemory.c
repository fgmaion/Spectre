int prep_c2c(){
  overdensity         = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(*overdensity));
  smooth_overdensity  = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(*overdensity)); // 1d: N input elements -> N/2 + 1 output elements.
  
  H_k                 = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(*H_k));

   plan               = fftw_plan_dft_3d(n0, n1, n2, overdensity, H_k, FFTW_FORWARD, FFTW_ESTIMATE); // FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, //FFTW_EXHAUSTIVE.

  return 0;
}


int prep_r2c(){
  overdensity         = (double*)       fftw_malloc(n0*n1*n2*sizeof(*overdensity));
  smooth_overdensity  = (double*)       fftw_malloc(n0*n1*n2*sizeof(*smooth_overdensity));
  
  H_k                 = (fftw_complex*) fftw_malloc((n2/2 + 1)*n1*n0*sizeof(*H_k)); // returns half the array, along the fastest memory change direction: x.

   plan               = fftw_plan_dft_r2c_3d(n0, n1, n2, overdensity, H_k, FFTW_ESTIMATE); // r2c is always forward. // FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE.
  
  return 0;
}


int prep_x2c(){
  fftw_import_wisdom_from_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  // prep_c2c();
  prep_r2c();

  fftw_export_wisdom_to_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  for(j=0; j<n0*n1*n2; j++)  overdensity[j] = 0.0;

  walltime("Wall time after array malloc");
  
  return 0;
}


int set_kbins_logkspacing(){
  logk_interval        = (logk_max - logk_min)/kbin_no;

  logk_limits          = (double *) malloc(kbin_no*sizeof(*logk_limits));

  for(j=0; j<kbin_no; j++)  logk_limits[j] = pow(10., logk_min + j*logk_interval);

  return 0;
}


int set_kbins_linearkspacing(){
  logk_limits          = (double *) malloc(kbin_no*sizeof(*logk_limits));

  double kmedia;

  kmedia               = 2.*pi/800.;

  for(j=0; j<kbin_no; j++)  logk_limits[j] = j*kmedia; // needs 800 bins.

  return 0;
}


int prep_pkRegression(){
  set_kbins_logkspacing();
  // set_kbins_linearkspacing();

  modes_perbin         = (int  *)    malloc(kbin_no*sizeof(*modes_perbin));

  mean_modk            = (double  *) malloc(kbin_no*sizeof(*mean_modk));
  binnedPk             = (double  *) malloc(kbin_no*sizeof(*binnedPk));
  Monopole             = (double  *) malloc(kbin_no*sizeof(*Monopole));
  Quadrupole           = (double  *) malloc(kbin_no*sizeof(*Quadrupole));
  Hexadecapole         = (double  *) malloc(kbin_no*sizeof(*Hexadecapole));
  // linearErrors      = (double  *) malloc(kbin_no*sizeof(*linearErrors));

  Sum_Li               = malloc(kbin_no*sizeof(double));  // Precompute coefficients for Legendre decomposition.
  Sum_Li2              = malloc(kbin_no*sizeof(double));

  Sum_Pi               = malloc(kbin_no*sizeof(double));
  Sum_PiLi             = malloc(kbin_no*sizeof(double));

  detA                 = malloc(kbin_no*sizeof(double));

  // index where fftw mode sits in pk binning array.
  kind                 = malloc(num_modes*sizeof(int));
  kLi                  = malloc(num_modes*sizeof(double));
  kM2                  = malloc(num_modes*sizeof(double));  // NGP/CIC square correction.

  for(k=0; k<kbin_no; k++){
    mean_modk[k]      = 0.0;

    Monopole[k]       = 0.0;
    Quadrupole[k]     = 0.0;
    modes_perbin[k]   =   0;

    Sum_Li[k]         = 0.0;
    Sum_Li2[k]        = 0.0;

    Sum_Pi[k]         = 0.0;
    Sum_PiLi[k]       = 0.0;
  }

  return 0;
}
