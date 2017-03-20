int prep_c2c(){
  nx = n2;

  overdensity         = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(*overdensity));
  smooth_overdensity  = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(*overdensity)); // 1d: N input elements -> N/2 + 1 output elements.
  
  H_k                 = (fftw_complex*) fftw_malloc(n0*n1*nx*sizeof(*H_k));

  plan                = fftw_plan_dft_3d(n0, n1, n2, overdensity, H_k, FFTW_FORWARD, FFTW_ESTIMATE); // FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, //FFTW_EXHAUSTIVE.

  return 0;
}


int prep_r2c(){
  nx = n2/2 + 1;

  overdensity         = (double*)       fftw_malloc(n0*n1*n2*sizeof(*overdensity));
  smooth_overdensity  = (double*)       fftw_malloc(n0*n1*n2*sizeof(*smooth_overdensity));
  
  H_k                 = (fftw_complex*) fftw_malloc(n0*n1*nx*sizeof(*H_k)); // returns half the array, along the fastest memory change direction: x.

  // FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE.
  plan                = fftw_plan_dft_r2c_3d(n0, n1, n2, overdensity, H_k, FFTW_MEASURE); // r2c is always forward.
  
  return 0;
}


int prep_x2c(){
  fftw_import_wisdom_from_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  //prep_c2c();
  prep_r2c();

  fftw_export_wisdom_to_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  for(j=0; j<n0*n1*n2; j++)  overdensity[j] = 0.0;

  walltime("Wall time after array malloc");
  
  return 0;
}


int set_kbins_logkspacing(regress* inst){
  logk_interval        = (logk_max - logk_min)/KBIN_NO;

  for(j=0; j<KBIN_NO; j++)  inst->logk_limits[j] = pow(10., logk_min + j*logk_interval);

  return 0;
}


int set_kbins_linearkspacing(regress* inst){
  double kmedia;

  kmedia               = 2.*pi/800.;

  for(j=0; j<KBIN_NO; j++)  inst->logk_limits[j] = j*kmedia; // needs 800 bins.

  return 0;
}


int regress_mem(regress* inst){
  inst->kind = calloc(n0*n1*nx, sizeof(int));

  inst->kLi  = calloc(n0*n1*nx, sizeof(double));
  inst->kM2  = calloc(n0*n1*nx, sizeof(double));

  return 0;
}


int prep_pkRegression(){
  regress_mem(&flat);
  regress_mem(&half);
  
  set_kbins_logkspacing(&flat);
  set_kbins_logkspacing(&half);
    
  return 0;
}
