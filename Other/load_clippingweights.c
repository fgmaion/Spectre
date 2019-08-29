int prep_clipping_calc(){
  // iplan            = fftw_plan_dft_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);
  iplan               = fftw_plan_dft_c2r_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_ESTIMATE);

  cell_metd0          = (double *)  calloc(n0*n1*n2, sizeof(double));  // highest d0 met by cell; could be integer.

  filter_factors      = (double *)  calloc(n0*n1*nx, sizeof(double));  // pre-compute (Fourier-space) Gaussian filter factors.

  set_randoccupied();

  prep_filterfactors();  // dx, dy, dz set by rand_occupied.

  return 0;
}

int oldload_clippingweights(){
  int line_no;

  if(data_mock_flag == 0){
    sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, (double) d0, lo_zlim, hi_zlim, loopCount);
  }

  if(data_mock_flag == 1){
    sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, (double) d0, lo_zlim, hi_zlim, loopCount);
  }

  inputfile = fopen(filepath, "r");

  line_count(inputfile, &line_no);

  for(j=0; j<line_no; j++)  fscanf(inputfile, "%le \n", &clip_galweight[j]);

  fclose(inputfile);

  printf("\n\nClipping weights: %s loaded", filepath);

  return 0;
}

int load_clippingweights(){
  int         line_no;

  // default ordering: double d0s[4] = {4., 6., 10., 1000.};
  if(data_mock_flag == 0)  sprintf(filepath, "%s/mocks_v1.7/clip_weights/W%d/mock_%03d_z_%.1lf_%.1lf_%d.dat", outputdir, fieldFlag, loopCount, lo_zlim, hi_zlim, fft_size);
  if(data_mock_flag == 1)  sprintf(filepath, "%s/data_v1.7/clip_weights/W%d/data_%.1lf_z_%.1lf_%d.dat",  outputdir, fieldFlag, lo_zlim, hi_zlim, fft_size);

  inputfile = fopen(filepath, "r");

  line_count(inputfile, &line_no);

  for(j=0; j<line_no; j++){
    if(d0 ==    2)  fscanf(inputfile, "%lf \t %*f \t %*f \t %*f \t %*f \t %*f", &clip_galweight[j]);
    if(d0 ==    4)  fscanf(inputfile, "%*f \t %lf \t %*f \t %*f \t %*f \t %*f", &clip_galweight[j]);
    if(d0 ==    6)  fscanf(inputfile, "%*f \t %*f \t %lf \t %*f \t %*f \t %*f", &clip_galweight[j]);
    if(d0 ==   10)  fscanf(inputfile, "%*f \t %*f \t %*f \t %lf \t %*f \t %*f", &clip_galweight[j]);
    if(d0 == 1000)  fscanf(inputfile, "%*f \t %*f \t %*f \t %*f \t %lf \t %*f", &clip_galweight[j]);
  }

  fclose(inputfile);

  // for(j=0; j<130; j++)  printf("\n%d \t %.6lf", j, clip_galweight[j]);

  return 0;
}

int set_clipping_weights(){
  if(d0 >= 1000){
    for(j=0; j<Vipers_Num; j++)  clip_galweight[j] = 1.0;

    return 0;
  }

  else{
    load_clippingweights();
  }

  return 0;
}
