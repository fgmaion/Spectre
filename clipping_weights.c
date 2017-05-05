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
  

int calc_clipping_weights(){
  walltime("Starting clipping calc.");
  
  for(j=0; j<n0*n1*n2; j++) overdensity[j] = 0.0; // for each mock.

  walltime("Overdensity calc.");
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + n2*ylabel + n2*n1*zlabel;

      overdensity[boxlabel]  +=  pow(dx*dy*dz*(*pt2nz)(rDist[j])*sampling[j], -1.); // N/<N>
    }
  }

  // printf("\nx min:  %.3f \t x max:  %.3f", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
  // printf("\ny min:  %.3f \t y max:  %.3f", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
  // printf("\nz min:  %.3f \t z max:  %.3f", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));

  walltime("Gaussian filter start.");
  
  // Smooth N/<N>.
  Gaussian_filter(); // assumes m0 = n0 etc. 

  walltime("Gaussian filter end.");

  // Scale (1 + delta) such that <1+ delta> = 1. within surveyed region; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0;
  double norm, snorm, frac_clip;

  norm = snorm = frac_clip = 0.0;

  #pragma omp parallel for reduction(+: norm, snorm) private(j) if(thread ==1)
  for(j=0; j<number_occupied; j++){
     norm +=        overdensity[occupied_indices[j]];
    snorm += smooth_overdensity[occupied_indices[j]];
  }
  
   norm /= number_occupied;
  snorm /= number_occupied;

  printf("\n\nMean renormalisation of delta: %.6lf, smooth delta: %.6lf", norm, snorm);

  #pragma omp parallel for private(i, j) if(thread ==1)
  for(j=0; j<number_occupied; j++){
             i = occupied_indices[j]; // only pick out unmasked cells; no need to reapply mask. 

           overdensity[i]  /=   norm;
    smooth_overdensity[i]  /=  snorm;

           overdensity[i]  -=   1.0;
    smooth_overdensity[i]  -=   1.0; // now it's delta.
  }

  #pragma omp parallel for private(j) if(thread == 1)
  for(j=0; j<n0*n1*n2; j++){
            cell_metd0[j]  = 0.0; // Holds the minimum d0 satisfied by cell (d_i > d0) if inside survey; 0 otherwise.

           overdensity[j] *= rand_occupied[j];
    smooth_overdensity[j] *= rand_occupied[j];
  }

  
  double td0; // this d0.
  double d0s[5] = {2., 4., 6., 10., 1000.};

  double  dc[5];
  
  int number_clipped;

  number_clipped = number_occupied;

  walltime("Thresholds volume calc.");
  
  for(k=0; k<5; k++){
    Index     =      0;
    frac_clip =    0.0;

    td0       = d0s[k];
    dc[k]     =    0.0;
    
    for(j=0; j<number_clipped; j++){
      i = occupied_indices[j];

      if(smooth_overdensity[i] >= td0){
        cell_metd0[i] = td0; // either smooth_overdensity or overdensity.
        // overdensity has not had <1 + delta> = 1 enforced in previous calc. 

        dc[k]     += overdensity[i] - td0;
        
        occupied_indices[Index] = i;

        Index     +=   1;
        frac_clip += 1.0;
      }
    }

    number_clipped = Index;

    frac_clip     /= number_occupied;
    dc[k]         /= number_occupied;

    dc[k]          = pow(dc[k], 2.); // shift in zeroth-mode power.
    
    printf("\nFor d0 of %.4lf, %lf%% of cells are clipped. DC shift in power is %.6lf", td0, 100.*frac_clip, dc[k]);
  }
  
  if(data_mock_flag == 0)  sprintf(filepath, "%s/mocks_v1.7/clip_weights/W%d/mock_%03d_z_%.1lf_%.1lf_%d.dat", outputdir, fieldFlag, loopCount, lo_zlim, hi_zlim, fft_size);
  if(data_mock_flag == 1)  sprintf(filepath, "%s/data_v1.7/clip_weights/W%d/data_%.1lf_z_%.1lf_%d.dat",  outputdir, fieldFlag, lo_zlim, hi_zlim, fft_size);

  output = fopen(filepath, "w");

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j]  == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + n2*ylabel + n2*n1*zlabel;

      // Ordering: d0s[4] = {4., 6., 10., 1000.};
      for(i=0; i<5; i++){
        if(d0s[i] <= cell_metd0[boxlabel])  fprintf(output, "%lf \t", (1. + d0s[i])/(1. + overdensity[boxlabel]));

        else{
          fprintf(output, "%lf \t", 1.0);
        }
      }

      fprintf(output, "%lf", cell_metd0[boxlabel]);
      // fprintf(output, "%.4lf \t %.4lf \t %.4lf", smooth_overdensity[boxlabel], overdensity[boxlabel], cell_metd0[boxlabel]);
    }

    else{
      for(i=0; i<6; i++)  fprintf(output, "%lf \t", 0.0);
    }

    fprintf(output, "\n");
  }

  fclose(output);

  // DC shift file. 
  if(data_mock_flag == 0){
    sprintf(filepath, "%s/mocks_v1.7/dc_shifts/W%d/mock_%.1lf_%.1lf_%d.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, fft_size);

    output = fopen(filepath, "a");

    for(i=0; i<5; i++)  fprintf(output, "%lf \t", dc[i]);

    fprintf(output, "\n");
  
    fclose(output);
  }
  
  return 0;
}


int print_xy(){
  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/overdensity/W%d/mock_%03d_z_%.1lf_%.1lf_xy.dat", root_dir, fieldFlag, loopCount, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true)  fprintf(output, "%.6le \t %.6le \t %.6le \n", xCoor[j], yCoor[j], zCoor[j]);
  }
  
  fclose(output);
  
  return 0;
}


int print_metd0(){
  int array[n0][n1];
  
  printf("\n\nPrinting met d0.");

  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/overdensity/W%d/mock_%03d_z_%.1lf_%.1lf.dat", root_dir, fieldFlag, loopCount, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      array[k][j] = 0;
      
      for(i=201; i<202; i++){
          Index = k*n1*n2 + j*n2 + i;

          if(cell_metd0[Index] ==    0.)  array[k][j] = 0;
          if(cell_metd0[Index] ==    4.)  array[k][j] = 1;
          if(cell_metd0[Index] ==    6.)  array[k][j] = 2;
          if(cell_metd0[Index] ==   10.)  array[k][j] = 3;                                   
          if(cell_metd0[Index] == 1000.)  array[k][j] = 4;
      }

      fprintf(output, "%d \t", array[k][j]);
    }

    fprintf(output, "\n");
  }
      
  fclose(output);
  
  return 0;
}
