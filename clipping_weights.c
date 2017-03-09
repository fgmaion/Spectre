int prep_clipping_calc(){
  if((d0 <= 1000) && (foldfactor == 1.0)){
    // iplan            = fftw_plan_dft_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);
    iplan               = fftw_plan_dft_c2r_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_ESTIMATE);

    cell_metd0          = (double *)  malloc(n0*n1*n2*sizeof(double)); // highest d0 met by cell; could be integer.
  
    filter_factors      = (double *)  malloc(n0*n1*nx*sizeof(double));  // pre-compute (Fourier-space) Gaussian filter factors. 

    prep_filterfactors();  // dx, dy, dz set by rand_occupied.

    set_randoccupied();
  }
  
  return 0;
}


int get_clipping_weights(){
 if(d0 >= 1000){
   for(j=0; j<Vipers_Num; j++)  clip_galweight[j] = 1.0;

   return 0;
 }

 else if(foldfactor > 1.0){
   load_clippingweights();  // %% calculation of clipping weights, must have no folding. %% 

   return 0;
 }

 else{
   walltime("\n\nStarting clipping calc. at");
   
  for(j=0; j<n0*n1*n2; j++) overdensity[j] = 0.0; // for each mock.
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + n2*ylabel + n2*n1*zlabel;

      overdensity[boxlabel]  +=  pow(dx*dy*dz*interp_nz(rDist[j])*sampling[j], -1.); // N/<N>
    }
  }

  printf("\nx min:  %.3f \t x max:  %.3f", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
  printf("\ny min:  %.3f \t y max:  %.3f", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
  printf("\nz min:  %.3f \t z max:  %.3f", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));
  
  // Smooth N/<N>.
  Gaussian_filter(); // assumes m0 = n0 etc. 
  
  // Scale (1 + delta) such that <1+ delta> = 1. within surveyed region; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0;
  double norm, snorm, frac_clip;

  norm = snorm = frac_clip = 0.0;

  // By rescaling (1 + smooth delta), don't slip below delta = -1.
  for(j=0; j<n0*n1*n2; j++){
    smooth_overdensity[j]  *=      rand_occupied[j]; // smoothing fills in gaps.

    norm                   +=        overdensity[j];
    snorm                  += smooth_overdensity[j];
  }

   norm /= number_occupied;
  snorm /= number_occupied;

  printf("\n\nMean renormalisation of delta: %.6lf, smooth delta: %.6lf", norm, snorm);
  
  // H_k was used in Gaussian filter. reassign to be used to hold clipping weights for each cell. 
  for(j=0; j<n0*n1*n2; j++)      cell_metd0[j] = 0.0; // Holds the minimum d0 satisfied by cell (d_i > d0) if inside survey; 0 otherwise. 

  for(j=0; j<number_occupied; j++){
                       i    =  occupied_indices[j];

           overdensity[i]  /=   norm;
    smooth_overdensity[i]  /=  snorm;

           overdensity[i]  -=   1.0;
    smooth_overdensity[i]  -=   1.0; // now it's delta.
  }

  double td0; // this d0.
  double d0s[4] = {4., 6., 10., 1000.};

  int number_clipped;

  number_clipped = number_occupied;

  printf("\n\nThresholds volume calc.");
  
  for(k=0; k<4; k++){
    Index     =      0;
    frac_clip =    0.0;

    td0       = d0s[k];
    
    for(j=0; j<number_clipped; j++){
      i = occupied_indices[j];

      if(smooth_overdensity[i] >= td0){
        cell_metd0[i] = td0; // either smooth_overdensity or overdensity.
                             // overdensity has not had <1 + delta> = 1 enforced in master. 

        occupied_indices[Index] = i;

        Index     +=   1;  

        frac_clip += 1.0;
      }
    }

    number_clipped = Index;

    frac_clip /= number_occupied;

    printf("\nFor d0 of %.4lf, %lf%% of cells are clipped", td0, 100.*frac_clip);
  }
  
  if(data_mock_flag == 0){
    sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/wghts_z_%.1lf_%.1lf_%d.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, loopCount);
  }
  
  if(data_mock_flag == 1){
    sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/wghts_%.1lf_z_%.1lf_%d.dat",  root_dir, fieldFlag, lo_zlim, hi_zlim, loopCount);
  }

  output = fopen(filepath, "w");
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j]  == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + n2*ylabel + n2*n1*zlabel;

      // Ordering: d0s[4] = {4., 6., 10., 1000.};
      for(i=0; i<4; i++){
        if(d0s[i] <= cell_metd0[boxlabel])  fprintf(output, "%.4lf \t", (1. + d0s[i])/(1. + overdensity[boxlabel]));

        else{
          fprintf(output, "%.4lf \t", 1.0);
        }
      }

      fprintf(output, "%.4lf", cell_metd0[boxlabel]);
      // fprintf(output, "%.4lf \t %.4lf \t %.4lf", smooth_overdensity[boxlabel], overdensity[boxlabel], cell_metd0[boxlabel]);
    }

    else{
      for(i=0; i<5; i++)  fprintf(output, "%.4lf \t", 0.0);
    }

    fprintf(output, "\n");
  }

  fclose(output);
 }

 // and now load. 
 load_clippingweights();
 
 return 0;
}


int load_clippingweights(){
  int line_no;

  // default ordering. 
  // double d0s[4] = {4., 6., 10., 1000.};
  
  if(data_mock_flag == 0){
    sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/wghts_z_%.1lf_%.1lf_%d.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, loopCount);
  }

  if(data_mock_flag == 1){
    sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/wghts_%.1lf_z_%.1lf_%d.dat",  root_dir, fieldFlag, lo_zlim, hi_zlim, loopCount);
  }
  
  inputfile = fopen(filepath, "r");  
    
  line_count(inputfile, &line_no);

  for(j=0; j<line_no; j++){
    if(d0 ==    4.)  fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 ==    6.)  fscanf(inputfile, "%*lf \t %lf \t %*lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 ==   10.)  fscanf(inputfile, "%*lf \t %*lf \t %lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 == 1000.)  fscanf(inputfile, "%*lf \t %*lf \t %*lf \t %lf \t %*lf", &clip_galweight[j]);
  }
  
  fclose(inputfile);

  // for(j=0; j<130; j++)  printf("\n%d \t %.6lf", j, clip_galweight[j]);
  
  return 0;
}
