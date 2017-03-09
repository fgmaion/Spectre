int prep_clipping_calc(){
  iplan               = fftw_plan_dft_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);

  cell_cweights       = (double *)      malloc(n0*n1*n2*sizeof(*cell_cweights));

  // pre-compute Gaussian filter factors. 
  filter_factors      = (double *)      malloc(n0*n1*n2*sizeof(*filter_factors));

  prep_filterfactors(dx, dy, dz);
  
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
   
  for(j=0; j<m0*m1*m2; j++) overdensity[j] = 0.0; // for each mock.
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + m2*ylabel + m2*m1*zlabel;

      overdensity[boxlabel]  +=  pow(dx*dy*dz*interp_nz(rDist[j])*sampling[j], -1.); // N/<N>
    }
  }

  printf("\nx min:  %.3f \t x max:  %.3f", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
  printf("\ny min:  %.3f \t y max:  %.3f", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
  printf("\nz min:  %.3f \t z max:  %.3f", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));
  
  // Smooth N/<N>.
  // Gaussian_filter(smooth_radius, dx, dy, dz); // assumes m0 = n0 etc. 
  
  // Scale (1 + delta) such that <1+ delta> = 1. within surveyed region; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0;
  double      norm = 0.0;
  double frac_clip = 0.0;

  // By rescaling (1 + smooth delta), don't slip below delta = -1.
  for(j=0; j<n0*n1*n2; j++){
    smooth_overdensity[j]  *=  rand_occupied[j]; // smoothing fills in gaps.

    norm                   += smooth_overdensity[j];
  }

  norm /= number_occupied;

  printf("\n\nMean renormalisation in clipping weights: %.4lf", norm);

  // H_k was used in Gaussian filter. reassign to be used to hold clipping weights for each cell. 
  for(j=0; j<n0*n1*n2; j++)      cell_cweights[j] = 1.0; // This will hold the clipping weight for the cell (initialise to no clipping)  

  for(j=0; j<number_occupied; j++){
    i                       =  occupied_indices[j];

    smooth_overdensity[i]  /=  norm;

    smooth_overdensity[i]  -=   1.0; // now it's delta.

    if(smooth_overdensity[i] > d0){  // either smooth_overdensity or overdensity.
      cell_cweights[i]      = (1. + d0)/(1. + overdensity[i]); // overdensity has not had <1 + delta> = 1 enforced. 

      frac_clip            += 1.0;    // must be surveyed, otherwise smooth_overdensity would be zero and cell would be bypassed.  
    }
  }

  frac_clip /= number_occupied;

  printf("\n\nd0: %lf, cells clipped: %lf \%", d0, 100.*frac_clip);
  
  for(j=0; j<Vipers_Num; j++){
    clip_galweight[j] = 0.0;

    if(Acceptanceflag[j]  == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + m2*ylabel + m2*m1*zlabel;
      
      clip_galweight[j]  =  cell_cweights[boxlabel];

      if(clip_galweight[j] > 1.)  printf("\nspurious weight: %.2lf", clip_galweight[j]);
    }
  }
 }

 return 0;
}


int load_clippingweights(){
  int line_no;
    
  if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/wghts_d0_%.2lf_z_%.1lf_%.1lf_%d.dat", root_dir, fieldFlag, d0, lo_zlim, hi_zlim, loopCount);  
  if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/wghts_d0_%.2lf_z_%.1lf_%.1lf_%d.dat",  root_dir, fieldFlag, d0, lo_zlim, hi_zlim, loopCount);   

  inputfile = fopen(filepath, "r");  
    
  line_count(inputfile, &line_no);
    
  for(j=0; j<line_no; j++)  fscanf(inputfile, "%le \n", &clip_galweight[j]);
    
  fclose(inputfile);
    
  return 0;
}
