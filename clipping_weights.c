int set_clippingweights(){
  if(appliedClippingThreshold >= 1000){
    for(j=0; j<Vipers_Num; j++)  clip_galweight[j] = 1.0;

    return 0;
  }

  else if((appliedClippingThreshold < 1000) && (Jenkins_foldfactor > 1.0)){
    load_clippingweights();  // clip folded measurements using precomputed (foldfactor == 1) weights.

    return 0;
  }

  else{       
    // ** calculation of clipping weights, must have no folding. **//  
    double                       chi;
    double*             gal_occupied;
    double*             cell_weights;
    double*            rand_occupied;

    double     number_occupied = 0.0;

     gal_occupied = malloc(n0*n1*n2*sizeof(*rand_occupied));
    rand_occupied = malloc(n0*n1*n2*sizeof(*rand_occupied));
    cell_weights  = malloc(n0*n1*n2*sizeof(*cell_weights));
    
    for(j=0; j<n0*n1*n2; j++){
      overdensity[j][0] = 0.0;
      overdensity[j][1] = 0.0;

       gal_occupied[j]  = 0.0;
      rand_occupied[j]  = 0.0;
    }

    for(j=0; j<n0*n1*n2; j++)  cell_weights[j] = 1.0;  // initialise to no clipping.
    
    for(j=0; j<rand_number; j++){
        boxlabel                 =  boxCoordinates(rand_x, rand_y, rand_z, j);

        rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
    }
    
    for(j=0; j<n0*n1*n2; j++){
      if(rand_occupied[j] == 1.){
        number_occupied += 1.0;  // Number of cells over which <1+ delta> is averaged.
      }
    }

    printf("\n\nvolume available: %.2le (h^-1 Mpc)^3.", CellVolume*number_occupied);

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){ 
	    chi                        =  interp_comovingDistance(zobs[j]);
    
            boxlabel                   =  boxCoordinates(xCoor, yCoor, zCoor, j);

	    gal_occupied[boxlabel]     =  1; // assign galaxies using NGP.

            overdensity[boxlabel][0]  +=  pow((*pt2nz)(chi)*sampling[j], -1.); // Galaxy counts. 
      }
    }

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]  /=  CellVolume;                             

    // Smooth n/<n>.  True/false flag for zero mean; apodises boundaries. 
    Gaussian_filter(clipping_smoothing_radius, 0);
    
    // Scale (1 + delta) such that <1+ delta> = 1.; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0; 
    double norm = 0.0;

    for(j=0; j<n0*n1*n2; j++){
      // even if randoms don't fully sample the field then still a volume average over a representative sample.  
      if(rand_occupied[j] > 0.){
	norm += overdensity[j][0];
      }
    }

    norm /= number_occupied;

    printf("\n\nMean renormalisation in clipping weights: %.2lf", norm);

    // Shouldn't need overdensity again.
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]  -=         1.0; // no homogeneous constraint. 
    
    // rescaling of <1+delta>.  by rescaling (1+d), don't slip below d=-1. 
    for(j=0; j<n0*n1*n2; j++){
      smooth_overdensity[j][0]  /=  norm;  
        
      smooth_overdensity[j][0]  -=   1.0; // now it's delta. 

      if(rand_occupied[j] == 0.0){
	       overdensity[j][0] = 0.0;
	smooth_overdensity[j][0] = 0.0;   // restore mask after smoothing. spurious empty cells where randoms should be?  
      }
    
      if(smooth_overdensity[j][0] > appliedClippingThreshold){  // either smooth_overdensity or overdensity. 	    	
        cell_weights[j]  = (1. + appliedClippingThreshold)/(1. + overdensity[j][0]);

        if(rand_occupied[j] > 0.)  fraction_clipped  += 1.0;    // clipped volume fraction.
      }
    }
    
    fraction_clipped /= number_occupied;
    
    printf("\n\napplied clipping threshold: %lf, fraction of cells clipped: %lf", appliedClippingThreshold, fraction_clipped);

    int count = 0;

    if(data_mock_flag == 0){
      sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount); 
    }

    if(data_mock_flag == 1){
      sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
    }

    output = fopen(filepath, "w");
    
    //  Clipping weights shouldn't be renormalised, multiplies (fkp_weight/sampling). 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j]  == true){ 
            boxlabel           = boxCoordinates(xCoor, yCoor, zCoor, j);

            clip_galweight[j]  =                 cell_weights[boxlabel];
        
	    if(clip_galweight[j] > 1.)  printf("\nspurious weight: %.2lf", clip_galweight[j]);
	}
        
        else{clip_galweight[j] = 0.0;}
    
        fprintf(output, "%le \n", clip_galweight[j]);
    }

    fclose(output);
  }
  
  return 0;
}


int load_clippingweights(){
  int line_no;
    
  if(data_mock_flag == 0){
    sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
  }
    
  if(data_mock_flag == 1){
    sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);      
  }

  inputfile = fopen(filepath, "r");  
    
  line_count(inputfile, &line_no);
    
  for(j=0; j<line_no; j++)  fscanf(inputfile, "%le \n", &clip_galweight[j]);
    
  fclose(inputfile);
    
  return 0;
}
