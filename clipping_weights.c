// New box with limited padding/no rotation.  
//
//
//

double rd(double arg, double base){
  return base*floor(arg/base);
}


double ru(double arg, double base){
  return base*(1. + floor(arg/base));
}


int calc_occupiedvol(){
  int m2, m1, m0;

  m0 =  n0;
  m1 =  n1;
  m2 =  n2;
  
  // Surveyed volume need only be done once; no rotation required.  
  int     number_occupied = 0;
  int*          rand_occupied;

  rand_occupied = malloc(m2*m1*m0*sizeof(*rand_occupied));

  for(j=0; j<m0*m1*m2; j++)  rand_occupied[j] = 0;

  // Reassign randoms to have constant \bar n(z); this will (much) better sample high z tail.
  pt2nz = &unity; // constant <n(z)>

  prep_inverseCumulative_nbar();

  // New basis for embedding volume.
  double min_x, max_x, dx, min_y, max_y, dy, min_z, max_z, dz, F, cos_dec;

  for(j=0; j<rand_number; j++){
    F           = gsl_rng_uniform(gsl_ran_r);  // Chi limits satisfied by construction.

    rand_chi[j] = inverse_cumulative_nbar(F);

    cos_dec     = cos(rand_dec[j]);

    rand_x[j]   =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
    rand_y[j]   =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
    rand_z[j]   = -rand_chi[j]*sin(rand_dec[j]);
  }

  // what if negative.  
  min_x = rd(arrayMin(rand_x, rand_number), 50.); min_y = rd(arrayMin(rand_y, rand_number), 50.); min_z = rd(arrayMin(rand_z, rand_number), 50.);
  max_x = ru(arrayMax(rand_x, rand_number), 50.); max_y = ru(arrayMax(rand_y, rand_number), 50.); max_z = ru(arrayMax(rand_z, rand_number), 50.);

  dx    = (max_x - min_x)/m2;
  dy    = (max_y - min_y)/m1;
  dz    = (max_z - min_z)/m0;

  printf("\n\nRandoms:");
  printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));

  printf("\n\nNew bounds:");
  printf("\nx: %.1lf \t %.1lf h^-1 Mpc, dx: %.6lf h^-1 Mpc", min_x, max_x, dx);
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc, dy: %.6lf h^-1 Mpc", min_y, max_y, dy);
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc, dz: %.6lf h^-1 Mpc", min_z, max_z, dz);

  for(j=0; j<rand_number; j++){
    xlabel      = (int) floor((rand_x[j] - min_x)/dx);
    ylabel      = (int) floor((rand_y[j] - min_y)/dy);
    zlabel      = (int) floor((rand_z[j] - min_z)/dz);

    boxlabel    = (int)    xlabel + m2*ylabel + m2*m1*zlabel;

    rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
  }

  printf("\n\nSurveyed vol. calc.");
  
  double oldvol, newvol, truvol, convergence;

  oldvol =        0.0;
  truvol = calc_vol();

  convergence = 100.;

  while(convergence > 1.){
    number_occupied = 0;

    for(j=0; j<rand_number; j++){
      F           = gsl_rng_uniform(gsl_ran_r);  // Chi limits satisfied by construction.

      rand_chi[j] = inverse_cumulative_nbar(F);

      cos_dec     = cos(rand_dec[j]);

      rand_x[j]   =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
      rand_y[j]   =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
      rand_z[j]   = -rand_chi[j]*sin(rand_dec[j]);

      xlabel      = (int) floor((rand_x[j] - min_x)/dx);
      ylabel      = (int) floor((rand_y[j] - min_y)/dy);
      zlabel      = (int) floor((rand_z[j] - min_z)/dz);

      boxlabel    = (int)    xlabel + m2*ylabel + m2*m1*zlabel;

      rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
    }

    for(j=0; j<m0*m1*m2; j++)  number_occupied   += rand_occupied[j];  // Number of cells over which <1+ delta> is averaged.

    newvol = dx*dy*dz*number_occupied/pow(10., 9.);

    convergence = 100.*(newvol-oldvol)/newvol;
    
    printf("\nExact volume: %.6lf (h^-1 Gpc)^3; randoms estimate: %.6lf (h^-1 Gpc)^3; %.6lf percent convergence; %.6lf percent accuracy", truvol, newvol, convergence, 100.*(truvol-newvol)/truvol);

    oldvol = newvol;
  }
  
  walltime("Walltime after randoms occupied");

  // %% calculation of clipping weights, must have no folding. %%
  double                       chi;
  double*             cell_weights;

  // Save cell weights in real part of H2_k to save memory/time.
  cell_weights  = malloc(n0*n1*n2*sizeof(*cell_weights));

  for(j=0; j<n0*n1*n2; j++){
     overdensity[j] = 0.0; // necessary? 
    cell_weights[j] = 1.0; // initialise to no clipping.
  }
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      chi         =  interp_comovingDistance(zobs[j]);

      xlabel      = (int)         floor((xCoor[j] - min_x)/dx);
      ylabel      = (int)         floor((yCoor[j] - min_y)/dy);
      zlabel      = (int)         floor((zCoor[j] - min_z)/dz);

      boxlabel    = (int)    xlabel + m2*ylabel + m2*m1*zlabel;

      overdensity[boxlabel]  +=  pow(dx*dy*dz*interp_nz(chi)*sampling[j], -1.); // Galaxy counts.
    }
  }

  // Smooth n/<n>.  True/false flag for zero mean; apodises boundaries.
  // Gaussian_filter(clipping_smoothing_radius, 0);

  // Scale (1 + delta) such that <1+ delta> = 1.; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0;
  double norm = 0.0;

  for(j=0; j<n0*n1*n2; j++){
    norm += overdensity[j];
  }
  
  norm /= number_occupied;

  printf("\n\nMean renormalisation in clipping weights: %.2lf", norm);
  
  
  return 0;
}


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
    // %% calculation of clipping weights, must have no folding. %% 
    double                       chi;
    double*             cell_weights;

    // Save cell weights in real part of H2_k to save memory/time. 
    cell_weights  = malloc(n0*n1*n2*sizeof(*cell_weights));
    
    for(j=0; j<n0*n1*n2; j++)  cell_weights[j] = 1.0; // initialise to no clipping.  

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){ 
        chi                     =  interp_comovingDistance(zobs[j]);
    
        boxlabel                =  boxCoordinates(xCoor, yCoor, zCoor, j);

        overdensity[boxlabel]  +=  pow(CellVolume*interp_nz(chi)*sampling[j], -1.); // Galaxy counts. 
      }
    }

    // Smooth n/<n>.  True/false flag for zero mean; apodises boundaries. 
    // Gaussian_filter(clipping_smoothing_radius, 0);
    
    // Scale (1 + delta) such that <1+ delta> = 1.; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0; 
    double norm = 0.0;

    for(j=0; j<n0*n1*n2; j++){
      // even if randoms don't fully sample the field then still a volume average over a representative sample.  
      if(rand_occupied[j] > 0.){
        norm += overdensity[j];
      }
    }

    norm /= number_occupied;

    printf("\n\nMean renormalisation in clipping weights: %.2lf", norm);

    for(j=0; j<n0*n1*n2; j++)  overdensity[j] -= 1.0;
    
    // By rescaling (1 + ddelta), don't slip below d=-1. 
    for(j=0; j<n0*n1*n2; j++){
      smooth_overdensity[j]  /=  norm;  
        
      smooth_overdensity[j]  -=   1.0; // now it's delta. 

      if(rand_occupied[j]    == 0.0){
               overdensity[j] = 0.0;
        smooth_overdensity[j] = 0.0;   // restore mask after smoothing. spurious empty cells where randoms should be?  
      }
    
      if(smooth_overdensity[j] > appliedClippingThreshold){  // either smooth_overdensity or overdensity. 	    	
        cell_weights[j]  = (1. + appliedClippingThreshold)/(1. + overdensity[j]);

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

    fclose(output); */
  }
  
  // walltime("Walltime after clipping weights");
  
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
