double rd(double arg, double base){
  return base*floor(arg/base);
}


double ru(double arg, double base){
  return base*(1. + floor(arg/base));
}


int prep_randoccupied(){
  m0 =  n0; // or other. 
  m1 =  n1;
  m2 =  n2;
  
  // Surveyed volume need only be done once; no rotation required.  
  rand_occupied = malloc(m2*m1*m0*sizeof(*rand_occupied));

  for(j=0; j<m0*m1*m2; j++)  rand_occupied[j] = 0;

  set_cnst_nbar();
  
  rand_newchi_newbasis();
  
  // what if negative.  Limited padding.   
  min_x = rd(arrayMin(rand_x, rand_number), 50.); min_y = rd(arrayMin(rand_y, rand_number), 50.); min_z = rd(arrayMin(rand_z, rand_number), 50.);
  max_x = ru(arrayMax(rand_x, rand_number), 50.); max_y = ru(arrayMax(rand_y, rand_number), 50.); max_z = ru(arrayMax(rand_z, rand_number), 50.);

  dx    = (max_x - min_x)/m2;  dy    = (max_y - min_y)/m1;  dz    = (max_z - min_z)/m0;

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
  
  double oldvol, newvol, truvol, convergence, accuracy;

  oldvol =        0.0;
  truvol = calc_vol();

  convergence = 100.;

  while(convergence > 10.0){
    number_occupied = 0;

    rand_newchi_newbasis();
    
    for(j=0; j<rand_number; j++){      
      xlabel      = (int) floor((rand_x[j] - min_x)/dx);
      ylabel      = (int) floor((rand_y[j] - min_y)/dy);
      zlabel      = (int) floor((rand_z[j] - min_z)/dz);

      boxlabel    = (int)    xlabel + m2*ylabel + m2*m1*zlabel;

      rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
    }

    for(j=0; j<m0*m1*m2; j++)  number_occupied   += rand_occupied[j];  // Number of cells over which <1+ delta> is averaged.
    
    newvol      = dx*dy*dz*number_occupied/pow(10., 9.);

    convergence = 100.*(newvol-oldvol)/newvol;

    accuracy    = 100.*(truvol-newvol)/truvol;
    
    printf("\nExact volume: %.6lf (h^-1 Gpc)^3; randoms estimate: %.6lf (h^-1 Gpc)^3; %+.6lf percent convergence; %+.6lf percent accuracy", truvol, newvol, convergence, accuracy);

    oldvol = newvol;
  }
  
  // roughly 20% of rand_occupied is empty. 
  occupied_indices = malloc(number_occupied*sizeof(int));  // array with indices of cells that are occupied by randoms.
  
  Index = 0;

  for(j=0; j<m0*m1*m2; j++){
    if(rand_occupied[j] == 1){
      occupied_indices[Index] = j;
      
      Index += 1;
    }
  }

  printf("\n\nDifference between all cells and occupied is %.4lf", 100.*(m0*m1*m2 - number_occupied)/((double) m0*m1*m2));

  walltime("Walltime after randoms occupied");
  
  return 0;
}


int use_delta(){                                                                                                                                                                             
  for(j=0; j<n0*n1*n2; j++){                                                                                                                                                                      if(rand_occupied[j] == 1)  norm += overdensity[j];                                                                                                                                        }

  norm /= number_occupied;

  printf("\n\nMean renormalisation in clipping weights: %.4lf", norm);                                                                                                                                                                                                                                                                                                                      // If no smoothing: rescaling (1 + delta), subtract off to delta, reapply mask.
  for(j=0; j<n0*n1*n2; j++){                                                                                                                                                                      overdensity[j]       /= norm;                                                                                                                                                                overdensity[j]       -=  1.0;                                                                                                                                                                overdensity[j]       *=  rand_occupied[j];

     smooth_overdensity[j] =    overdensity[j];                                                                                                                                              
  }                                                                                                                                                                                           
  return 0;
}


int calc_clipping_weights(){
 if(d0 >= 1000){
   for(j=0; j<Vipers_Num; j++)  clip_galweight[j] = 1.0;

   return 0;
 }

 else if(foldfactor > 1.0){
   load_clippingweights();  // clip folded measurements using precomputed (foldfactor == 1) weights.

   return 0;
 }
  
 else{
  // %% calculation of clipping weights, must have no folding. %%
  double  chi;

  // Save cell weights in real part of H2_k to save memory/time.
  for(j=0; j<n0*n1*n2; j++) overdensity[j] = 0.0; // for each mock.
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + m2*ylabel + m2*m1*zlabel;

      overdensity[boxlabel]  +=  pow(dx*dy*dz*interp_nz(chi)*sampling[j], -1.); // N/<N>
    }
  }

  printf("\nx min:  %.3f \t x max:  %.3f", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
  printf("\ny min:  %.3f \t y max:  %.3f", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
  printf("\nz min:  %.3f \t z max:  %.3f", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));
  
  // Smooth N/<N>.  True/false flag for zero mean; apodises boundaries. Padding required?
  Gaussian_filter(smooth_radius, dx, dy, dz); // assumes m0 = n0 etc. 
  
  // Scale (1 + delta) such that <1+ delta> = 1. within surveyed region; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1.0;
  double      norm = 0.0;
  double frac_clip = 0.0;

  // By rescaling (1 + smooth delta), don't slip below delta = -1.
  for(j=0; j<n0*n1*n2; j++){
    smooth_overdensity[j]  *=  rand_occupied[j];

    norm                   += smooth_overdensity[j];
  }

  norm /= number_occupied;

  printf("\n\nMean renormalisation in clipping weights: %.4lf", norm);

  // H_k was used in Gaussian filter. reassign to be used to hold clipping weights for each cell. 
  for(j=0; j<n0*n1*n2; j++)      H_k[j][0] = 1.0; // This will hold the clipping weight for the cell (initialise to no clipping)  

  for(j=0; j<number_occupied; j++){
    i                       =  occupied_indices[j];

    smooth_overdensity[i]  /=  norm;

    smooth_overdensity[i]  -=   1.0; // now it's delta.

    if(smooth_overdensity[i] > d0){  // either smooth_overdensity or overdensity.
      H_k[i][0]             = (1. + d0)/(1. + overdensity[i]);

      frac_clip            += 1.0;    // must be surveyed, otherwise smooth_overdensity would be zero and cell would be bypassed.  
    }
  }

  frac_clip /= number_occupied;

  printf("\n\nd0: %lf, percentage of cells clipped: %lf", d0, 100.*frac_clip);
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j]  == true){
      // chi        =  interp_comovingDistance(zobs[j]);

      // cos_dec    =  cos(dec[j]*pi/180.0);

      // xCoor[j]   =  chi*cos( ra[j]*pi/180.0)*cos_dec;
      // yCoor[j]   =  chi*sin( ra[j]*pi/180.0)*cos_dec;
      // zCoor[j]   = -chi*sin(dec[j]*pi/180.0); // Includes reflection.

      xlabel     = (int)  floor((xCoor[j] - min_x)/dx);
      ylabel     = (int)  floor((yCoor[j] - min_y)/dy);
      zlabel     = (int)  floor((zCoor[j] - min_z)/dz);

      boxlabel   = (int)  xlabel + m2*ylabel + m2*m1*zlabel;
      
      clip_galweight[j]  =  H_k[boxlabel][0];

      if(clip_galweight[j] > 1.)  printf("\nspurious weight: %.2lf", clip_galweight[j]);
    }

    else{
      clip_galweight[j] = 0.0;
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
