double safe_round(double min, double max, double space, double* base_min, double* base_max){
  *base_min = floor(min) - space;                       
  *base_max =  ceil(max) + space;                     
                       
  return 0;
}


int print_randoccupied(){
  int array[n0][n1];

  printf("\n\nPrinting rand. occupied.");

  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/overdensity/W%d/rand_z_%.1lf_%.1lf.dat", root_dir, fieldFlag, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      array[k][j] = 0;
      
      for(i=201; i<202; i++){  
        Index = k*n1*n2 + j*n2 + i;
        
        if(rand_occupied[Index] < 1)  array[k][j] = -1;
      }
     
      fprintf(output, "%d \t", array[k][j]);
    }

    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}


int set_randoccupied(){
  double      x2, z2;
  double c_ra, c_dec;

  c_ra    =  CentreRA*(pi/180.);
  c_dec   = CentreDec*(pi/180.);
    
  number_occupied = 0;
  
  // Surveyed volume need only be done once; no rotation required.
  rand_occupied = malloc(n2*n1*n0*sizeof(*rand_occupied));

  for(j=0; j<n0*n1*n2; j++)  rand_occupied[j] = 0;

  set_cnst_nbar(); // set pt2nz and recalculate spline for inverse nbar.  

  #pragma omp parallel for private(j, F, x2, z2) if(thread == 1)
  for(j=0; j<rand_number; j++){
    // thread-safe 
    drand48_r(&randBuffers[omp_get_thread_num()], &F);  // No rotation. Chi limits satisfied by construction. 
    
    rand_chi[j]    = inverse_cumulative_nbar(F);

    x2             =  rand_chi[j]*cos(rand_ra[j])*cos(rand_dec[j]);
    z2             = -rand_chi[j]*sin(rand_dec[j]);

    rand_x[j]      = -sin(c_dec)*x2 - cos(c_dec)*z2;
    rand_z[j]      =  cos(c_dec)*x2 - sin(c_dec)*z2;

    rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos(rand_dec[j]);

    rand_x[j]     +=  stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
    rand_y[j]     +=  stefano_trans_y;
    rand_z[j]     +=  stefano_trans_z;
  }

  // min_x = min_y = min_z =   0.0;  same volume - will need centering, and rotated. 
  // max_x = max_y = max_z = 800.0; 

  safe_round(arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number), 10., &min_x, &max_x);
  safe_round(arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number), 10., &min_y, &max_y);
  safe_round(arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number), 10., &min_z, &max_z);
  
  dx    = (max_x - min_x)/n2;  dy    = (max_y - min_y)/n1;  dz    = (max_z - min_z)/n0;
  
  printf("\n\nRandoms:");
  printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));

  printf("\n\nNew bounds:");
  printf("\nx: %.1lf \t %.1lf h^-1 Mpc, dx: %.6lf h^-1 Mpc", min_x, max_x, dx);
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc, dy: %.6lf h^-1 Mpc", min_y, max_y, dy);
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc, dz: %.6lf h^-1 Mpc", min_z, max_z, dz);

  // #pragma omp parallel for private(j, xlabel, ylabel, zlabel, boxlabel) if(thread == 1)
  for(j=0; j<rand_number; j++){
    xlabel      = (int) floor((rand_x[j] - min_x)/dx);
    ylabel      = (int) floor((rand_y[j] - min_y)/dy);
    zlabel      = (int) floor((rand_z[j] - min_z)/dz);

    boxlabel    = xlabel + n2*ylabel + n2*n1*zlabel;

    rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
  }

  // original estimate. 
  for(j=0; j<n0*n1*n2; j++)  number_occupied += rand_occupied[j];
  
  printf("\n\nSurveyed vol. calc.");

  double oldvol, newvol, truvol, convergence, accuracy;

  oldvol =        0.0;
  truvol = calc_vol();

  convergence  = 100.;

  while(convergence > 0.01){   
    // #pragma omp parallel for private(j, F, x2, z2, xlabel, ylabel, zlabel, boxlabel) if(thread == 1)
    for(j=0; j<rand_number; j++){
      drand48_r(&randBuffers[omp_get_thread_num()], &F);
    
      rand_chi[j]    = inverse_cumulative_nbar(F);
      
      x2             =  rand_chi[j]*cos(rand_ra[j])*cos(rand_dec[j]);
      z2             = -rand_chi[j]*sin(rand_dec[j]);

      rand_x[j]      = -sin(c_dec)*x2  - cos(c_dec)*z2;
      rand_z[j]      =  cos(c_dec)*x2  - sin(c_dec)*z2;

      rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos(rand_dec[j]);

      rand_x[j]     +=  stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
      rand_y[j]     +=  stefano_trans_y;
      rand_z[j]     +=  stefano_trans_z;
      
      xlabel         = (int) floor((rand_x[j] - min_x)/dx);
      ylabel         = (int) floor((rand_y[j] - min_y)/dy);
      zlabel         = (int) floor((rand_z[j] - min_z)/dz);

      boxlabel       = (int) xlabel + n2*ylabel + n2*n1*zlabel;
    
      if(rand_occupied[boxlabel] <  1){      
        rand_occupied[boxlabel]  = 1;  // binary array, is surveyed or not

        number_occupied         += 1;
      }
    }
  
    newvol      = dx*dy*dz*number_occupied/pow(10., 9.);
    
    convergence = 100.*(newvol-oldvol)/newvol;

    accuracy    = 100.*(truvol-newvol)/truvol;

    printf("\nExact volume: %.6lf (h^-1 Gpc)^3; randoms estimate: %.6lf (h^-1 Gpc)^3; %+.6lf percent convergence; %+.6lf percent accuracy", truvol, newvol, convergence, accuracy);
    
    oldvol = newvol;
  }
  
  // print_randoccupied();

  // roughly 20% of rand_occupied is empty.
  occupied_indices = malloc(number_occupied*sizeof(int));  // array with indices of cells that are occupied by randoms.

  Index = 0;

  for(j=0; j<n0*n1*n2; j++){
    if(rand_occupied[j] > 0){
      occupied_indices[Index] = j;

      Index += 1;
    }
  }

  printf("\n\nDifference between all cells and occupied is %.4lf", 100.*(n0*n1*n2 - number_occupied)/((double) n0*n1*n2));
  
  walltime("Walltime after randoms occupied");

  return 0;
}


