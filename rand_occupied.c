double rd(double arg, double base){
  return base*floor(arg/base);
}


double ru(double arg, double base){
  return base*(1. + floor(arg/base));
}


int set_randoccupied(){
  number_occupied = 0;
  
  // Surveyed volume need only be done once; no rotation required.
  rand_occupied = malloc(n2*n1*n0*sizeof(*rand_occupied));

  for(j=0; j<n0*n1*n2; j++)  rand_occupied[j] = 0;

  set_cnst_nbar(); // set pt2nz and recalculate spline for inverse nbar.  

  for(j=0; j<rand_number; j++){
    // no rotation.
    F              = gsl_rng_uniform(gsl_ran_r);  // Chi limits satisfied by construction.

    rand_chi[j]    = inverse_cumulative_nbar(F);

    cos_dec        = cos(rand_dec[j]);

    rand_x[j]      =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
    rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
    rand_z[j]      = -rand_chi[j]*sin(rand_dec[j]);
  }

  // min_x = min_y = min_z =   0.0;  same volume - will need centering, and rotated. 
  // max_x = max_y = max_z = 800.0; 

  // what if negative, i.e. W4; limited padding.
  min_x = rd(arrayMin(rand_x, rand_number), 50.); min_y = rd(arrayMin(rand_y, rand_number), 50.); min_z = rd(arrayMin(rand_z, rand_number), 50.);
  max_x = ru(arrayMax(rand_x, rand_number), 50.); max_y = ru(arrayMax(rand_y, rand_number), 50.); max_z = ru(arrayMax(rand_z, rand_number), 50.);

  dx    = (max_x - min_x)/n2;  dy    = (max_y - min_y)/n1;  dz    = (max_z - min_z)/n0;

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

    boxlabel    = (int)    xlabel + n2*ylabel + n2*n1*zlabel;

    rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not
  }

  // original estimate. 
  for(j=0; j<n0*n1*n2; j++)  number_occupied += rand_occupied[j];
  
  printf("\n\nSurveyed vol. calc.");

  double oldvol, newvol, truvol, convergence, accuracy;

  oldvol =        0.0;
  truvol = calc_vol();

  convergence = 100.;

  while(convergence > 1.0){
    for(j=0; j<rand_number; j++){
      // no rotation.
      F              = gsl_rng_uniform(gsl_ran_r);  // Chi limits satisfied by construction.

      rand_chi[j]    = inverse_cumulative_nbar(F);

      cos_dec        = cos(rand_dec[j]);

      rand_x[j]      =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
      rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
      rand_z[j]      = -rand_chi[j]*sin(rand_dec[j]);
    
      xlabel         = (int) floor((rand_x[j] - min_x)/dx);
      ylabel         = (int) floor((rand_y[j] - min_y)/dy);
      zlabel         = (int) floor((rand_z[j] - min_z)/dz);

      boxlabel       = (int)    xlabel + n2*ylabel + n2*n1*zlabel;

      if(rand_occupied[boxlabel] == 0){      
        rand_occupied[boxlabel]  =  1;  // binary array, is surveyed or not

        number_occupied         +=  1;
      }
    }

    newvol      = dx*dy*dz*number_occupied/pow(10., 9.);

    convergence = 100.*(newvol-oldvol)/newvol;

    accuracy    = 100.*(truvol-newvol)/truvol;

    printf("\nExact volume: %.6lf (h^-1 Gpc)^3; randoms estimate: %.6lf (h^-1 Gpc)^3; %+.6lf percent convergence; %+.6lf percent accuracy", truvol, newvol, convergence, accuracy);

    oldvol = newvol;
  }

  // roughly 20% of rand_occupied is empty.
  occupied_indices = malloc(number_occupied*sizeof(int));  // array with indices of cells that are occupied by randoms.

  Index = 0;

  for(j=0; j<n0*n1*n2; j++){
    if(rand_occupied[j] == 1){
      occupied_indices[Index] = j;

      Index += 1;
    }
  }

  printf("\n\nDifference between all cells and occupied is %.4lf", 100.*(n0*n1*n2 - number_occupied)/((double) n0*n1*n2));

  walltime("Walltime after randoms occupied");

  return 0;
}
