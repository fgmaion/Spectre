int RandCoorCalc(){
  rand_rshift = malloc(NuRandoms*sizeof(*rand_rshift));
  rand_polar  = malloc(NuRandoms*sizeof(*rand_polar));
  rand_x      = malloc(NuRandoms*sizeof(*rand_x));
  rand_y      = malloc(NuRandoms*sizeof(*rand_y));
  rand_z      = malloc(NuRandoms*sizeof(*rand_z));

  printf("\nMemory for randoms successfully assigned.");

  for(j=0; j<NuRandoms; j++){                                                                                       
    // derived parameters 
    rand_polar[j] =  (pi/2.0) - (pi/180.0)*rand_dec[j];  // Converted to radians. 
    rand_rA[j]   *=  (pi/180.0);                         // Converted to radians.

    // Cosmology dependent, HOD mock parameters assumed, see header.h
    rand_x[j]      = rand_chi[j]*sin(rand_polar[j])*cos(rand_rA[j]);
    rand_y[j]      = rand_chi[j]*sin(rand_polar[j])*sin(rand_rA[j]);
    rand_z[j]      = rand_chi[j]*cos(rand_polar[j]);
    rand_rshift[j] = interp_inverseComovingDistance(rand_chi[j]);
    rand_rA[j]    /= (pi/180.0);                            // Converted to degrees
  }
  
  printf("\nRandoms coordinate calculation completed.");
  return 0;
}
