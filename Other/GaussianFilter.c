int prep_filterfactors(void){
  double xNy     = pi/dx; // Clipping weight grid has higher resolution than that for P(k) estimates.
  double yNy     = pi/dy; // CellSize > dx.
  double zNy     = pi/dz;

  double funda_x = 2.*pi*pow(n2, -1.)*pow(dx, -1.);
  double funda_y = 2.*pi*pow(n1, -1.)*pow(dy, -1.);
  double funda_z = 2.*pi*pow(n0, -1.)*pow(dz, -1.);

  walltime("Gaussian filter start");
  
  #pragma omp parallel for private(Index, k, j, i, k_z, k_y, k_x, kSq) if(thread == 1)
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      k_z = funda_z*k;
      k_y = funda_y*j;

      if(k_z > zNy)    k_z -= n0*funda_z;
      if(k_y > yNy)    k_y -= n1*funda_y;

      for(i=0; i<nx; i++){
        if(k_x > xNy)  k_x -= n2*funda_x;

        k_x   = funda_x*i;

        Index = k*n1*nx + j*nx + i;

        kSq   = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
        
        filter_factors[Index] = exp(-kSq*pow(smooth_radius, 2.)/2.)/(n0*n1*n2);

        // if(Index<20)  printf("\n%.6lf \t %.9lf", kSq, filter_factors[Index]);
      }
    }
  }

  walltime("Gaussian filter finish");

  return 0;
}


int Gaussian_filter(void){  
  // Gaussian filter the array overdensity, filter radius set by global variable GaussianFilter_radius.
  fftw_execute(plan); // plan works off n0, n1, n2.  
  
  #pragma omp parallel for private(k) if(thread == 1)
  for(k=0; k<n0*n1*nx; k++){
    H_k[k][0] *= filter_factors[k];
    H_k[k][1] *= filter_factors[k];
  }

  fftw_execute(iplan);

  return 0;
}
