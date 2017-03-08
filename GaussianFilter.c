int Gaussian_filter(double radius, double dx, double dy, double dz){
  // Gaussian filter the array overdensity, filter radius set by global variable GaussianFilter_radius.

  fftw_execute(plan);
  
  // Gaussian smooth the counts.
  double factor;

  double xNy = pi/dx;                                     
  double yNy = pi/dy;                                       
  double zNy = pi/dz;                                       

  double funda_x = 2.*pi*pow(n2, -1.)*pow(dx, -1.);
  double funda_y = 2.*pi*pow(n1, -1.)*pow(dy, -1.);
  double funda_z = 2.*pi*pow(n0, -1.)*pow(dz, -1.);
  
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      for(i=0; i<n2; i++){
        k_x = funda_x*i;
        k_y = funda_y*j;
        k_z = funda_z*k;

        if(k_x>xNy)  k_x   -= n2*funda_x;
        if(k_y>yNy)  k_y   -= n1*funda_y;
        if(k_z>zNy)  k_z   -= n0*funda_z;
        
        Index               = k*n1*n2 + j*n2 + i;

        kSq                 = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
        
        factor              = exp(-kSq*pow(radius, 2.)/2.);

        H_k[Index][0]      *= factor;
        H_k[Index][1]      *= factor;
        
        H_k[Index][0]      /= n0*n1*n2;
        H_k[Index][1]      /= n0*n1*n2;
      }
    }
  }
  
  fftw_plan iplan;

  iplan = fftw_plan_dft_c2r_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_ESTIMATE);
  
  fftw_execute(iplan);

  return 0;
}
