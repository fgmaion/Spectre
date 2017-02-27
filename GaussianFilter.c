int Gaussian_filter(double GaussianFilter_radius, int zero_mean){
  // Gaussian filter the array overdensity, filter radius set by global variable GaussianFilter_radius.
  prep_fftw();

  fftw_execute(p);

  // Gaussian smooth the counts.
  double GaussianFilter;

  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      for(i=0; i<n2; i++){
	k_x = kIntervalx*i;
	k_y = kIntervaly*j;
	k_z = kIntervalz*k;

	if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
	if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
	if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

	Index                              = k*n1*n2 + j*n2 + i;

	kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

	kmodulus                           = pow(kSq, 0.5);

	mu                                 = k_z/kmodulus;
	if(kmodulus < 0.000001)       mu   = 0.0;

	// Radius 1.
	GaussianFilter                     = exp(-1.*kSq*0.5*pow(GaussianFilter_radius, 2.));

	H_k[Index][0]                     *= GaussianFilter;
	H_k[Index][1]                     *= GaussianFilter;

	H_k[Index][0]                     /= n0*n1*n2;
	H_k[Index][1]                     /= n0*n1*n2;
      }
    }
  }

  if(zero_mean == 1){
    // Return a zero mean field.
    H_k[0][0] = 0.0;
    H_k[0][1] = 0.0;
  }

  int negkIndex;

  // Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
  for(k=n0-1; k>=n0/2; k--){
    for(j=0; j<n1; j++){
      for(i=0; i<n2; i++){
	negkIndex          = k*n1*n2 + j*n2 + i;

	Index              = 0;

	// zero maps to zero on reflection through the origin.
	if(i!=0)  Index   += (n2 - i);
	if(j!=0)  Index   += (n1 - j)*n2;
	Index   += (n0 - k)*n1*n2;

	H_k[negkIndex][0]  =     H_k[Index][0];
	H_k[negkIndex][1]  = -1.*H_k[Index][1];

	if(negkIndex == Index)   H_k[Index][1] = 0.0; // purely real
      }
    }
  }

  // in the k_z=0 plane one semi-circle is independent, k_y>0.
  for(j=n1-1; j>=n1/2; j--){
    for(i=0; i<n2; i++){
      negkIndex          = j*n2 + i;

      Index              = 0;

      // zero maps to zero on reflection through the origin.
      if(i!=0)  Index   += (n2 - i);
      Index   += (n1 - j)*n2;

      H_k[negkIndex][0]  =     H_k[Index][0];
      H_k[negkIndex][1]  = -1.*H_k[Index][1];

      if(negkIndex == Index)   H_k[Index][1] = 0.0;
    }
  }

  // on the line k_z=k_y=0, one half is independent, k_x>=0.
  for(i=n2-1; i>=n2/2; i--){
    negkIndex          = i;

    Index              = 0;

    // zero maps to zero on reflection through the origin.
    Index             += (n2 - i);

    H_k[negkIndex][0]  =      H_k[Index][0];
    H_k[negkIndex][1]  =  -1.*H_k[Index][1];

    if(negkIndex == Index)    H_k[Index][1] = 0.0;
  }
  
  iplan              = fftw_plan_dft_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(iplan);

  // should be unnecessary.
  // for(j=0; j<n0*n2*n1; j++)  smooth_overdensity[j][1] = 0.0;

  return 0;
}
