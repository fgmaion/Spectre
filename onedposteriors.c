int set_minChiSq(){
  minChiSq = pow(10., 12.);

  // Data races currently:  #pragma omp parallel for private(k, j, i) reduction(min : minChiSq) if(thread == 1)
  for(k=0; k<Res; k++){
    for(j=0; j<Res; j++){
      for(i=0; i<Res; i++){        
        // printf("\nchi sq: %.6lf", ChiSqGrid[k][j][i][00][00]);
        // printf("\nthread id = %d", omp_get_thread_num());
        
        if(ChiSqGrid[k][j][i][0][0] < minChiSq){
          minChiSq    =            ChiSqGrid[k][j][i][0][0];

          minX2_fsig8 = min_fsigma8     + fsigma8Interval*k;
          minX2_bsig8 = min_bsigma8     + bsigma8Interval*j;
          minX2_sigp  = min_velDisperse +   sigmaInterval*i;
        }
      }
    }
  }

  // printf("\n\nMin. chi sq. = %.6lf", minChiSq);  

  return 0;
}


double calc_onedposteriors(double* maxL_fsig8, double* maxL_bsig8, double* maxL_sigv){
  double   bpost_norm = -99.; // norm of bsig8 posterior. 
  double   fpost_norm = -99.; // norm of fsig8 posterior.
  double   spost_norm = -99.; // norm of sigma posterior. 
  
  int     bpost_peak_ind = 0;
  int     fpost_peak_ind = 0;
  int     spost_peak_ind = 0;
  
  double   sigp_post[Res];
  double  bsig8_post[Res];
  double  fsig8_post[Res];
   
  for(k=0; k<Res; k++){
     sigp_post[k] = 0.0;
    fsig8_post[k] = 0.0;
    bsig8_post[k] = 0.0;
  }
  
  // int nthreads, tid;  
  #pragma omp parallel for private(k, j, i) if(thread == 1)
  for(k=0; k<Res; k++){
    for(j=0; j<Res; j++){
      for(i=0; i<Res; i++){
        ChiSqGrid[k][j][i][0][0] = exp(-ChiSqGrid[k][j][i][0][0]/2.);  // Turn into likelihood.
      }
    }
  }

  for(k=0; k<Res; k++){
    fsig8_post[k] = 0.0;
    bsig8_post[k] = 0.0;
     sigp_post[k] = 0.0;
  }

  #pragma omp parallel for reduction(+: fsig8_post[:Res], bsig8_post[:Res], sigp_post[:Res]) private(k, j, i) if(thread == 1)
  for(j=0; j<Res; j++){
    for(i=0; i<Res; i++){
      for(k=0; k<Res; k++){
        fsig8_post[k] += ChiSqGrid[k][j][i][0][0];
        bsig8_post[k] += ChiSqGrid[j][k][i][0][0];
         sigp_post[k] += ChiSqGrid[i][j][k][0][0];
      }
    }
  }

  // double sum = 0.0;
  // for(k=0; k<Res; k++)  sum += fsig8_post[k];
  // printf("\n\nSUM: %lf \t %lf", sum, 1.*Res*Res*Res);

  for(k=0; k<Res; k++){
    if(fsig8_post[k] > fpost_norm){
      fpost_norm     = fsig8_post[k];

      fpost_peak_ind = k;
    }

    if(bsig8_post[k] > bpost_norm){
      bpost_norm     = bsig8_post[k];

      bpost_peak_ind = k;
    }

    if(sigp_post[k] > spost_norm){
      spost_norm      = sigp_post[k];

      spost_peak_ind = k;
    }
  }

  *maxL_fsig8 = min_fsigma8      + fsigma8Interval*fpost_peak_ind;
  *maxL_bsig8 = min_bsigma8      + bsigma8Interval*bpost_peak_ind;
  *maxL_sigv  = min_velDisperse  +   sigmaInterval*spost_peak_ind;
  
  return 0.0;
}
