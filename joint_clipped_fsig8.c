double calc_fs8(double z){
  void*                  null;
  double  ln_a, f, fs8;

  ln_a  = -log(1. + z);
  
  f     = f_Om_545(ln_a, null);

  fs8   = f*camb_sig8; // camb_sig8 already for given redshift, not z=0;

  printf("\nHERE: %.3lf \t %.3lf \t %.3lf", z, f, camb_sig8);
  
  return fs8;
}

double calc_clipped_fsig8_cov(int d0i, int d0j, int field, double kmax, double z_eff, double fs8_zeff){
  // calc. d0_i, d0_j element of the covariance matrix. 
  double result = 0.0;

  double mock_fs8_d0i;
  double mock_fs8_d0j;

  double mean_fs8_d0i = 0.0;
  double mean_fs8_d0j = 0.0;

  double cal_d0i, cal_d0j;
  
  getmockmean_params(d0i);  mean_fs8_d0i = fsigma8;
  getmockmean_params(d0j);  mean_fs8_d0j = fsigma8;
  
  cal_d0i = fs8_zeff/mean_fs8_d0i;
  cal_d0j = fs8_zeff/mean_fs8_d0j;

  // printf("\n\nCalibration factors: %.2lf \t %.2lf", cal_d0i, cal_d0j);
  
  double fsig8s_d0i[CatalogNumber];
  double fsig8s_d0j[CatalogNumber];
  
  // Load mock fsig8s. 
  sprintf(filepath, "%s/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mocks_%.1lf_%.1lf.dat", outputdir, d0i, fieldFlag, ChiSq_kmax, lo_zlim, hi_zlim);
  
  inputfile = fopen(filepath, "r");

  for(j=0; j<CatalogNumber; j++)  fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*lf \t %*lf \t %*lf \n", &fsig8s_d0i[j]);
  // for(j=0; j<CatalogNumber; j++)  printf("\n%.4lf", fsig8s_d0i[j]);
  
  fclose(inputfile);
  
  // and for d0_j.
  sprintf(filepath, "%s/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mocks_%.1lf_%.1lf.dat", outputdir, d0j, fieldFlag, ChiSq_kmax, lo_zlim, hi_zlim);
  
  inputfile = fopen(filepath, "r");

  for(j=0; j<CatalogNumber; j++)  fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*lf \t %*lf \t %*lf \n", &fsig8s_d0j[j]);
  // for(j=0; j<CatalogNumber; j++)  printf("\n%.4lf", fsig8s_d0j[j]);
  
  fclose(inputfile);
  
  // Calculate matrix element. 
  for(j=0; j<CatalogNumber; j++)  result += (cal_d0i*fsig8s_d0i[j] - fs8_zeff)*(cal_d0j*fsig8s_d0j[j] - fs8_zeff);  // covariance element for calibrated fsig8.
  
  result /= (double) CatalogNumber;

  // printf("\n%.6lf", result);
  
  return result;
}


int calc_bestfit_fsig8(int field, double kmax, double z_eff){
  int  thresholds[4] = {4, 6, 10, 1000};

  int      s, n = 4; // Number of thresholds. 
  
  double       D[n]; // Data for four thresholds.
  double     var[n]; // Sq. error on data.
  double   cal_D[n]; // Calibration factors.

  double  cov[n][n]; // 4 x 4 covariance matrix. 
  double cov2[n][n]; // Copy of covariance, as inversion seems to be in place. 

  double  inva[n*n];

  gsl_matrix_view   m = gsl_matrix_view_array(&cov2,  n, n);
  gsl_matrix_view inv = gsl_matrix_view_array(inva,  n, n);
  gsl_permutation*  p = gsl_permutation_alloc(n);

  double fs8_zeff;

  
  fs8_zeff = calc_fs8(z_eff);

  printf("\n\nf * sigma_8 predictions: %.4lf", fs8_zeff);
  
  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
       cov[ii][jj] = calc_clipped_fsig8_cov(thresholds[ii], thresholds[jj], field, kmax, z_eff, fs8_zeff);  // Covariance matrix calc. for calibrated fsig8.

      cov2[ii][jj] = cov[ii][jj];
    }

    var[ii] = cov[ii][ii];
  }
  
  // Assume diagonal covariance //  i.e. inverse variance weighting.                                                                                                                           //                        
  // for(ii=0; ii<n; ii++){
  //   for(jj=0; jj<n; jj++){
  //     if(ii != jj) cov[ii][jj] = cov2[ii][jj] = 0.0;
  //   }
  // }                                                                                                                                                                                                                                        
  // Test covariance. //
  // cov[0][0] = var[0] = 0.50;
  // cov[0][1] = 0.25;  
  // cov[1][0] = 0.25;
  // cov[1][1] = var[1] = 0.75;
  
  printf("\n\nCovariance matrix of (calibrated) fsig8 with d0.\n");

  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      printf("%.6lf \t ", cov[ii][jj]);
    }

    printf("\n");
  }
  
  printf("\n\nCorrelation matrix of (calibrated) fsig8 with d0.\n");

  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      printf("%.6lf \t ", cov[ii][jj]/(sqrt(var[ii])*sqrt(var[jj])));
    }

    printf("\n");
  }
  
  gsl_linalg_LU_decomp(&m.matrix, p, &s);    
  gsl_linalg_LU_invert(&m.matrix, p, &inv.matrix);
  
  printf("\n\n(Calibrated) Precision matrix: \n");
      
  for(ii=0; ii<n; ii++){
   for(jj=0; jj<n; jj++)  printf("%+.4lf \t", gsl_matrix_get(&inv.matrix, ii, jj));

   printf("\n");
  }

  gsl_permutation_free(p);
  
  double numer, denom;  // F = numer/denom;

  numer = 0.0;
  denom = 0.0;

  for(i=0; i<n; i++){
    for(j=0; j<n; j++)  denom += gsl_matrix_get(&inv.matrix, i, j);  // Sum of all precision matrix elements.
  }

  printf("\n\nd0 \t mock mean fsig8 \t calibration \n");
  
  for(ii=0; ii<n; ii++){
    getmockmean_params(thresholds[ii]); // Sets fsigma8.

    cal_D[ii] = fs8_zeff/fsigma8;       // Calibration factors.

    printf("%d \t %.4lf \t %.4lf \n", thresholds[ii], fsigma8, cal_D[ii]);
  }
  
  // Load PDR-2 clipped f \sigma_8 data, and calibrate. 
  for(i=0; i<n; i++){
    sprintf(filepath, "%s/data_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/data_%.1lf_%.1lf.dat", outputdir, thresholds[i], field, kmax, lo_zlim, hi_zlim);

    inputfile = fopen(filepath, "r");

    fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*lf \t %*lf \t %*lf \n", &D[i]);

    fclose(inputfile);

    D[i] *= cal_D[i];  // Calibrate data.
  }

  
  printf("\n\nd0 \t D \t\t D_calib \t err_calib: \t calibrated frac err (%): \n");

  for(ii=0; ii<n; ii++)  printf("%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf% \n", thresholds[ii], D[ii]/cal_D[ii], D[ii], sqrt(cov[ii][ii]), 100.*sqrt(cov[ii][ii])/D[ii]);
  
  for(i=0; i<n; i++){
    for(j=0; j<n; j++)  numer     += gsl_matrix_get(&inv.matrix, i, j)*D[j];    
  }

  double F, var_F;
  
  F = numer/denom;

  // variance calc. 
  var_F = 0.0;

  for(i=0; i<n; i++){
    for(j=0; j<n; j++)  var_F += gsl_matrix_get(&inv.matrix, i, j);
  }

  var_F = 1.0/var_F;

  printf("\nFIELD: W%d \t \t KMAX: %.2lf \t \t RESULT: %.4lf \t \t ERROR:  %.4lf \t FRACTIONAL ERROR:  %.4lf%", field, kmax, F, sqrt(var_F), 100.*sqrt(var_F)/F);

  // Write to file. 
  sprintf(filepath, "%s/data_v1.7/fsig8/joint/W%d/kmax_%.1lf/data_%.1lf_%.1lf.dat", outputdir, field, kmax, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  fprintf(output, "%.6lf \t %.6lf", F, sqrt(var_F));
  
  return 0;
}


