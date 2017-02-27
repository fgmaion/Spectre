double mock_mean_fsig8(int d0, int field, double kmax){
  // Given d0, field and k_max, calculate the mean fsig8 of the available mocks.
  double mock_fsig8;
  double mock_mean_fsig8 = 0.0;

  for(j=1; j<26; j++){
    // mocks 1 to 25.                                                                                                                                                                                                                        
    if((d0 == 4) && (field == 1))  sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new_mock_%d_0.6_0.9_snipCorrected.dat",     d0, field, kmax, j);
    else if(field == 4)            sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mock_%d_0.6_0.9_snipCorrected.dat",     d0, field, kmax, j);
    else                           sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new/new_mock_%d_0.6_0.9_snipCorrected.dat", d0, field, kmax, j);

    // printf("\n%s", filepath);
    
    inputfile = fopen(filepath, "r");

    fscanf(inputfile, "%le", &mock_fsig8);

    mock_mean_fsig8 += mock_fsig8;

    fclose(inputfile);
  }

  mock_mean_fsig8 /= 25.;

  return mock_mean_fsig8;
}


double calc_clipped_fsig8_cov(int d0i, int d0j, int field, double kmax){
  // calc d0_i, d0_j element of the covariance matrix. 
  const double fsig8_loz = 0.462;
  const double fsig8_hiz = 0.430;

  double result = 0.0;

  double mock_fsig8_d0i;
  double mock_fsig8_d0j;

  double mock_mean_fsig8_d0i = 0.0;
  double mock_mean_fsig8_d0j = 0.0;

  double cal_d0i, cal_d0j;

  mock_mean_fsig8_d0i = mock_mean_fsig8(d0i,  field, kmax);
  mock_mean_fsig8_d0j = mock_mean_fsig8(d0j,  field, kmax);

  cal_d0i = fsig8_loz/mock_mean_fsig8_d0i;
  cal_d0j = fsig8_loz/mock_mean_fsig8_d0j;

  // printf("\n\nCalibration factors: %.2lf \t %.2lf", cal_ii, cal_jj);

  for(j=1; j<26; j++){
    // mocks 1 to 26.                                                                                                                                                                                                                        
    if((d0i == 4) && (field == 1))    sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new_mock_%d_0.6_0.9_snipCorrected.dat",     d0i, field, kmax, j);
    else if(field == 4)               sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mock_%d_0.6_0.9_snipCorrected.dat",         d0i, field, kmax, j);
    else                              sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new/new_mock_%d_0.6_0.9_snipCorrected.dat", d0i, field, kmax, j);

    inputfile = fopen(filepath, "r");

    fscanf(inputfile, "%le", &mock_fsig8_d0i);

    fclose(inputfile);


    if((d0j == 4) && (field == 1))    sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new_mock_%d_0.6_0.9_snipCorrected.dat",     d0j, field, kmax, j);
    else if(field == 4)               sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mock_%d_0.6_0.9_snipCorrected.dat",         d0j, field, kmax, j);
    else                              sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/new/new_mock_%d_0.6_0.9_snipCorrected.dat", d0j, field, kmax, j);

    inputfile = fopen(filepath, "r");

    fscanf(inputfile, "%le", &mock_fsig8_d0j);

    // covariance element for calibrated fsig8. 
    result += (cal_d0i*mock_fsig8_d0i - fsig8_loz)*(cal_d0j*mock_fsig8_d0j - fsig8_loz);

    fclose(inputfile);

    // printf("\n%.2lf \t %.2lf \t %.6lf", mock_fsig8_ii, mock_fsig8_jj, result);
  }

  result /= 25.;

  // printf("\n%.6lf", result);

  return result;
}


int calc_bestfit_fsig8(int field, double kmax){
  // theory prediction. 
  const double fsig8_loz = 0.462;
  const double fsig8_hiz = 0.430;

  int thresholds[4] = {4, 6, 10, 1000};

  int    s, n = 4;

  double     D[n];   // Data for four thresholds.
  double   var[n];   // Sq. error on data.
  double cal_D[n];   // Calibration factors.

  double  cov[n][n]; // 4 x 4 covariance matrix. 
  double cov2[n][n]; // Copy of covariance, inversion seems to be in place. 

  double inva[n*n];

  gsl_matrix_view   m = gsl_matrix_view_array(cov2,  n, n);
  gsl_matrix_view inv = gsl_matrix_view_array(inva,  n, n);
  gsl_permutation*  p = gsl_permutation_alloc(n);

  // Covariance matrix calc. for calibrated fsig8.
  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      cov[ii][jj]  = calc_clipped_fsig8_cov(thresholds[ii], thresholds[jj], field, kmax);

      cov2[ii][jj] = cov[ii][jj];
    }

    var[ii] = cov[ii][ii];
  }
  /*
  // Assume diagonal covariance //  i.e. inverse variance weighting.                                                                                                                                                                        
  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      if(ii != jj) cov[ii][jj] = cov2[ii][jj] = 0.0;
    }
  }
  */                                                                                                                                                                                                                                         

  //** Test covariance. **//
  //cov[0][0] = var[0] = 0.50;
  //cov[0][1] = 0.25;  
  //cov[1][0] = 0.25;
  //cov[1][1] = var[1] = 0.75;
  
  printf("\n\nCovariance matrix of fsig8 with d0.\n");

  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      printf("%.6lf \t ", cov[ii][jj]);
    }

    printf("\n");
  }

  printf("\n\nCorrelation matrix of fsig8 with d0.\n");

  for(ii=0; ii<n; ii++){
    for(jj=0; jj<n; jj++){
      printf("%.6lf \t ", cov[ii][jj]/(sqrt(var[ii])*sqrt(var[jj])));
    }

    printf("\n");
  }
  
  gsl_linalg_LU_decomp(&m.matrix, p, &s);    
  gsl_linalg_LU_invert(&m.matrix, p, &inv.matrix);
  
  printf("\n\nInverse covariance: \n");
      
  for(ii=0; ii<n; ii++){
   for(jj=0; jj<n; jj++)  printf("%+.4lf \t", gsl_matrix_get(&inv.matrix, ii, jj));

   printf("\n");
  }
  
  //
  // Check inverse is right here. 
  //

  gsl_permutation_free(p);

  // F = numer/denom;
  double denom, numer;

  denom = 0.0;
  numer = 0.0;

  // Sum of all precision matrix elements.
  for(i=0; i<n; i++){
    for(j=0; j<n; j++)  denom += gsl_matrix_get(&inv.matrix, i, j);
  }

  // Calibration factors.                                                                                                                                                                                                                   
  for(ii=0; ii<n; ii++)  cal_D[ii] = fsig8_loz/mock_mean_fsig8(thresholds[ii],  field, kmax);

  // Load PDR-2 clipped f \sigma_8 data, and calibrate. 
  for(i=0; i<n; i++){
    sprintf(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3/data_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/data_0.6_0.9_snipCorrected.dat", thresholds[i], field, kmax);

    inputfile = fopen(filepath, "r");

    fscanf(inputfile, "%le", &D[i]);

    fclose(inputfile);

    // Calibrate data.                                                                                                                                                                                                                      
    D[i] *= cal_D[i];
  }
  
  printf("\n\nd0 \t D \t\t D_calib \t err_calib: \t frac err_calib: \n");

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
   
  return 0;
}


