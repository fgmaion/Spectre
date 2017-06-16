int Fisher_matrix(double dbsigma8, double dvelDispersion){
  double dfsigma8;

  if(lo_zlim == 0.6){
    dfsigma8  =   0.4618;
  }

  else if(lo_zlim == 0.9){
    dfsigma8  =   0.4300;
  }

  else{
    printf("\n\nChange z limits for forecast");
  }

  int         nn;
  int nparam = 3;
  int npoint = 4;
  
  double  dtheta; 
  
  double  params[3]  = {dfsigma8, dbsigma8, dvelDispersion};
  double sparams[3];
  
  // eqn. (15) of Blot 2016. 
  double shifts[4] = { 1.02,   0.98,    1.04,    0.96};
  double   amps[4] = {2./3., -2./3., -1./12., +1./12.};

  double dmodel_dparam[nparam][order];

  prep_ctype_ChiSq(); // matches kvals in each mock multipoles to FFTlog (rather than to all kVals).
  
  for(i=0; i<nparam; i++){    
    for(j=0; j<order; j++)  dmodel_dparam[i][j] = 0.0;
  }
  
  for(int jjj=0; jjj<nparam; jjj++){
    for(int iii=0; iii<npoint; iii++){
      for(j=0; j<nparam; j++)  sparams[j] = params[j];
    
      sparams[jjj]    *=        shifts[iii];
      dtheta           =   0.02*params[jjj];
      
      fsigma8          =         sparams[0];
      bsigma8          =         sparams[1];
      velDispersion    =         sparams[2];
         
      epsilon_pad      =                0.0;     
      alpha_pad        =                1.0; 
 
      model_compute(0, 0, 0, 0, 0);

      for(j=0; j<mono_order; j++)  dmodel_dparam[jjj][j]              += amps[iii]*convlmonoCorr->pk[fftlog_indices[j]][0]/dtheta;
      for(j=0; j<mono_order; j++)  dmodel_dparam[jjj][j + mono_order] += amps[iii]*convlquadCorr->pk[fftlog_indices[j]][0]/dtheta;
    }
  }

  int signum;
  
  gsl_matrix* inverse = gsl_matrix_alloc(order, order);
  gsl_permutation* p  =   gsl_permutation_alloc(order);

  gsl_matrix* FisherMat = gsl_matrix_alloc(nparam, nparam);
  gsl_matrix* FisherInv = gsl_matrix_alloc(nparam, nparam);
  gsl_permutation* q    =    gsl_permutation_alloc(nparam);
  
  gsl_linalg_LU_decomp(Covariance, p, &signum);  
  gsl_linalg_LU_invert(Covariance, p, inverse);
    
  for(j=0; j<nparam; j++){
    for(i=0; i<nparam; i++){
      Interim = 0.0;

      gsl_matrix_set(FisherMat, j, i, 0.0);
      
      for(int k=0; k<order; k++){
        for(int l=0; l<order; l++){
            Interim += dmodel_dparam[j][k]*dmodel_dparam[i][l]*gsl_matrix_get(inverse, k, l);
        }
      }
      
      gsl_matrix_set(FisherMat, j, i, Interim);
    }
  }

  // Fisher matrix
  sprintf(filepath, "%s/fisher/matrix_W%d_%.1lf_%.1lf_kmax_%.1lf.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);

  output = fopen(filepath, "w");

  for(j=0; j<nparam; j++){
    for(i=0; i<nparam; i++)  fprintf(output, "%.6lf \t", gsl_matrix_get(FisherMat, j, i));

    fprintf(output, "\n");
  }

  fclose(output);

  // Marinalised errors
  gsl_linalg_LU_decomp(FisherMat, q,   &signum);
  gsl_linalg_LU_invert(FisherMat, q, FisherInv);

  printf("\n\nMarginalised errors");
  
  for(j=0; j<nparam; j++){
    printf("\nFor kmax of %.1lf, %.8lf +- %.8lf", ChiSq_kmax, params[j], sqrt(gsl_matrix_get(FisherInv, j, j)));
  }

  sprintf(filepath, "%s/fisher/margerrs_W%d_%.1lf_%.1lf_kmax_%.1lf.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);

  output = fopen(filepath, "w");

  for(j=0; j<nparam; j++)  fprintf(output, "%.6lf \t %.6lf \n", params[j], sqrt(gsl_matrix_get(FisherInv, j, j)));

  fclose(output);
  
  printf("\n\n");
  
  return 0;
}

