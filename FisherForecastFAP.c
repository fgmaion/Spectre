int Fisher_matrix_fap(double dbsigma8, double dvelDispersion){
  double dFAP, dfsigma8;

  if(lo_zlim == 0.6){
    dFAP      =   0.9533;
    dfsigma8  =   0.4618;
  }

  else if(lo_zlim == 0.9){ 
    dFAP      = 1.456636;
    dfsigma8  =   0.4300; 
  }

  else{
    printf("\n\nChange z limits for forecast.");

    return 1;
  }
  
  int         nn;
  int nparam = 4;
  int npoint = 4;
  
  double  dtheta; 
  
  double  params[4]  = {dfsigma8, dbsigma8, dvelDispersion, dFAP};
  double sparams[4];
  
  // eqn. (15) of Blot 2016. 
  double shifts[4] = { 1.02,   0.98,    1.04,    0.96};
  double   amps[4] = {2./3., -2./3., -1./12., +1./12.};

  double     F, FAP;
    
  double dmodel_dparam[nparam][order];
  
  for(i=0; i<nparam; i++){    
    for(j=0; j<order; j++)  dmodel_dparam[i][j] = 0.0;
  }

  prep_ctype_ChiSq();
  
  // Calculate parameter derivatives with five point stencil. 
  for(int jjj=0; jjj<nparam; jjj++){
    for(int iii=0; iii<npoint; iii++){
      for(j=0; j<nparam; j++)  sparams[j] = params[j];
    
      sparams[jjj]    *=        shifts[iii];
      dtheta           =   0.02*params[jjj];
      
      fsigma8          =         sparams[0];
      bsigma8          =         sparams[1];
      velDispersion    =         sparams[2];
      FAP              =         sparams[3];

      F                =           dFAP/FAP; 
      
      epsilon_pad      = pow(F, 1./3.) - 1.;  // F    = (1. + epsilon)**3.     
      alpha_pad        =                1.0; 
 
      model_compute(0, 0, 0, 0, 0);

      for(j=0; j<mono_order; j++)  dmodel_dparam[jjj][j]              += amps[iii]*convlmonoCorr->pk[fftlog_indices[j]][0]/dtheta;
      for(j=0; j<mono_order; j++)  dmodel_dparam[jjj][j + mono_order] += amps[iii]*convlquadCorr->pk[fftlog_indices[j]][0]/dtheta;
    }

    printf("\n\n");

    for(j=0; j<mono_order; j++)  printf("\n%.6le \t %.6le", dmodel_dparam[jjj][j], dmodel_dparam[jjj][j + mono_order]);
  }
  
  int signum;

  // inverse covariance. 
  gsl_matrix* inverse   = gsl_matrix_alloc(order, order);
  gsl_permutation* p    =   gsl_permutation_alloc(order);

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
  sprintf(filepath, "%s/fisher/fap_matrix_W%d_%.1lf_%.1lf_kmax_%.1lf.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);
  
  output = fopen(filepath, "w");

  for(j=0; j<nparam; j++){
    for(i=0; i<nparam; i++)  fprintf(output, "%.6lf \t", gsl_matrix_get(FisherMat, j, i));

    fprintf(output, "\n");
  }

  fclose(output);


  // Margnialised errors. 
  gsl_linalg_LU_decomp(FisherMat, q,   &signum);
  gsl_linalg_LU_invert(FisherMat, q, FisherInv);

  printf("\n\nMarginalised errors.");
  
  for(j=0; j<nparam; j++){
    printf("\nFor kmax of %.1lf, %.8lf +- %.8lf", ChiSq_kmax, params[j], sqrt(gsl_matrix_get(FisherInv, j, j)));
  }
  
  sprintf(filepath, "%s/fisher/fap_margerrs_W%d_%.1lf_%.1lf_kmax_%.1lf.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);

  output = fopen(filepath, "w");

  for(j=0; j<nparam; j++)  fprintf(output, "%.6lf \t %.6lf \n", params[j], sqrt(gsl_matrix_get(FisherInv, j, j)));
  
  fclose(output);
    
  printf("\n\n");
  
  return 0;
}
