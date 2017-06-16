#define min(a, b) (((a) < (b)) ? (a) : (b))

int initialiseCovariance(int mocks){
  for(j=0; j<order; j++){
    for(k=0; k<order; k++){
      Interim = 0.0;

      // Wrong for backwards compatability. i should begin at 0 and denominator should be: mocks - 1. 
      for(i=1; i<mocks; i++)  Interim += dMultipoles[i][j]*dMultipoles[i][k]/mocks;
      
      // Correct
      // for(i=0; i<mocks; i++)  Interim += dMultipoles[i][j]*dMultipoles[i][k]/(mocks - 1);
      
      gsl_matrix_set(Covariance, j, k, Interim);
    }
  }

  return 0;
}

int prewhitenCov(void){
  for(j=0; j<order; j++){
    for(k=0; k<order; k++)  gsl_matrix_set(sigma_norm, j, k, 0.0);

    gsl_matrix_set(sigma_norm, j, j, pow(gsl_matrix_get(Covariance, j, j), -0.5));  // N.B. sigma_norm is (1./sigma)
  }

  for(j=0; j<order; j++){
    for(k=0; k<order; k++)  gsl_matrix_set(Covariance, j, k, gsl_matrix_get(Covariance, j, k)*gsl_matrix_get(sigma_norm, j, j)*gsl_matrix_get(sigma_norm, k, k));  // pre-whitening.
  }

  return 0;
}

int load_CovarianceMatrix_withoutfolding(int mocks, int start, char filepath[]){
  char     Nthfilepath[200];
  char   firstfilepath[200];

  sprintf(firstfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, start, lo_zlim, hi_zlim);
  // sprintf(firstfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, start, lo_zlim, hi_zlim);
  
  printf("\n\n%s", firstfilepath);
  
  inputfile  = fopen(firstfilepath, "r");

  line_count(inputfile, &lineNo);
  
  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim<ChiSq_kmin)  chiSq_kminIndex = i + 1;

    if(Interim>ChiSq_kmax){
      chiSq_kmaxIndex = i;

      break;
    }
  }

  fclose(inputfile);
  
  mono_order = (chiSq_kmaxIndex - chiSq_kminIndex);  // Total number of data points.
  order      =  mono_order*hiMultipoleOrder;

  assignCovMat(mocks);
  
  // Be careful with 0 or 1 for the mock numbering.
  for(k=0; k<mocks; k++){
    sprintf(Nthfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, k + start, lo_zlim, hi_zlim);   //  Be careful with 0 or 1 for initial mock.
    // sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);
    
    inputfile = fopen(Nthfilepath, "r");

    for(i=0; i<chiSq_kmaxIndex; i++){
      if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

      else{
        fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
      }
    }

    fclose(inputfile);
  }
  
  for(kk=0; kk<mocks; kk++){
    for(ii=0; ii<mono_order; ii++)  Multipoles[kk][ii] -=  shotnoise_instances[kk]; // uses fold factor of 4.0 currently; subtract off monopole.
    // for(ii=0; ii<order;      ii++)  Multipoles[kk][ii] *=         mean_suppression;
  }
  
  for(k=0; k<order; k++){
    for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;
    for(i=0; i<mocks; i++)  dMultipoles[i][k]  = Multipoles[i][k] - MeanMultipoles[k];  // Zero mean.
  }

  // Covariance is an N x N matrix, where N corresponds to order, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are ignored.
  initialiseCovariance(mocks);

  printf("\n\nMean multipoles:");

  printf("\nk \t \t mean mono \t mean quad \t sig_mono \t sig_quad \n");

  for(k=0; k<mono_order; k++){
    printf("\n%e \t %e \t %e \t %e \t %e", kVals[k], MeanMultipoles[k], MeanMultipoles[k + mono_order], sqrt(gsl_matrix_get(Covariance, k, k)), sqrt(gsl_matrix_get(Covariance, k+mono_order, k+mono_order)));
  }
  
  return 0;
}


int get_kindices(int start, char filepath[]){
  char   firstfilepath[200];
  char  foldedfilepath[200];
  
  sprintf(firstfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, start, lo_zlim, hi_zlim);
  // sprintf(firstfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, start, lo_zlim, hi_zlim);
  
  printf("\n\n%s", firstfilepath);
  
  inputfile  = fopen(firstfilepath, "r");
    
  line_count(inputfile, &lineNo);

  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim<ChiSq_kmin){
      chiSq_kminIndex = i + 1;
    }

    if(Interim > jenkins_fold_kjoin){
      jenkins_foldIndex_unfoldedfile = i;

      break;
    }
  }

  fclose(inputfile);

  printf("\n\nk for switching to folded measurement: %.3lf (%d)", jenkins_fold_kjoin, jenkins_foldIndex_unfoldedfile);

  // again, but for folded measurement.
  sprintf(foldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);
  // sprintf(foldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);
  
  inputfile  = fopen(foldedfilepath, "r");

  line_count(inputfile, &lineNo);

  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim < jenkins_fold_kjoin){
      jenkins_foldIndex_foldedfile = i + 1;
    }

    if(Interim>ChiSq_kmax){
      chiSq_kmaxIndex = i;

      break;
    }
  }

  fclose(inputfile);

  mono_order = (jenkins_foldIndex_unfoldedfile - chiSq_kminIndex) + (chiSq_kmaxIndex - jenkins_foldIndex_foldedfile);
  order      = mono_order*hiMultipoleOrder;
  
  return 0;
}


int get_Multipoles(int mocks, int start, char filepath[]){
  char     Nthfilepath[200];
  char Nfoldedfilepath[200];

  for(k=0; k<mocks; k++){
    sprintf(Nthfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, k + start, lo_zlim, hi_zlim);
    // sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);
    
    inputfile = fopen(Nthfilepath, "r");

    for(i=0; i<jenkins_foldIndex_unfoldedfile; i++){
      if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

      else{
        fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
      }
    }

    fclose(inputfile);

    sprintf(Nfoldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, k + start, lo_zlim, hi_zlim);  // add in folded measurements, e.g. at k_join = 0.2;
    // sprintf(Nfoldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, k + start, lo_zlim, hi_zlim);
    
    inputfile = fopen(Nfoldedfilepath, "r");

    for(i=0; i<chiSq_kmaxIndex; i++){
      if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

      else{
        fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
      }
    }

    fclose(inputfile);
  }    

 return 0;
}


int get_kmaxes(int mocks, int start, char filepath[]){                                                                                                                                         int min_index = 0;                                                                                                                                                                                                                                                                                                                                                                     
  for(i=0; i<mono_order; i++){                                                                                                                                                               
    for(j=min_index; j<ChiSq_nkmaxes; j++){                                                                                                                                                  
      if(kVals[i] > ChiSq_kmaxes[j]){                                                                                                                                                        
        ChiSq_ikmaxes[j] = i;                                                                                                                                                                                                                                                                                                                                                            
        min_index = j + 1;

        break;                                                                                                                                                                              
      }                                                                                                                                                                                      
    }                                                                                                                                                                                        
  }
  
  return 0; 
}


int load_CovarianceMatrix_withfolding(int mocks, int start, char filepath[]){
  get_kindices(start, filepath);
  
  assignCovMat(mocks);
  
  get_Multipoles(mocks, start, filepath);
  
  // get_kmaxes(mocks, start, filepath);
  
  for(kk=0; kk<mocks; kk++){
    for(ii=0; ii<mono_order; ii++)  Multipoles[kk][ii] -= shotnoise_instances[kk];      // uses fold factor of 4.0 currently; subtract off monopole.
    // for(ii=0; ii<order;      ii++)  Multipoles[kk][ii] *=        mean_suppression;
  }
  
  for(k=0; k<order; k++){
    for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;                // Be careful with 0 or 1 for the mock numbering.
    for(i=0; i<mocks; i++)  dMultipoles[i][k]  = Multipoles[i][k] - MeanMultipoles[k];  // New variables have zero mean.
  }
  
  // Covariance is an N x N matrix, where N corresponds to order, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements.  Hex-blah elements are ignored. 
  initialiseCovariance(mocks);
    
  printf("\n\nMean multipoles:");
    
  printf("\nk \t \t mean mono \t mean quad \t sig_mono \t sig_quad \n");
    
  for(k=0; k<mono_order; k++){
    printf("\n%e \t %e \t %e \t %e \t %e", kVals[k], MeanMultipoles[k], MeanMultipoles[k + mono_order], sqrt(gsl_matrix_get(Covariance, k, k)),
                                                                                                        sqrt(gsl_matrix_get(Covariance, k + mono_order, k + mono_order)));
  }
  
  return 0;
}


int load_CovarianceMatrix(int mocks, int start){
  assign_chisq_kmaxes();
  
  sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock", covariance_mocks_path, d0, fieldFlag);

  if(ChiSq_kmax <= jenkins_fold_kjoin){
    printf("\n\nLoading covariance without folding: %s", filepath);

    load_CovarianceMatrix_withoutfolding(mocks, start, filepath);
  }
  
  else{
    printf("\n\nLoading covariance with folding: %s", filepath);

    load_CovarianceMatrix_withfolding(mocks, start, filepath);
  }

  printCov();

  print_meanMultipoes();
  
  return 0;
}

int printCov(){
  sprintf(filepath, "%s/mocks_v1.7/pk_derivedprops/d0_%d/W%d/cov_zlim_%.1lf_%.1lf_kmax_%.1lf.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);

  output = fopen(filepath, "w");

  for(j=0; j<order; j++){
    for(k=0; k<order; k++){
      fprintf(output, "%e \t", gsl_matrix_get(Covariance, j, k));
    }
    
    fprintf(output, "\n");
  }
  
  fclose(output);

  return 0;
}
