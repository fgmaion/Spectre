#define min(a, b) (((a) < (b)) ? (a) : (b))

int load_CovarianceMatrix(int mocks, int start){
    sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock", covariance_mocks_path, (int) round(appliedClippingThreshold), fieldFlag); 
  
    if(ChiSq_kmax <= jenkins_fold_kjoin){
        printf("\n\nLoading covariance without folding");
            
        load_CovarianceMatrix_withoutfolding(mocks, start, filepath);
    }
    
    else{
        printf("\n\nLoading covariance with folding");
        
        load_CovarianceMatrix_withfolding(mocks, start, filepath);
    }

    return 0;
}


int load_CovarianceMatrix_withoutfolding(int mocks, int start, char filepath[]){
    char     Nthfilepath[200];    
    char   firstfilepath[200];

    sprintf(firstfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, start, lo_zlim, hi_zlim);

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
    order      = mono_order*hiMultipoleOrder;
    
    assignCovMat(mocks); 
    
    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<mocks; k++){
        sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);	//  Be careful with 0 or 1 for initial mock. 
	
        inputfile = fopen(Nthfilepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++){  
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
            }
        }
        
        fclose(inputfile);  
    }

    snipping_amplitudeCorrection(Multipoles[k], mono_order);

    for(k=0; k<order; k++){
        for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;
    }
    
    for(k=0; k<order; k++){
      for(i=0; i<mocks; i++)  dMultipoles[i][k]  = Multipoles[i][k] - MeanMultipoles[k];  // Zero mean. 
    }
    
    // Covariance is an N x N matrix, where N corresponds to order, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
    
    initialiseCovariance(mocks);
    
    printf("\n\nMean multipoles:");
    
    printf("\nk \t \t mean mono \t mean quad \t sig_mono \t sig_quad \n");
    
    for(k=0; k<mono_order; k++){
        printf("\n%e \t %e \t %e \t %e \t %e", kVals[k], MeanMultipoles[k], MeanMultipoles[k + mono_order], sqrt(gsl_matrix_get(Covariance, k, k)), sqrt(gsl_matrix_get(Covariance, k+mono_order, k+mono_order)));
    }
    
    // Pre-whiten data and covariance. 
    prewhitenCov();
    
    // fprintf_Cov();
    
    // fprintf_meanMultipoes();
    
    Covariance_eigenVecs(mocks);
    
    return 0;
}


int load_CovarianceMatrix_withfolding(int mocks, int start, char filepath[]){
    char     Nthfilepath[200];    
    char   firstfilepath[200];
    char  foldedfilepath[200];
    char Nfoldedfilepath[200];

    sprintf(firstfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, start, lo_zlim, hi_zlim);
    
    printf("\n\n%s", firstfilepath);
    
    inputfile  = fopen(firstfilepath, "r");
    
    line_count(inputfile, &lineNo);

    for(i=0; i<lineNo; i++){
        fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);
        
        if(Interim<ChiSq_kmin){
          chiSq_kminIndex = i + 1;
        }
     
        if(Interim>jenkins_fold_kjoin){
            jenkins_foldIndex_unfoldedfile = i;
        
            break;
        }
    }
    
    fclose(inputfile);
    
    printf("\n\nk for switching to folded measurement: %.3lf (%d)", jenkins_fold_kjoin, jenkins_foldIndex_unfoldedfile);
    
    // again, but for folded measurement. 
    sprintf(foldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);
    
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
    
    printf("\nkmax limit for ChiSq: %.3lf (%d) \n\n", ChiSq_kmax, chiSq_kmaxIndex);
    
    mono_order = (jenkins_foldIndex_unfoldedfile - chiSq_kminIndex) + (chiSq_kmaxIndex - jenkins_foldIndex_foldedfile);
    order      = mono_order*hiMultipoleOrder;
    
    printf("\n\nmono order: %d, order: %d", mono_order, order);
    
    assignCovMat(mocks);
    
    for(k=0; k<mocks; k++){        
	sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);	

        inputfile = fopen(Nthfilepath, "r");

        for(i=0; i< jenkins_foldIndex_unfoldedfile; i++){  
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
            }
        }
        
        fclose(inputfile);
	 
	sprintf(Nfoldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, k + start, lo_zlim, hi_zlim);  // add in folded measurements, e.g. at k_join = 0.2;

        inputfile = fopen(Nfoldedfilepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++){  
            if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
            
            }
        }
        
        fclose(inputfile); 

	snipping_amplitudeCorrection(Multipoles[k], mono_order);
    }

    for(k=0; k<order; k++){
      for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;  // Be careful with 0 or 1 for the mock numbering.
      for(i=0; i<mocks; i++)  dMultipoles[i][k]  = Multipoles[i][k] - MeanMultipoles[k];  // New variables have zero mean.
    }
        
    // Covariance is an N x N matrix, where N corresponds to order, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
    
    initialiseCovariance(mocks);
    
    printf("\n\nMean multipoles:");
    
    printf("\nk \t \t mean mono \t mean quad \t sig_mono \t sig_quad \n");
    
    for(k=0; k<mono_order; k++){
        printf("\n%e \t %e \t %e \t %e \t %e", kVals[k], MeanMultipoles[k], MeanMultipoles[k + mono_order], sqrt(gsl_matrix_get(Covariance, k, k)), sqrt(gsl_matrix_get(Covariance, k+mono_order, k+mono_order)));
    }
    
    // fprintf_meanMultipoes();

    // Pre-whiten data and covariance. 
    prewhitenCov();
    
    // fprintf_Cov();
        
    Covariance_eigenVecs(mocks);
    
    return 0;
}


int initialiseCovariance(int mocks){    
    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
          Interim = 0.0;
        
          for(i=1; i<mocks; i++)  Interim += dMultipoles[i][j]*dMultipoles[i][k]/mocks;  
        
          gsl_matrix_set(Covariance, j, k, Interim);
        }
    }

    return 0;
}


int prewhitenCov(){
    for(j=0; j<order; j++){
        for(k=0; k<order; k++)  gsl_matrix_set(sigma_norm, j, k, 0.0);
        
        gsl_matrix_set(sigma_norm, j, j, pow(gsl_matrix_get(Covariance, j, j), -0.5));  // NB sigma_norm is (1./sigma)
    }
    
    for(j=0; j<order; j++){
      for(k=0; k<order; k++)  gsl_matrix_set(Covariance, j, k, gsl_matrix_get(Covariance, j, k)*gsl_matrix_get(sigma_norm, j, j)*gsl_matrix_get(sigma_norm, k, k));  // Pre-whitening. 
    }
    
    return 0;
}


int fprintf_meanMultipoes(){
    char meanmultipoles_filepath[200];

    if(ChiSq_kmax >= 0.799){
      sprintf(meanmultipoles_filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/mean_multipoles/d0_%d/W%d/meanmultipoles_zlim_%.1lf_%.1lf_old.dat", root_dir, (int) round(appliedClippingThreshold), fieldFlag, lo_zlim, hi_zlim);  
      
      output = fopen(meanmultipoles_filepath, "w"); 

      for(k=0; k<mono_order; k++)  fprintf(output, "%e \t %e \t %e \t %e \t %e \n", kVals[k], MeanMultipoles[k], sqrt(gsl_matrix_get(Covariance, k, k)), 
					                                                      MeanMultipoles[k + mono_order], sqrt(gsl_matrix_get(Covariance, k+mono_order, k+mono_order)));

      fclose(output);
    }    

    return 0;
}


int fprintf_Cov(){
    char cov_filepath[200];

    if(ChiSq_kmax >= 0.799){
      sprintf(cov_filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/cov_W%d_%.1lf_%.1lf_d0_%d", root_dir, fieldFlag, lo_zlim, hi_zlim, (int) round(appliedClippingThreshold));

      printf("\n\n%s", cov_filepath);

      output = fopen(cov_filepath, "w"); 
      
      for(j=0; j<order; j++){
        for(k=0; k<order; k++){        
          fprintf(output, "%e \t", gsl_matrix_get(Covariance, j, k));                  
        }
    
        fprintf(output, "\n");
      }

      fclose(output);
    }

    return 0;
}


int scale_Cov(int N){
  for(j=0; j<order; j++)  gsl_matrix_set(sigma_norm, j, j, sqrt((double) N)*gsl_matrix_get(sigma_norm, j, j));

  return 0;
}
