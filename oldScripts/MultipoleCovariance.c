#define min(a, b) (((a) < (b)) ? (a) : (b))

int load_CovarianceMatrix(int mocks, int start){
    // Retrieve the number of rows necessary before k value is greater than kmax for chi sq. evaluation. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar/clipped_fullCube_500_noCosVar_kaiserLorentz_%d.dat", root_dir, 0);
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, 0); 
    // 500s mask
    // sprintf(filepath,"%s/Data/500s/HOD_mocks_zobs_allgals_clipped/HOD_mock_512_%d.dat", root_dir, start);    
    // sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_specmask_%d.dat", root_dir, start);
    // sprintf(filepath, "%s/W1_Spectro_V5_0/mocks_pk/mocks_W1_v8.0_500_001_Nagoya_v5_Samhain_multipoles.dat", root_dir, loopCount);
    
                                // ----------------------------------------------------------------------------//   
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_10/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    
    // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk_MunichSuccess_clipped/clipped_d0_1000/mock_1_256_pk_d0_1000.00.dat", root_dir);

    // sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/pk_smoothed_%.1lf_reflected_2field_nbar/mock", root_dir, fieldFlag, nz_smoothRadius);
    
    // Change in loadmock() in chisq_minimisation.c; pre-chapter 7 results. 
    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/old_d0_1000/d0_%d/W%d/mock", root_dir, (int) round(appliedClippingThreshold), fieldFlag);
    
    // Clipping results. 
    sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/pk/d0_%d/W%d/mock", root_dir, (int) round(appliedClippingThreshold), fieldFlag);

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


int get_lineNumber(){
  // Retrieve the number of lines in the input file.                                                                                                        
  ch         = 0;
  lineNo     = 0;

  do{
    ch = fgetc(inputfile);

    if(ch == '\n')
      lineNo += 1;
  } while(ch != EOF);

  rewind(inputfile);

  return 0;
}


int snipping_amplitudeCorrection(){
  // Note I:  As a difference estimator, covariance does not have to explicity corrected for shotnoise.
  // Note II: Assumes same rescaling for each field. 

  if(hi_zlim == 1.20){		 
    if(appliedClippingThreshold ==   10.){  
      Multipoles[k][i - chiSq_kminIndex]              *=   1.30;                                             
      Multipoles[k][mono_order + i - chiSq_kminIndex] *=   1.30;  
    }

   if(appliedClippingThreshold ==    6.){  
      Multipoles[k][i - chiSq_kminIndex]              *=   1.85;                   
      Multipoles[k][mono_order + i - chiSq_kminIndex] *=   1.85;
   }

   if(appliedClippingThreshold ==    4.){
      Multipoles[k][i - chiSq_kminIndex]              *=   3.00;                 
      Multipoles[k][mono_order + i - chiSq_kminIndex] *=   3.00;
   }
  }
		
  else if(hi_zlim == 1.0){
   if(appliedClippingThreshold ==    4.){
      Multipoles[k][i - chiSq_kminIndex]              *=   5.00;                                                                                
      Multipoles[k][mono_order + i - chiSq_kminIndex] *=   5.00;
   }
  }

  return 0;
}

int load_CovarianceMatrix_withoutfolding(int mocks, int start, char filepath[]){
    char     Nthfilepath[200];    
    char   firstfilepath[200];

    sprintf(firstfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, start, lo_zlim, hi_zlim);

    printf("\n\n%s", firstfilepath);
    
    inputfile  = fopen(firstfilepath, "r");
    
    get_lineNumber();

    for(i=0; i<lineNo; i++){
        fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);
        
        if(Interim<ChiSq_kmin){
            // printf("\n %d \t %e", i, Interim);
            chiSq_kminIndex = i + 1;
        }
     
        if(Interim>ChiSq_kmax){
            chiSq_kmaxIndex = i;
        
            break;
        }
    }
    
    fclose(inputfile);
    
    printf("\n\nkmin limit for ChiSq: %.3lf (%d)",    ChiSq_kmin, chiSq_kminIndex);
    printf("\nkmax limit for ChiSq: %.3lf (%d) \n\n", ChiSq_kmax, chiSq_kmaxIndex);
    
    // total number of data points.
    mono_order = (chiSq_kmaxIndex - chiSq_kminIndex);
    order      = mono_order*hiMultipoleOrder;
    
    printf("\n\nmono order: %d, order: %d \n\n", mono_order, order);
    
    // Number of mocks, number of bins to kmax, number of multipoles to be fit. 
    assignCovMat(mocks);
    
    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<mocks; k++){
        // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/new_pk_15_06_15/clipped_d0_1000/mock_%d_256_pk_d0_1000.00.dat", root_dir, k + start);
       
	sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);	
	
        inputfile = fopen(Nthfilepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++){  
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);

		snipping_amplitudeCorrection();
            }
        }
        
        fclose(inputfile);  
    }

    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<order; k++){
        for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;
    }
    
    for(k=0; k<order; k++){
        // New variables have zero mean. 
        for(i=0; i<mocks; i++)  dMultipoles[i][k]  = Multipoles[i][k] - MeanMultipoles[k];
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
    
    get_lineNumber();

    for(i=0; i<lineNo; i++){
        fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);
        
        if(Interim<ChiSq_kmin){
            // printf("\n %d \t %e", i, Interim);
            chiSq_kminIndex = i + 1;
        }
     
        if(Interim>jenkins_fold_kjoin){
            jenkins_foldIndex_unfoldedfile = i;
        
            break;
        }
    }
    
    fclose(inputfile);
    
    printf("\n\nkmin limit for ChiSq: %.3lf (%d)",    ChiSq_kmin, chiSq_kminIndex);
    
    printf("\n\nk for switching to folded measurement: %.3lf (%d)", jenkins_fold_kjoin, jenkins_foldIndex_unfoldedfile);
    
    // Similar again but for folded measurement. 
    sprintf(foldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);
    
    inputfile  = fopen(foldedfilepath, "r");
    
    get_lineNumber();

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
    
    // total number of data points.
    mono_order = (jenkins_foldIndex_unfoldedfile - chiSq_kminIndex) + (chiSq_kmaxIndex - jenkins_foldIndex_foldedfile);
    order      = mono_order*hiMultipoleOrder;
    
    printf("\n\nmono order: %d, order: %d", mono_order, order);
    
    // Number of mocks, number of bins to kmax, number of multipoles to be fit. 
    assignCovMat(mocks);
    
    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<mocks; k++){
        // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/new_pk_15_06_15/clipped_d0_1000/mock_%d_256_pk_d0_1000.00.dat", root_dir, k + start);
        
	sprintf(Nthfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, k + start, lo_zlim, hi_zlim);	

        inputfile = fopen(Nthfilepath, "r");

        for(i=0; i< jenkins_foldIndex_unfoldedfile; i++){  
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
       
                snipping_amplitudeCorrection();
            }
        }
        
        fclose(inputfile);  
        
        // printf("\n\n");
        
        // add in folded measurements from e.g. jenkins_fold_kjoin = 0.1; 	
	sprintf(Nfoldedfilepath, "%s_%d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, k + start, lo_zlim, hi_zlim);

        inputfile = fopen(Nfoldedfilepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++){  
            if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &Multipoles[k][mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
            
		snipping_amplitudeCorrection();
            }
        }
        
        fclose(inputfile); 
        
        // printf("\n\n\n");
    }

    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<order; k++){
        for(i=0; i<mocks; i++)  MeanMultipoles[k] += Multipoles[i][k]/mocks;
    }
    
    // fprintf_meanMultipoes();
    
    for(k=0; k<order; k++){
        // New variables have zero mean. 
        for(i=0; i<mocks; i++)  dMultipoles[i][k] = Multipoles[i][k] - MeanMultipoles[k];
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
    // assumes dMultipoles have zero mean. 
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
        
        // NB sigma_norm is   1./sigma 
        gsl_matrix_set(sigma_norm, j, j, pow(gsl_matrix_get(Covariance, j, j), -0.5));

	// printf("\n %e", gsl_matrix_get(sigma_norm, j, j));
    }
    
    // Pre-whitening of covariance. 
    for(j=0; j<order; j++){
        for(k=0; k<order; k++)  gsl_matrix_set(Covariance, j, k, gsl_matrix_get(Covariance, j, k)*gsl_matrix_get(sigma_norm, j, j)*gsl_matrix_get(sigma_norm, k, k));
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
