int load_CovarianceMatrix(int mocks, int start){
    // Retrieve the number of rows necessary before k value is greater than kmax for chi sq. evaluation. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar/clipped_fullCube_500_noCosVar_kaiserLorentz_%d.dat", root_dir, 0);
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, 0); // 500s mask
    // sprintf(filepath,"%s/Data/500s/HOD_mocks_zobs_allgals_clipped/HOD_mock_512_%d.dat", root_dir, start);
    
    // sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_specmask_%d.dat", root_dir, start);
    
    // sprintf(filepath, "%s/W1_Spectro_V5_0/mocks_pk/mocks_W1_v8.0_500_001_Nagoya_v5_Samhain_multipoles.dat", root_dir, loopCount);
    
                                // ----------------------------------------------------------------------------//   
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_10/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_001_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir);
    
    // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk_MunichSuccess_clipped/clipped_d0_1000/mock_1_256_pk_d0_1000.00.dat", root_dir);

    sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/pk_smoothed_%.1lf_reflected_2field_nbar/mock", root_dir, fieldFlag, nz_smoothRadius);

    char   Nthfilepath[200];    
    char firstfilepath[200];

    sprintf(firstfilepath, "%s_%d_256_pk.dat", filepath, start);

    printf("\n\n%s", firstfilepath);

    inputfile  = fopen(firstfilepath, "r");
    
    // Retrieve the number of lines in the input file. 
    ch         = 0;
    lineNo     = 0;
    
    do{
        ch = fgetc(inputfile);        
    
        if(ch == '\n')
       	  lineNo += 1;
    } while(ch != EOF);

    rewind(inputfile);

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
    mono_order = chiSq_kmaxIndex - chiSq_kminIndex;
    
    order      = mono_order*hiMultipoleOrder;
    
    printf("\n\nmono order: %d", mono_order);
    
    // Number of mocks, number of bins to kmax, number of multipoles to be fit. 
    assignCovMat(mocks);
    
    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<mocks; k++){
        // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/new_pk_15_06_15/clipped_d0_1000/mock_%d_256_pk_d0_1000.00.dat", root_dir, k + start);
        
        sprintf(Nthfilepath, "%s_%d_256_pk.dat", filepath, k + start);
	
        inputfile = fopen(Nthfilepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++){  
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i - chiSq_kminIndex], &Multipoles[k][i - chiSq_kminIndex], &Multipoles[k][mono_order + i - chiSq_kminIndex]);
            
                // printf("\n%d \t %d \t %e \t %e \t %e", k, i, kVals[i - chiSq_kminIndex], Multipoles[k][i - chiSq_kminIndex], Multipoles[k][mono_order + i - chiSq_kminIndex]);
            }
        }
        
        fclose(inputfile);    
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
    
    // Pre-whiten data and covariance. 
    prewhitenCov();
    
    fprintf_Cov();
    
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
    }
    
    // Pre-whitening of covariance. 
    for(j=0; j<order; j++){
        for(k=0; k<order; k++)  gsl_matrix_set(Covariance, j, k, gsl_matrix_get(Covariance, j, k)*gsl_matrix_get(sigma_norm, j, j)*gsl_matrix_get(sigma_norm, k, k));
    }
    
    return 0;
}


int fprintf_meanMultipoes(){
    sprintf(filepath, "%s/Data/500s/spoc_zobs_meanmultipoles.dat", root_dir);

    output = fopen(filepath, "w"); 

    for(k=0; k<mono_order; k++)  fprintf(output, "%e \t %e \t %e \n", kVals[k], MeanMultipoles[k], MeanMultipoles[k + mono_order]);

    fclose(output);
    
    return 0;
}


int fprintf_Cov(){
    sprintf(filepath, "%s/W1_Spectro_V7_0/multipoles_cov_W%d.dat", root_dir, fieldFlag);

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
