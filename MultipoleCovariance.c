int CovarianceMatrix(int mocks){
    // Retrieve the number of rows necessary before k value is greater than kmax for chi sq. evaluation. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar/clipped_fullCube_500_noCosVar_kaiserLorentz_%d.dat", root_dir, 0);
    sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, 0); // 500s mask
    
    inputfile  = fopen(filepath, "r");
    
    // Retrieve the number of lines in the input file. 
    ch         = 0;
    lineNo     = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch    == '\n')
       	  lineNo += 1;
    } while (ch  != EOF);

    rewind(inputfile);

    for(i=0; i<lineNo; i++){
        fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);
    
        if(Interim>ChiSq_kmax){
            chiSq_kmaxIndex = (i -1);
        
            break;
        }
    }
    
    fclose(inputfile);
    
    printf("\n\nkmax limit for ChiSq: %e (%d) \n\n", ChiSq_kmax, chiSq_kmaxIndex);
    
    // total number of data points.
    order = chiSq_kmaxIndex*hiMultipoleOrder;
    
    // Number of mocks, number of bins to kmax, number of multipoles to be fit. 
    assignCovMat(mocks);
    
    // Be careful with 0 or 1 for the mock numbering. 
    for(k=0; k<mocks; k++){
        // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar/clipped_fullCube_500_noCosVar_kaiserLorentz_%d.dat", root_dir, k);
        sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, k); // 500s mask
            
        inputfile = fopen(filepath, "r");

        for(i=0; i<chiSq_kmaxIndex; i++)  fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &kVals[i], &Multipoles[k][i], &Multipoles[k][i+chiSq_kmaxIndex]);
    
        fclose(inputfile);    
    }  
    
    for(k=0; k<order; k++){
        for(i=0; i<mocks; i++){
            MeanMultipoles[k] += Multipoles[i][k]/mocks;
        }
    }
    
    // printf("\n\nMean multipoles:");
    
    // for(k=0; k<chiSq_kmaxIndex; k++)  printf("\n%e \t %e \t %e", kVals[k], MeanMultipoles[k], MeanMultipoles[k+chiSq_kmaxIndex]);
    
    for(k=0; k<order; k++){
        // New variables have zero mean. 
        for(i=0; i<mocks; i++)  dMultipoles[i][k] = Multipoles[i][k] - MeanMultipoles[k];
    }
    
    // Covariance is an N x N matrix, where N corresponds to order, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
    
    initialiseCovariance(mocks);
    
    // Pre-whiten data and covariance. 
    prewhitenCov(mocks);
    
    fprintf_Cov();
    
    Covariance_eigenVecs(mocks);
    
    return 0;
}


int initialiseCovariance(int mocks){    
    // assumes dMultipoles have zero mean. 
    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
          Interim = 0.0;
        
          for(i=0; i<mocks; i++)  Interim += dMultipoles[i][j]*dMultipoles[i][k]/mocks;  
        
          gsl_matrix_set(Covariance, j, k, Interim);
        }
    }

    return 0;
}


int prewhitenCov(int mocks){
    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
          gsl_matrix_set(sigma_norm, j, k, 0.0);
        }
        
        gsl_matrix_set(sigma_norm, j, j, pow(gsl_matrix_get(Covariance, j, j), -0.5));
    }
    
    // Pre-whitening of covariance. 
    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
            gsl_matrix_set(Covariance, j, k, gsl_matrix_get(Covariance, j, k)*gsl_matrix_get(sigma_norm, j, j)*gsl_matrix_get(sigma_norm, k, k));
        }
    }
    
    return 0;
}


int fprintf_Cov(){
    sprintf(filepath, "%s/Data/likelihood/Clipped_mask500s_Covariance.dat", root_dir);

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
