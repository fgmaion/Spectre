int prepChiSq_minimisation(){
    minChiSq = pow(10., 12.);

    // Observed data. 
    xdata    = malloc(order*sizeof(double));
    xtheory  = malloc(order*sizeof(double));
     
    // Decorrelated data. 
    ydata    = malloc(order*sizeof(double));
    ytheory  = malloc(order*sizeof(double));
    
    return 0;
}


int kvals_matchup(){
    // Given the k values of the measured P(k) multipoles, find closest fftlog modes. 
    double     diff;
    double min_diff;
    
    fftlog_indices = malloc(chiSq_kmaxIndex*sizeof(*fftlog_indices));
    
    printf("Match up between observed k vals and FFTlog vals complete.\n");
    
    for(i=0; i<chiSq_kmaxIndex; i++){  
        min_diff = pow(10., 99.);
    
        for(j=0; j<FFTlogRes; j++){  
            diff = fabs(mono_config->krvals[j][0] - kVals[i]);
            
            if(diff<min_diff){
                min_diff = diff;
            
                fftlog_indices[i]                 = j;
            }
        }
        
        if(100.*min_diff/kVals[i] > 1.)  printf("\nError matchup in k, %.2e percentage match in k of %.2e", 100.*min_diff/kVals[i], kVals[i]);
    }
    
    return 0;
}


int load_mock(int mockNumber){
    // load first mock again. Multipoles was used to store dMultipoles previously. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields/clipped_fullCube_500_kaiserLorentz_%d.dat", root_dir, mockNumber);
    sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, mockNumber);
            
    inputfile = fopen(filepath, "r");

    for(i=0; i<chiSq_kmaxIndex; i++)  fscanf(inputfile, "%*le \t %le \t %le \t %*d \n", &xdata[i], &xdata[i+chiSq_kmaxIndex]);
    
    fclose(inputfile);    
    
    return 0;
}


int input_check(){
    pt2Pk    = &linearPk_Gaussian;
    pt2RSD_k = &kaiserLorentz_multipole;
    
    // prep. clipping correction. 
    FFTlogRes        = 4096;
        
    // initialisation. 
    FFTlog_memory(FFTlogRes, beta, velDispersion);

    // Given current parameter combination, set P(k).
    FFTLog_setInput(mono_config, beta, velDispersion);
    FFTLog_setInput(quad_config, beta, velDispersion);
    
    // Transform to correlation functions, mono and quad. 
    xi_mu(mono_config);
    xi_mu(quad_config);

    varCalc(mono_config, &variance, &u0);
    
    // u0 calc. in the absence of the value of \delta_0.
    u0  = appliedClippingThreshold/sqrt(2.*variance);
    
    printf("\nu0, expected: %e, measured: %e", u0, inverse_erf(2.*sqrt(0.6125) - 1.));
        
    // Currently evaluated at second order. 
    clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // pk_mu(clipmono_config);
    // pk_mu(clipquad_config);
    
    // could include hex terms of correlation fn. 
    cnvldmonoCorr(convlmonoCorr, clipmono_config, clipquad_config, zero_config);
    cnvldquadCorr(convlquadCorr, clipmono_config, clipquad_config, zero_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    
    /*
    double cnvldpk_zero;

    // convolved pk zero point. 
    for(i=0; i<convlmonoCorr->N;   i++){ 
        if((convlmonoCorr->krvals[i][0]) >= pow(10., -2.)){
            cnvldpk_zero = convlmonoCorr->pk[i][0];
            
            printf("\n\nConvolved P(k) zero point: %e", cnvldpk_zero);
            
            break;
        }   
    }
    */
    
    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = convlmonoCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = convlquadCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[j]                      =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[chiSq_kmaxIndex + j]    =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
    
    sprintf(filepath, "%s/Data/likelihood/input_clipped_convld_Multipoles.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<FFTlogRes; j++){ 
        if((mono_config->krvals[j][0] > 0.01) && (mono_config->krvals[j][0] < 1.)){ 
           fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
            // fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0], quad_config->pk[j][0], clipmono_config->pk[j][0], clipquad_config->pk[j][0]);
        }
    }
    
    fclose(output);

    return 0;
}


int ChiSq_minimisation(){    
    
    prepChiSq_minimisation();
    
    // assigns memory for chi sq. grid, likelihood grid, and 1 and 2 parameter posteriors. 
    LikelihoodMemory();
    
    pt2Pk    = &linearPk_Gaussian;
    pt2RSD_k = &kaiserLorentz_multipole;
    
    // prep. clipping correction. 
    FFTlogRes        = 4096;
    
    // Must have mask multipoles available, uncomment prep_VIPERS_maskMultipoles().
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    // connect measured k vals with their nearest FFTlog counterpart. 
    kvals_matchup();    
    
    load_mock(3);
    
    // new decorrelated, unit variance, zero mean variables. 
    for(j=0; j<order; j++){
        ydata[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
    }
    
    // input_check();
    
    printf("\n\nBeginning Chi sq. minimisation.");
    
    printf("\n\nbeta: %.2f, sigma: %.2f, A11 sq.: %.2f, kmax: %.2f", beta, velDispersion, A11Sq, ChiSq_kmax);
  
    printf("\n\nPriors:");
    printf("\n%.2f < beta  < %.2f", min_beta, max_beta);  
    printf("\n%.2f < sigma < %.2f", min_velDisperse, max_velDisperse);
    printf("\n%.2f < A11 sq. < %.2f", min_A11Sq, max_A11Sq);
    printf("\n");
    
    double   betaInterval = (max_beta - min_beta)/dRes;
    double  sigmaInterval = (max_velDisperse - min_velDisperse)/dRes; 
    double  A11SqInterval = (max_A11Sq - min_A11Sq)/dRes;
    
    
    sprintf(filepath, "%s/Data/likelihood/clipped_mask500s_likelihoodGrid.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(jj=0; jj<Res; jj++){
        beta = min_beta + betaInterval*jj;

        for(kk=0; kk<Res; kk++){
            velDispersion = min_velDisperse + sigmaInterval*kk;

            for(ii=0; ii<Res; ii++){
                A11Sq = min_A11Sq + A11SqInterval*ii;
      
                ChiSqGrid[jj][kk][ii]        =      ChiSqEval();
                          
                lnLikelihoodGrid[jj][kk][ii] = -0.5*ChiSqGrid[jj][kk][ii]; 
                
                fprintf(output, "%e \t %e \t %e \t %e \t %e \t %d \t %d \t %d \n", beta, velDispersion, A11Sq, u0, ChiSqGrid[jj][kk][ii], jj, kk, ii);
                
                // printf("\n%e \t %e \t %e \t %e \t %e", beta, velDispersion, A11Sq, u0, ChiSqGrid[jj][kk][ii]);
                
                if(ChiSqGrid[jj][kk][ii]<minChiSq){
                    minChiSq       = ChiSqGrid[jj][kk][ii];
                    
                    minChiSq_beta  = beta;
                    minChiSq_A11Sq = A11Sq;
                    minChiSq_sigma = velDispersion;
                    
                    printf("\n%e \t %e \t %e \t %e \t %e", beta, velDispersion, A11Sq, u0, minChiSq);
                }
	       }
       }
    }
    
    fclose(output);
    
    // Mono and quad. fit
    double dof = 2.*chiSq_kmaxIndex - paramNumber;
    
    printf("\n\nMinimum chi sq.: %e, degrees of freedom %e +- %e", minChiSq, dof, sqrt(2.*dof));
    
    fprintfBestfit();
    
    // calc_betaPosterior();
    
    calc_betaSigmaPosterior();
    
    return 0;
}


double ChiSqEval(){
    double ChiSq = 0.0;
    
    // Given current parameter combination, set P(k).
    FFTLog_setInput(mono_config, beta, velDispersion);
    FFTLog_setInput(quad_config, beta, velDispersion);
    
    // Transform to correlation functions, mono and quad. 
    xi_mu(mono_config);
    xi_mu(quad_config);

    varCalc(mono_config, &variance, &u0);
    
    // u0 calc. in the absence of the value of \delta_0.
    u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);
    // u0  = appliedClippingThreshold/sqrt(2.*variance);
    
    // printf("\nu0: %e", u0);
        
    // Currently evaluated at second order. 
    clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // pk_mu(clipmono_config);
    // pk_mu(clipquad_config);
    
    // could include hex terms of correlation fn. 
    cnvldmonoCorr(convlmonoCorr, clipmono_config, clipquad_config, zero_config);
    cnvldquadCorr(convlquadCorr, clipmono_config, clipquad_config, zero_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    // pk_mu(convlhexCorr);
    
    /*
    double cnvldpk_zero;

    // convolved pk zero point. 
    for(i=0; i<convlmonoCorr->N;   i++){ 
        if((convlmonoCorr->krvals[i][0]) >= pow(10., -2.)){
            cnvldpk_zero = convlmonoCorr->pk[i][0];
            
            printf("\n\nConvolved P(k) zero point: %e", cnvldpk_zero);
            
            break;
        }   
    }
    */
    
    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = convlmonoCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = convlquadCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[j]                      =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[chiSq_kmaxIndex + j]    =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
    // new decorrelated, unit variance, zero mean variables. 
    for(j=0; j<order; j++){
        ytheory[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ytheory[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];
    }

    /*
    sprintf(filepath, "%s/Data/likelihood/Clippedbestfit.dat", root_dir);

    output = fopen(filepath, "w");

    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%d \t %e \t %e \t %e \t %e \n", j, 
                                                         ydata[j], 
                                                         ytheory[j], 
                                                         sqrt(gsl_vector_get(eval, j)), 
                                                         pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j));
        
    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%d \t %e \t %e \t %e \t %e \n", j+chiSq_kmaxIndex, 
                                                         ydata[j+chiSq_kmaxIndex], 
                                                         ytheory[j+chiSq_kmaxIndex], 
                                                         sqrt(gsl_vector_get(eval, j+chiSq_kmaxIndex)), 
                                                         pow(ydata[j+chiSq_kmaxIndex] - ytheory[j+chiSq_kmaxIndex], 2.)/gsl_vector_get(eval, j+chiSq_kmaxIndex));
    
    fclose(output);
    */
    
    for(j=0; j<order; j++){  
        if(gsl_vector_get(eval, j) > 1.3){
            ChiSq += pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j); 
        }
    }

    return ChiSq;
}


int fprintfBestfit(){
    beta          = minChiSq_beta;
    velDispersion = minChiSq_sigma;
    A11Sq         = minChiSq_A11Sq;

    // Given current parameter combination, set P(k).
    FFTLog_setInput(mono_config, beta, velDispersion);
    FFTLog_setInput(quad_config, beta, velDispersion);
    
    // Transform to correlation functions, mono and quad. 
    xi_mu(mono_config);
    xi_mu(quad_config);

    varCalc(mono_config, &variance, &u0);
    
    // u0 calc. in the absence of the value of \delta_0.
    u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);
        
    // Currently evaluated at second order. 
    clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // pk_mu(clipmono_config);
    // pk_mu(clipquad_config);
    
    // could include hex terms of correlation fn. 
    cnvldmonoCorr(convlmonoCorr, clipmono_config, clipquad_config, zero_config);
    cnvldquadCorr(convlquadCorr, clipmono_config, clipquad_config, zero_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    // pk_mu(convlhexCorr);
    
    /*
    double cnvldpk_zero;

    // convolved pk zero point. 
    for(i=0; i<convlmonoCorr->N;   i++){ 
        if((convlmonoCorr->krvals[i][0]) >= pow(10., -2.)){
            cnvldpk_zero = convlmonoCorr->pk[i][0];
            
            printf("\n\nConvolved P(k) zero point: %e", cnvldpk_zero);
            
            break;
        }   
    }
    */
    
    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = convlmonoCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = convlquadCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[j]                      =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[chiSq_kmaxIndex + j]    =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
    // new decorrelated, unit variance, zero mean variables. 
    for(j=0; j<order; j++){
        ytheory[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ytheory[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];
    }

    // Best fit correlated variables. 
    sprintf(filepath, "%s/Data/likelihood/Clipped_mask500s_bestfit_x.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j], 
                                                         xtheory[j], 
                                                         sqrt(gsl_matrix_get(Covariance, j, j)), 
                                                         pow(xdata[j] - xtheory[j], 2.)/gsl_matrix_get(Covariance, j, j));
        
    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j+chiSq_kmaxIndex], 
                                                         xtheory[j+chiSq_kmaxIndex], 
                                                         sqrt(gsl_matrix_get(Covariance, j+chiSq_kmaxIndex, j+chiSq_kmaxIndex)), 
                                                         pow(xdata[j+chiSq_kmaxIndex] - xtheory[j+chiSq_kmaxIndex], 2.)/gsl_matrix_get(Covariance, j+chiSq_kmaxIndex, j+chiSq_kmaxIndex));
    
    fclose(output);


    // Best fit decorrelated variables. 
    sprintf(filepath, "%s/Data/likelihood/Clipped_mask500s_bestfit_y.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%d \t %e \t %e \t %e \t %e \n", j, 
                                                         ydata[j], 
                                                         ytheory[j], 
                                                         sqrt(gsl_vector_get(eval, j)), 
                                                         pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j));
        
    for(j=0; j<chiSq_kmaxIndex; j++) 
        fprintf(output, "%d \t %e \t %e \t %e \t %e \n", j+chiSq_kmaxIndex, 
                                                         ydata[j+chiSq_kmaxIndex], 
                                                         ytheory[j+chiSq_kmaxIndex], 
                                                         sqrt(gsl_vector_get(eval, j+chiSq_kmaxIndex)), 
                                                         pow(ydata[j+chiSq_kmaxIndex] - ytheory[j+chiSq_kmaxIndex], 2.)/gsl_vector_get(eval, j+chiSq_kmaxIndex));
    
    fclose(output);
    
    return 0;
}
