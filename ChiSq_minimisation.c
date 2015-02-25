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
    
    fftlog_indices = malloc(mono_order*sizeof(*fftlog_indices));
    
    for(i=0; i<mono_order; i++){  
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
    
    printf("\n\nk vals match up.\n");
    
    for(j=0; j<mono_order; j++)  printf("%e \t %e \n", kVals[j], mono_config->krvals[fftlog_indices[j]][0]);
    
    printf("\nMatch up between observed k vals and FFTlog vals complete.\n");
    
    return 0;
}


int load_mock(int mockNumber){
    // load first mock again. Multipoles was used to store dMultipoles previously. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields/clipped_fullCube_500_kaiserLorentz_%d.dat", root_dir, mockNumber);
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, mockNumber);
    
    // sprintf(filepath,"%s/Data/500s/HOD_mocks_zobs_allgals_clipped/HOD_mock_512_%d.dat", root_dir, mockNumber);
    // sprintf(filepath,"%s/Data/500s/zobs_meanmultipoles.dat", root_dir);
            
    sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_specmask_%d.dat", root_dir, mockNumber);
            
    inputfile = fopen(filepath, "r");

    printf("\nCorrelated data, mock %d", mockNumber);

    for(i=0; i<chiSq_kmaxIndex; i++){   
            if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
            else{
                fscanf(inputfile, "%*le \t %le \t %le \t %*d \n", &xdata[i - chiSq_kminIndex], &xdata[chiSq_kmaxIndex + i - 2*chiSq_kminIndex]);
            
                printf("\n%e \t %e", xdata[i - chiSq_kminIndex], xdata[i - chiSq_kminIndex + mono_order]);
            }
               
            // mock avg. multipoles. 
            // fscanf(inputfile, "%*le \t %le \t %le \n", &xdata[i], &xdata[i + mono_order]);
            
            // printf("\n%e \t %e", xdata[i], xdata[i + mono_order]);
    }
    
    fclose(inputfile);    
    
    return 0;
}


int ChiSq_minimisation(){    
    prepChiSq_minimisation();
    
    // assigns memory for chi sq. grid, likelihood grid, and 1 and 2 parameter posteriors. 
    LikelihoodMemory();
    
    // pt2RSD_k         = &kaiserLorentz_multipole;
    
    // prep. clipping correction. 
    FFTlogRes        = 4096;
    
    // Must have mask multipoles available, uncomment prep_VIPERS_maskMultipoles().
    FFTlog_memory(FFTlogRes, 1., velDispersion);
    
    // connect measured k vals with their nearest FFTlog counterpart. 
    kvals_matchup();    
    
    load_mock(1);
    
    // new decorrelated, unit variance, zero mean variables. 
    for(j=0; j<order; j++){
        ydata[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
    }
    
    // input_check();
    
    printf("\n\n\nBeginning Chi sq. minimisation.");
  
    printf("\n\nPriors:");
    printf("\n%.2f < fsigma8    < %.2f",   min_fsigma8,     max_fsigma8);  
    printf("\n%.2f < vel. disp. < %.2f",   min_velDisperse, max_velDisperse);
    printf("\n%.2f < bsigma8    <  %.2f",  min_bsigma8,     max_bsigma8);
    printf("\n%.2f < A11 sq. < %.2f",      min_A11Sq,       max_A11Sq);
    printf("\n");
    
    double       sigmaInterval = (max_velDisperse - min_velDisperse)/dRes; 
    double       A11SqInterval = (max_A11Sq       - min_A11Sq)/dRes;
    
    double     fsigma8Interval = (max_fsigma8     - min_fsigma8)/dRes;
    double     bsigma8Interval = (max_bsigma8     - min_bsigma8)/dRes;
    
    
    // Mono and quad. fit
    double dof = order - paramNumber;
    
    printf("\n\nexpected chi sq. %.2f +- %.2f\n\n", dof, sqrt(2.*dof));
    
    
    sprintf(filepath, "%s/Data/500s/likelihoodGrid.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(jj=0; jj<Res; jj++){
        fsigma8 = min_fsigma8 + fsigma8Interval*jj;

        for(kk=0; kk<Res; kk++){
            bsigma8 = min_bsigma8 + bsigma8Interval*kk;

            for(ii=0; ii<Res; ii++){
                velDispersion = min_velDisperse + sigmaInterval*ii;
                
                // A11Sq = min_A11Sq + A11SqInterval*ii;
                A11Sq = 1.0;
                
                ChiSqGrid[jj][kk][ii]        =      ChiSqEval();
                    
                lnLikelihoodGrid[jj][kk][ii] = -0.5*ChiSqGrid[jj][kk][ii]; 
                
                fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %d \t %d \t %d \n", fsigma8, velDispersion, bsigma8, A11Sq, u0, ChiSqGrid[jj][kk][ii], jj, kk, ii);
                
                // printf("\n%e \t %e \t %e \t %e \t %e \t %e", fsigma8, bsigma8, velDispersion, A11Sq, u0, ChiSqGrid[jj][kk][ii]);
                
                if(ChiSqGrid[jj][kk][ii]<minChiSq){
                    minChiSq          = ChiSqGrid[jj][kk][ii];
                    
                    minChiSq_fsigma8  = fsigma8;
                    minChiSq_A11Sq    = A11Sq;
                    minChiSq_sigma    = velDispersion;
                    minChiSq_bsigma8  = bsigma8;
                    
                    printf("\n%e \t %e \t %e \t %e \t %e \t %e", fsigma8, velDispersion, bsigma8, A11Sq, u0, minChiSq);
                }
	       }
       }
    }
    
    fclose(output);
    /*
    fprintfBestfit();
    
    calc_fsigma8Posterior();
    */
    // calc_betaBiasPosterior();
    
    // calc_betaSigmaPosterior();
    
    return 0;
}


double ChiSqEval(){
    double ChiSq = 0.0;
    
    // printf("\n%e \t %e", fsigma8, bsigma8);
    
    FFTLog_setInput(mono_config, fsigma8/bsigma8, velDispersion);
    FFTLog_setInput(quad_config, fsigma8/bsigma8, velDispersion);
    
    // Transform to correlation functions, mono and quad. 
    xi_mu(mono_config);
    xi_mu(quad_config);

    // varCalc(mono_config, &variance, &u0);
    
    // u0 calc. in the absence of the value of \delta_0.
    // u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);
    // u0  = appliedClippingThreshold/sqrt(2.*variance);
    
    // printf("\nu0: %e", u0);
        
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config);
    
    // Currently evaluated at second order. 
    // clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    // clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // pk_mu(clipmono_config);
    // pk_mu(clipquad_config);
    
    // could include hex terms of correlation fn. 
    // cnvldmonoCorr(convlmonoCorr, clipmono_config, clipquad_config, zero_config);
    // cnvldquadCorr(convlquadCorr, clipmono_config, clipquad_config, zero_config);
    
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
    
    for(j=0; j<mono_order; j++)  xtheory[j]                                      =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    for(j=0; j<mono_order; j++)  xtheory[j + mono_order]                         =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
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
    
    for(j=0; j<order; j++)  ChiSq += pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j); 
    
    /*
    for(j=0; j<order; j++){  
        ChiSq += pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.); 
    }*/

    return ChiSq;
}


int fprintfBestfit(){
    fsigma8         = minChiSq_fsigma8;
    bsigma8         = minChiSq_bsigma8;
    velDispersion   = minChiSq_sigma;
    
    A11Sq           = 1.;
    // A11Sq        = minChiSq_A11Sq;

    printf("\n\nBest fit parameters, fsigma8: %.4e, vel. disp.: %.4e, bsigma8: %.4e, A11: %.4e", fsigma8, velDispersion, bsigma8, A11Sq);

    ChiSqEval();


    // Best fit correlated variables. 
    sprintf(filepath, "%s/Data/500s/clipped_bestfit_correlatedMonopole_kmax_%.2f.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_order; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j], 
                                                         xtheory[j], 
                                                         1./gsl_matrix_get(sigma_norm, j, j),  
                                                         pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.));
    
    fclose(output); 
    
    
    sprintf(filepath, "%s/Data/500s/clipped_bestfit_correlatedQuadrupole_kmax_%.2f.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_order; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j + mono_order], 
                                                         xtheory[j + mono_order], 
                                                         1./gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), 
                                                         pow(xdata[j+mono_order] - xtheory[j+mono_order], 2.)*pow(gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), 2.));
    
    fclose(output);
    
    return 0;
}


int fprintf_model(){
    fsigma8         = minChiSq_fsigma8;
    bsigma8         = minChiSq_bsigma8;
    velDispersion   = minChiSq_sigma;
    
    A11Sq           = 1.;
    // A11Sq        = minChiSq_A11Sq;
    
    ChiSqEval();

    // Best fit correlated variables. 
    sprintf(filepath, "%s/Data/500s/model_correlatedMonopole_kmax_%.2f.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<FFTlogRes; j++)  fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);

    fclose(output); 
    
    return 0;
}
