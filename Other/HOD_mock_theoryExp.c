int HOD_mock_theoryExp(){
    pt2Pk    = &splintHODpk;
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
    /*
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
    */
    
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, zero_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, zero_config);
    
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

    // for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[j]                      =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    // for(j=0; j<chiSq_kmaxIndex; j++)  xtheory[chiSq_kmaxIndex + j]    =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
    
    sprintf(filepath, "%s/Data/likelihood/HOD_mock_theoryExp.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<FFTlogRes; j++){ 
        if((mono_config->krvals[j][0] > 0.01) && (mono_config->krvals[j][0] < 1.)){ 
           fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0], quad_config->pk[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
        }
    }
    
    fclose(output);

    return 0;
}
