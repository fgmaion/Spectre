int prepChiSq_minimisation(){
    minChiSq =                pow(10., 12.);

    // Observed data. 
    xdata    = malloc(order*sizeof(double));
    xtheory  = malloc(order*sizeof(double));
     
    // Decorrelated data. 
    ydata    = malloc(order*sizeof(double));
    ytheory  = malloc(order*sizeof(double));
    
    return 0;
}


int kvals_matchup(){
    // Given the k values of the measured P(k) multipoles, find nearest fftlog modes. 
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
       
        // if(100.*min_diff/kVals[i] > 1.)  printf("\nError matchup in k, %.2e percentage match in k of %.2e", 100.*min_diff/kVals[i], kVals[i]);
    }
    
    printf("\n\nk vals match up.\n");
    
    printf("\nMatch up between observed k vals and FFTlog vals complete.\n");
    
    for(j=0; j<mono_order; j++)  printf("%e \t %e \n", kVals[j], mono_config->krvals[fftlog_indices[j]][0]);
    
    return 0;
}


int load_aptest_multipoles(){
  sprintf(filepath, "%s/W1_Spectro_V7_2/AP_Multipoles_Ballinger_pk_prediction_1.03_0.03_0.50_1.00_6.00.dat", root_dir);

  inputfile = fopen(filepath, "r");
  
  printf("\n\nAP test multipoles");

  for(i=0; i<chiSq_kmaxIndex; i++){
    if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%*le \t %le \t %le \t %*d \n", &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);

      printf("\n%le \t %le", xdata[i - chiSq_kminIndex], xdata[mono_order + i - chiSq_kminIndex]);
    }
  }
  
  fclose(inputfile);
  
  /*
  printf("\n\nWith hypothetical errors.");

  // Add hypothetical errors. Diagonal covariance, NB sigma_norm is 1./sigma
  for(j=0; j<mono_order; j++){  
    xdata[j]              += rand_gaussian(gsl_ran_r, 1./gsl_matrix_get(sigma_norm, j, j));

    xdata[j + mono_order] += rand_gaussian(gsl_ran_r, 1./gsl_matrix_get(sigma_norm, j + mono_order, j + mono_order));

    // printf("\n%le \t %le \t %le \t %le \t %le \t %le \t %le", kVals[j], xdata[j], 1./gsl_matrix_get(sigma_norm, j, j), xdata[j + mono_order], 1./gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), 
    //                                                                                                  AP_P0(kVals[j], 0.5, 6.0, 0.03228, 1.03228), AP_P2(kVals[j], 0.5, 6.0, 0.03228, 1.03228));
  }

  printf("\n\n");
  */
  /*
  sprintf(filepath, "%s/W1_Spectro_V7_2/aptest_multipoles.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<mono_order; j++){
    fprintf(output, "\n%le \t %le \t %le \t %le \t %le \t %le \t %le", kVals[j], xdata[j], 1./gsl_matrix_get(sigma_norm, j, j), xdata[j + mono_order], 1./gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order),
	                                                                         AP_P0(kVals[j], 0.5, 6.0, 0.03228, 1.03228), AP_P2(kVals[j], 0.5, 6.0, 0.03228, 1.03228));
														      
  }

  fclose(output);
  */
  return 0;
}


int set_meanmultipoles(){
  for(i=0; i<order; i++){
    xdata[i] = MeanMultipoles[i];

    // printf("xdata: %.3le\n", xdata[i]);
  }

  return 0;
}


int load_meanmultipoles(){
  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/mean_multipoles/d0_%d/W%d/meanmultipoles_zlim_%.1lf_%.1lf.dat", root_dir, (int) round(appliedClippingThreshold), fieldFlag, lo_zlim, hi_zlim);

  inputfile = fopen(filepath, "r");

  for(i=0; i<chiSq_kmaxIndex; i++){
    if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %le \t %*le \t %le \t %*le \n");

    else{
      fscanf(inputfile, "%*le \t %le \t %*le \t %le \t %*le \n", &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);

      printf("\n%le \t %le", xdata[i - chiSq_kminIndex], xdata[mono_order + i - chiSq_kminIndex]);
    }
  }

  fclose(inputfile);

  return 0;
}


int load_mock(int mockNumber){
    // load first mock again. Multipoles was used to store dMultipoles previously. 
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields/clipped_fullCube_500_kaiserLorentz_%d.dat", root_dir, mockNumber);
    // sprintf(filepath,"%s/Data/likelihood/ClippedGaussian_fields_noCosVar_window500s_zeromean/clipped_fullCube_500_noCosVar_window500s_zeromean_kaiserLorentz_%d.dat", root_dir, mockNumber);
    
    // sprintf(filepath,"%s/Data/500s/HOD_mocks_zobs_allgals_clipped/HOD_mock_512_%d.dat", root_dir, mockNumber);
    // sprintf(filepath,"%s/Data/500s/zobs_meanmultipoles.dat", root_dir);
            
    // sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_specmask_%d.dat", root_dir, mockNumber);
                                                                  // ----------------------------------------------------------------------------//   
    // if(mockNumber<10)        sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/mocks_W1_v8.0_500_00%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir, mockNumber);
    // else if(mockNumber<100)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/mocks_W1_v8.0_500_0%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",  root_dir, mockNumber);
    // else                     sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",   root_dir, mockNumber);
       
    // if(mockNumber<10)        sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_10/mocks_W1_v8.0_500_00%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir, mockNumber);
    // else if(mockNumber<100)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_10/mocks_W1_v8.0_500_0%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",  root_dir, mockNumber);
    // else                     sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_10/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",   root_dir, mockNumber);
       
    // if(mockNumber<10)        sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_00%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir, mockNumber);
    // else if(mockNumber<100)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_0%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",  root_dir, mockNumber);
    // else                     sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",   root_dir, mockNumber);
       
    // new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight.
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/new_pk_15_06_15/clipped_d0_1000/mock_%d_256_pk_d0_1000.00.dat", root_dir, mockNumber);
    
    /// sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/old_d0_1000/d0_%d/W%d/mock_%d", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber);

    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/true_nbar/d0_%d/W%d/mock_%d", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber);
   
    //  Chapter 7, snipping results.
    sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/pk/d0_%d/W%d/mock_%d", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber);

    if(ChiSq_kmax <= jenkins_fold_kjoin){
        printf("\n\nLoading mock without folding");

	// amplitude rescaling in load_withoutfolding.
        load_withoutfolding(filepath);
    }

    else{
        printf("\n\nLoading mock with folding");
        
	// amplitude recsaling in load_withfolding. 
        load_withfolding(filepath);
    }

    // Snipped results are raw power, need shot noise corrected. 
    // d0 = {1000, 10, 6}, Ps = {277.16, 243.944925, 208.99925}.  Determined from Jf_8 results. 

    if(fieldFlag == 1){
      // W1

      if(hi_zlim = 1.20){
	if(appliedClippingThreshold      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 277.16;  // Correct monopole. 
	else if(appliedClippingThreshold ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 243.94;  // Note: as a difference estimator, covariance does not have to explicity corrected. 
	else if(appliedClippingThreshold ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 208.99;
	else if(appliedClippingThreshold ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 167.90;
      }

      else if(hi_zlim = 1.00){
        if(appliedClippingThreshold      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 700.16;  // 700. seems spurious. aliasing?
	else if(appliedClippingThreshold ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 230.30;
      }
    }


    if(fieldFlag == 4){
      // W4
      if(appliedClippingThreshold      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 289.16;  // Correct monopole.                                                                                       
      else if(appliedClippingThreshold ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 249.26;  // Note: as a difference estimator, covariance does not have to explicity corrected.                       
      else if(appliedClippingThreshold ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 208.00;
      else if(appliedClippingThreshold ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 168.59;
    }

    return 0;
}


int load_withoutfolding(char filepath[]){
    char  firstfilepath[200];

    printf("\n\nCorrelated data.");

    sprintf( firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, lo_zlim, hi_zlim);

    inputfile = fopen(firstfilepath, "r");
    
    for(i=0; i<chiSq_kmaxIndex; i++){  
        if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
        else{
            fscanf(inputfile, "%*le \t %le \t %le \t %*d \n", &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);
        
            printf("\n%le \t %le", xdata[i - chiSq_kminIndex], xdata[mono_order + i - chiSq_kminIndex]);
        }
    }
        
    fclose(inputfile);  
    
    return 0;
}


int load_withfolding(char filepath[]){
    double loc_k;

    char  firstfilepath[200];
    char foldedfilepath[200];
    
    printf("\n\nCorrelated data.");

    sprintf( firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, lo_zlim, hi_zlim);
    sprintf(foldedfilepath, "%s_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, lo_zlim, hi_zlim);
    
    inputfile = fopen(firstfilepath, "r");
    
    for(i=0; i<jenkins_foldIndex_unfoldedfile; i++){  
        if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
        else{
  	    fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &loc_k, &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);
        
            printf("\n%le \t %le \t %le", loc_k, xdata[i - chiSq_kminIndex], xdata[mono_order + i - chiSq_kminIndex]);
        }
    }
        
    fclose(inputfile);  
    
    // Add in folded measurements from e.g. jenkins_fold_kjoin = 0.1; 
    inputfile = fopen(foldedfilepath, "r");

    for(i=0; i<chiSq_kmaxIndex; i++){  
        if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");
        
        else{
	  fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &loc_k, &xdata[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &xdata[mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
            
	  printf("\n%le \t %le \t %le", loc_k, xdata[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], xdata[mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
        }
    }
        
    fclose(inputfile); 
    
    return 0;
}


int load_data(){            
    // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/d0_%d/W%d/data", root_dir, (int) round(appliedClippingThreshold), fieldFlag);
    // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/old_d0_%d/W%d/data", root_dir, (int) round(appliedClippingThreshold), fieldFlag);   

    sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/pk/d0_%d/W%d/data", root_dir, (int) round(appliedClippingThreshold), fieldFlag);

    if(ChiSq_kmax <= jenkins_fold_kjoin){
      printf("\n\nLoading data without folding");

      // Ignore reference to mock, should work just fine. 
      load_withoutfolding(filepath);
    } 

    else{
      printf("\n\nLoading data with folding");

      // Ignore reference to mock, should work just fine.
      load_withfolding(filepath);
    }

    if(fieldFlag == 1){
      if(appliedClippingThreshold      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 335.538;  // Correct monopole.                                                                                      
      else if(appliedClippingThreshold ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 283.993;  // Note: as a difference estimator, covariance does not have to explicity corrected.                       
      else if(appliedClippingThreshold ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 229.602;
      else if(appliedClippingThreshold ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 180.127;
    }

    if(fieldFlag == 4){
      if(appliedClippingThreshold      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 331.20;  // Correct monopole.                                                                                                                     
      else if(appliedClippingThreshold ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 278.29;  // Note: as a difference estimator, covariance does not have to explicity corrected.                                                     
      else if(appliedClippingThreshold ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 227.51;
      else if(appliedClippingThreshold ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 177.75;
    }

    // Amplitude corrected in ChiSq_Eval. 
    // if(appliedClippingThreshold      == 1000.)  for(j=0; j<order; j++) xdata[j]      *=   1.00;  // Correct monopole.                                                                                                                    
    // else if(appliedClippingThreshold ==   10.)  for(j=0; j<order; j++) xdata[j]      *=   1.30;  // Note: as a difference estimator, covariance does not have to explicity corrected.                                                    
    // else if(appliedClippingThreshold ==    6.)  for(j=0; j<order; j++) xdata[j]      *=   1.70;

    /*
    if(ChiSq_kmax >= 0.8){
      // write combined unfolded+folded p(k).
      sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk_withfolding/old_d0_%d/W%d/data_loz_%.1lf_hiz_%.1lf.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, lo_zlim, hi_zlim);

      output = fopen(filepath, "w");

      for(j=0; j<mono_order; j++){
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], xdata[j], 1./gsl_matrix_get(sigma_norm, j, j), xdata[j + mono_order], 1./gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order));
      }

      fclose(output);
    }
    */
    return 0;
}


int ChiSq_minimisation(int mockNumber){    
    int ll, mm;
    
    if(data_mock_flag == 0)  load_mock(mockNumber);
    
    else if(data_mock_flag == 1)  load_data();
    
    
    else{
      // scale_Cov();  

      // load_meanmultipoles();

      // load_aptest_multipoles();
    }

    // new zero mean, unit variance, decorrelated variables. 
    for(j=0; j<order; j++){
        ydata[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
    }
    
    // printf("\n\nxdata: \t evec 1: \t diag sigma: \t ydata:\n");
    // for(j=0; j<order; j++)  printf("%le \t %le \t %le \t %le \n", xdata[j], gsl_vector_get(col, j), gsl_matrix_get(sigma_norm, j, j), ydata[j]);
    // input_check();
    // prep_dlnPR_dlnk();
    
    printf("\n\n\nBeginning Chi sq. minimisation.");
  
    printf("\n\nPriors:");
    printf("\n%.2f < fsigma8       < %.2f",      min_fsigma8,      max_fsigma8);  
    printf("\n%.2f < vel. disp.    < %.2f",      min_velDisperse,  max_velDisperse);
    printf("\n%.2f < bsigma8       <  %.2f",     min_bsigma8,      max_bsigma8);
    printf("\n%.2f < alpha   pad.  < %.2f",      min_alpha_pad,    max_alpha_pad);
    printf("\n%.2f < epsilon pad.  < %.2f",      min_epsilon_pad,  max_epsilon_pad);
    printf("\n%.2f < A11 sq.       < %.2f",      min_A11Sq,        max_A11Sq);
    printf("\n");
    
    double     fsigma8Interval  = (max_fsigma8     - min_fsigma8)/dRes;
    double     bsigma8Interval  = (max_bsigma8     - min_bsigma8)/dRes;
    double       sigmaInterval  = (max_velDisperse - min_velDisperse)/dRes; 
    double  alpha_padInterval   = (max_alpha_pad   - min_alpha_pad)/dRes_ap;
    double  epsilon_padInterval = (max_epsilon_pad - min_epsilon_pad)/dRes_ap;
    double       A11SqInterval  = (max_A11Sq       - min_A11Sq)/dRes;
    
    // mono and quad. fit
    double dof = order - paramNumber;
    
    ChiSq_expected = dof; 

    printf("\n\nexpected chi sq. %.2f +- %.2f\n\n", dof, sqrt(2.*dof));
        
    // Alocock-Paczynski, epsilon dependent correction terms.
    double mono_epsilonCorrection_1, mono_epsilonCorrection_2;
    double quad_epsilonCorrection_1, quad_epsilonCorrection_2;
    
    printf("\nf\sig_8  \sig_v  b\sig_8  A_11  u0 alpha_p  eps_p  min \chi^2\n");    
    
    //if(data_mock_flag == 0){    
    //  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/max_likelihoodChains/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_kmax_%.1lf.dat", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber, lo_zlim, hi_zlim, ChiSq_kmax);
    //}

    //else if(data_mock_flag == 1){
    //  sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/max_likelihoodChains/d0_%d/W%d/data_zlim_%.1lf_%.1lf_kmax_%.1lf.dat", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber, lo_zlim, hi_zlim, ChiSq_kmax);
    //}

    //else{
    //  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/max_likelihoodChains/d0_%d/W%d/supervipers_zlim_%.1lf_%.1lf_kmax_%.1lf.dat", root_dir, (int) round(appliedClippingThreshold), fieldFlag, mockNumber, lo_zlim, hi_zlim, ChiSq_kmax);
    //}

    // output = fopen(filepath, "w");
    
    // time_t mytime;
    
    for(jj=0; jj<Res; jj++){    
        fsigma8 = min_fsigma8 + fsigma8Interval*jj;

        // mytime  = time(NULL);
        
        // printf("\n\n%s", ctime(&mytime));

        for(kk=0; kk<Res; kk++){
            bsigma8 = min_bsigma8 + bsigma8Interval*kk;

            for(ii=0; ii<Res; ii++){
                velDispersion = min_velDisperse + sigmaInterval*ii;
                
                for(ll=0;ll<Res_ap; ll++){
                    alpha_pad = min_alpha_pad + alpha_padInterval*ll;
                
                    // FFTLog_setInput(mono_config, fsigma8/bsigma8, velDispersion);  // Can set Res_ap = 1, and alpha_pad = 1.0, epsilon_pad = 0.0 to neglect AP.
                    // FFTLog_setInput(quad_config, fsigma8/bsigma8, velDispersion);  // Can set Res_ap = 1, and alpha_pad = 1.0, epsilon_pad = 0.0 to neglect AP.
                
                    // Alcock - Paczynski corrected multipoles, stripped of the epsilon 'corrections'.
                    // FFTLog_updateInput_monoAP(mono_config, fsigma8/bsigma8, velDispersion, alpha_pad, &mono_epsilonCorrection_1, &mono_epsilonCorrection_2);
                    // FFTLog_updateInput_quadAP(quad_config, fsigma8/bsigma8, velDispersion, alpha_pad, &quad_epsilonCorrection_1, &quad_epsilonCorrection_2);
                
                    for(mm=0; mm<Res_ap; mm++){
                        epsilon_pad = min_epsilon_pad + epsilon_padInterval*mm;
                    
			alpha_pad   = 1.0;
			epsilon_pad = 0.0; 

                        // Alcock - Paczynski corrected multipoles, add the epsilon 'corrections'.
                        // FFTLog_updateInput_monoAP_epsiloncorrections(mono_config_ap, mono_config, epsilon_pad, mono_epsilonCorrection_1, mono_epsilonCorrection_2);
                        // FFTLog_updateInput_quadAP_epsiloncorrections(quad_config_ap, quad_config, epsilon_pad, quad_epsilonCorrection_1, quad_epsilonCorrection_2);
                    
                        // no marginalisation over clipping. 
                        A11Sq                                 =                                          1.0;
                                    
                        ChiSqGrid[jj][kk][ii][ll][mm]         =                                  ChiSqEval();
                    
                        lnLikelihoodGrid[jj][kk][ii][ll][mm]  =           -0.5*ChiSqGrid[jj][kk][ii][ll][mm]; 
                    
                        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mocks_fsig8/no_clipping/kmax_%.1lf/no_marginalisation/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_fsig8.dat", 
                        // root_dir,     ChiSq_kmax, mockNumber);
    
                        // output = fopen(filepath, "w");
            
                        // printf("%e \t %e \t %e \t %e \t %e \t %e \n", fsigma8, velDispersion, bsigma8, A11Sq, u0, ChiSqGrid[jj][kk][ii][ll][mm]);
    
                        // fclose(output);
              
                        if(ChiSqGrid[jj][kk][ii][ll][mm] < minChiSq){
                            minChiSq             = ChiSqGrid[jj][kk][ii][ll][mm];
                    
                            minChiSq_fsigma8     = fsigma8;
                            minChiSq_A11Sq       = A11Sq;
                            minChiSq_sigma       = velDispersion;
                            minChiSq_bsigma8     = bsigma8;
                            minChiSq_alpha_pad   = alpha_pad;
                            minChiSq_epsilon_pad = epsilon_pad;
                    
                            // fprintf(output, "%.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \n", fsigma8, velDispersion, bsigma8, A11Sq, 
			    //	                                                                                                u0, alpha_pad, epsilon_pad, minChiSq);

			    printf("%.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \t %.2lf \n", fsigma8, velDispersion, bsigma8, A11Sq, u0, alpha_pad, epsilon_pad, minChiSq);
                        }
	            }
            }
        }
      }
    }
    
    // fclose(output);
    
    fsigma8       = minChiSq_fsigma8;
    A11Sq         = minChiSq_A11Sq;
    velDispersion = minChiSq_sigma;
    bsigma8       = minChiSq_bsigma8;
    epsilon_pad   = 0.0;
    alpha_pad     = 1.0; 
    
    ChiSqEval();
    
    printf("\n\nBest fit theory: ");
    
    double mono_chi2, quad_chi2;
    double tot_chi2 = 0.0;

    for(i=0; i<mono_order; i++){  
      mono_chi2 = pow(xdata[i]-xtheory[i], 2.)*pow(gsl_matrix_get(sigma_norm, i, i), 2.);

      quad_chi2 = pow(xdata[i+mono_order]-xtheory[i+mono_order], 2.)*pow(gsl_matrix_get(sigma_norm, i+mono_order, i+mono_order), 2.);

      tot_chi2 += mono_chi2 + quad_chi2; 	

      printf("\n%e \t %e \t %e \t %e \t %e \t %e \t %e \t %.2lf \t %.2lf \t %.2lf", kVals[i], xdata[i], xdata[mono_order + i], xtheory[i], xtheory[i+mono_order], 1./gsl_matrix_get(sigma_norm, i, i), 1./gsl_matrix_get(sigma_norm, i+mono_order, i+mono_order), mono_chi2, quad_chi2, tot_chi2);
    
      // fprintfBestfit();
    }
    return  0;
}


int W1_Spectro_V7_2_chiSq_minimisation(){
  double maxlike_fsig8, maxlike_sigv, maxlike_bsig8;
  
  prepChiSq_minimisation(); 
  
  // assigns memory for chi sq. grid, likelihood grid, and 1 and 2 parameter posteriors.                                     
  LikelihoodMemory();
  
  // low, but sufficient resolution. 
  FFTlogRes  =     768;
  // FFTlogRes  = 4096;
  
  pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

  // Must have mask multipoles available, uncomment prep_VIPERS_maskMultipoles().                                            
  FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);
  
  // connect measured k vals with their nearest FFTlog counterpart.                                                          
  kvals_matchup();

  minChiSq = pow(10., 12.);

  /*
  printf("\n\nBeginning MCMC.");

  // For AP correction.  
  prep_dlnPR_dlnk();

  metropolis_mcmc();
  */
  
  ChiSq_minimisation(666);

  double fsig8_error;
  
  maxlike_fsig8 = calc_fsigma8Posterior();        
  maxlike_sigv  = calc_velDispPosterior();        
  maxlike_bsig8 = calc_bsigma8Posterior();
  
  // W1_Spectro_V7_3 -> Clipping
  sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/data_%.1lf_%.1lf_snipCorrected.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, ChiSq_kmax, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", maxlike_fsig8,      maxlike_sigv,    maxlike_bsig8, ChiSq_expected);
  fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", minChiSq_fsigma8, minChiSq_sigma, minChiSq_bsigma8,       minChiSq);

  fclose(output);

  // printf("%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", maxlike_fsig8, fsig8_error, maxlike_sigv, maxlike_bsig8, ChiSq_expected, minChiSq);

  // calc_fsigma8_bsigma8_posterior();

  /*
  sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/data_%.1lf_%.1lf.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, ChiSq_kmax, lo_zlim, hi_zlim);
    
  output = fopen(filepath, "w");

  fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", maxlike_fsig8, fsig8_error, maxlike_sigv, maxlike_bsig8, ChiSq_expected, minChiSq);
  
  fclose(output);
  */
  return 0;
}


int ensemble_fsig8(int mockNumber, int start, int totalCats){
    int jjj, rand;
    
    double maxlike_fsig8, maxlike_sigv, maxlike_bsig8;
    
    prepChiSq_minimisation();
    
    // assigns memory for chi sq. grid, likelihood grid, and 1 and 2 parameter posteriors. 
    LikelihoodMemory();
    
    // prep. clipping correction. 
    FFTlogRes  = 768;
    // FFTlogRes = 4096;

    pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

    // Must have mask multipoles available, uncomment prep_VIPERS_maskMultipoles().
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);
    
    // connect measured k vals with their nearest FFTlog counterpart. 
    kvals_matchup();    
    
    for(jjj=0; jjj<mockNumber; jjj++){    
        minChiSq = pow(10., 12.);
    
        // Random sub-selection of mocks. 
        // rand     = gsl_rng_uniform_int(gsl_ran_r, totalCats);
            
        // printf("\n\nReducing mock %d", start+rand);
	
        ChiSq_minimisation(start+jjj);
	
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mocks_fsig8/no_clipping/kmax_%.1lf/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_fsig8.dat", root_dir, ChiSq_kmax, start + rand);
        
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mocks_fsig8/linear_matterpk/clipping/kmax_%.1lf/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_fsig8.dat", root_dir, ChiSq_kmax, 
        // start + rand);
        
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mocks_fsig8/linear_rsd/clipping/kmax_%.1lf/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_fsig8.dat", root_dir, ChiSq_kmax, 
        // start+rand);
        
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mocks_fsig8/linear_matterpk/clipping_d0_6/kmax_%.1lf/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_fsig8.dat", root_dir, 
        // ChiSq_kmax, start+rand);
	
        maxlike_fsig8 = calc_fsigma8Posterior();        
        maxlike_sigv  = calc_velDispPosterior();
        maxlike_bsig8 = calc_bsigma8Posterior();
	
        //-- new method of clipped p(k) calculation. Shot noise does not scale with A11.  mocks for W1, (v8 500s) + Nagoya v6 + Samhain + specweight + fkpweight. --//
        // sprintf(filepath, "%s/W1_Spectro_V7_1/mocks_fsig8/d0_1000/kmax_%.1lf/mock_%d_256_pk_d0_1000.00_fsig8.dat", root_dir, ChiSq_kmax, start+rand);
        // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mock_%d_%.1lf_%.1lf.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, ChiSq_kmax, start+jjj, lo_zlim, hi_zlim);
	// sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/fsig8/true_nbar/d0_%d/W%d/kmax_%.1lf/mock_%d_%.1lf_%.1lf.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, ChiSq_kmax, start+jjj, lo_zlim, hi_zlim);
        sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mock_%d_%.1lf_%.1lf_snipCorrected.dat", root_dir, (int) ceil(appliedClippingThreshold), fieldFlag, ChiSq_kmax, start+jjj, lo_zlim, hi_zlim);

        output = fopen(filepath, "w");

        fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", maxlike_fsig8,      maxlike_sigv,    maxlike_bsig8, ChiSq_expected);
	fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", minChiSq_fsigma8, minChiSq_sigma, minChiSq_bsigma8,       minChiSq);

        fclose(output);
    }
    
    return 0;
} 


double ChiSqEval_ap(){
  // Chi sq. evaluation, accounting for ap scaling between model and fiducial cosmology. 
  double dummy        = 0.0;
  double fieldArea    = 0.0;
  double ChiSq        = 0.0;
  double cnvldpk_zero = 0.0;

  // FFTLog_setInput(mono_config, fsigma8/bsigma8, velDispersion);
  // FFTLog_setInput(quad_config, fsigma8/bsigma8, velDispersion);
  // FFTLog_setInput( hex_config, fsigma8/bsigma8, velDispersion);
  
  // Alcock - Paczynski corrected multipoles.
  // NOTE:  P(k) model in AP_P0 & AP_P2 must match that in FFTlog.c SetInput().
  // FFTLog_updateInput_monoAP(mono_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);                                                                     
  // FFTLog_updateInput_quadAP(quad_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);                                                                      
  // FFTLog_updateInput_hexAP(  hex_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad); 
  
  FFTLog_updateInput_multipolesAP(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);
  
  // Transform to correlation functions, mono and quad.                                                                                                                                                          
  xi_mu(mono_config);
  xi_mu(quad_config);
  xi_mu( hex_config);
  
  // window multipoles set to joint field by FFTLog_setInput.                                                                                                                                                     
  pt2maskMultipoles = &splint_VIPERS_jmaskMultipoles;

  // calls (*pt2maskMultipoles) directly.
  cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);

  pk_mu(convlmonoCorr);

  for(j=0; j<mono_config->N; j++){
    if(convlmonoCorr->krvals[j][0] > 0.001){
      cnvldpk_zero = convlmonoCorr->pk[j][0];

      break;
    }
  }
  
  pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

  cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  
  pk_mu(convlmonoCorr);
  pk_mu(convlquadCorr);
    
  if(fieldFlag == 1) fieldArea = W1area;
  if(fieldFlag == 4) fieldArea = W4area;

  for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] -= cnvldpk_zero*(fieldArea/TotalW1W4area)*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
  for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] -= cnvldpk_zero*(fieldArea/TotalW1W4area)*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);
  
  for(j=0; j<mono_order; j++)  xtheory[j]                   =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];                                                            
  for(j=0; j<mono_order; j++)  xtheory[j + mono_order]      =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];                                                            
  
  // for(j=0; j<mono_order; j++)  xtheory[j]                   =  AP_P0(mono_config->krvals[j][0], fsigma8/bsigma8, velDispersion, epsilon_pad, alpha_pad);
  // for(j=0; j<mono_order; j++)  xtheory[j + mono_order]      =  AP_P2(mono_config->krvals[j][0], fsigma8/bsigma8, velDispersion, epsilon_pad, alpha_pad);

  // for(j=0; j<mono_order; j++){  
  //  xtheory[j]                   =  mono_config->pk[fftlog_indices[j]][0];
  //  xtheory[j + mono_order]      =  quad_config->pk[fftlog_indices[j]][0];

    // printf("\n%.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf", xtheory[j], xtheory[j + mono_order], xdata[j], xdata[j + mono_order], gsl_matrix_get(sigma_norm, j, j), gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.));
  //}

  // new decorrelated, unit variance, zero mean variables.                                                                                                                                                       
  for(j=0; j<order; j++){
    ytheory[j] = 0.0;

    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++)  ytheory[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];
  }
  
  for(j=0; j<order; j++)  ChiSq += pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j);
  
  // for(j=0; j<order; j++)  ChiSq += pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.);
  
  return ChiSq;
}


double ChiSqEval(){
    double fieldArea    = 0.0;

    double ChiSq        = 0.0;
    double cnvldpk_zero = 0.0;

    // printf("\n%e \t %e", fsigma8, bsigma8);

    // A11Sq               =       1.;
    // fsigma8             = 0.9*0.79;
    // bsigma8             = 1.0*0.79;
    // velDispersion       =       3.;
    
    FFTLog_setInput(mono_config, fsigma8/bsigma8, velDispersion);
    FFTLog_setInput(quad_config, fsigma8/bsigma8, velDispersion);
    FFTLog_setInput( hex_config, fsigma8/bsigma8, velDispersion);

    // Transform to correlation functions, mono and quad.
    xi_mu(mono_config);
    xi_mu(quad_config);
    xi_mu( hex_config);
    
    varCalc(mono_config, &variance);

    if(data_mock_flag == 0){
      if(appliedClippingThreshold      == 1000.)  A11Sq  =      1.00;
      else if(appliedClippingThreshold ==   10.)  A11Sq  =   1./1.30; // calculated from mock mean or per mock estimate?
      else if(appliedClippingThreshold ==    6.)  A11Sq  =   1./1.85;
      else if(appliedClippingThreshold ==    4.)  A11Sq  =   1./3.00;
    }
    
    if(data_mock_flag == 1){
      if(fieldFlag == 1){
        if(appliedClippingThreshold      == 1000.)  A11Sq  =      1.00;  // Correct monopole.                                                                                                                     
        else if(appliedClippingThreshold ==   10.)  A11Sq  =   1./1.30;  // Note: as a difference estimator, covariance does not have to explicity corrected.                                                     
        else if(appliedClippingThreshold ==    6.)  A11Sq  =   1./1.70;
        else if(appliedClippingThreshold ==    4.)  A11Sq  =   1./2.70;  
      }

      if(fieldFlag == 4){
        if(appliedClippingThreshold      == 1000.)  A11Sq  =      1.00;  // Correct monopole.                                                                                                                                               
        else if(appliedClippingThreshold ==   10.)  A11Sq  =   1./1.30;  // Note: as a difference estimator, covariance does not have to explicity corrected.                                                                               
        else if(appliedClippingThreshold ==    6.)  A11Sq  =   1./1.90;
        else if(appliedClippingThreshold ==    4.)  A11Sq  =   1./2.70;
      }
    }

    // u0 calc. in the absence of the value of \delta_0.                                                                                                                                                                                    
    u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);

    // printf("\nu0: %e", u0);
                                                                                                                                                                                                                                            
    clipmono(clipmono_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
    clipquad(clipquad_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
    
    if(appliedClippingThreshold < 1000.){     
      for(j=0; j<mono_config->N; j++){
        mono_config->xi[j][0] = clipmono_config->xi[j][0];
        quad_config->xi[j][0] = clipquad_config->xi[j][0];
      }
    }

    pt2maskMultipoles = &splint_VIPERS_jmaskMultipoles;

    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    
    pk_mu(convlmonoCorr);

    // zero point determined for P'(0) for the joint fields. 
    for(j=0; j<mono_config->N; j++){
      if(convlmonoCorr->krvals[j][0] > 0.001){
        cnvldpk_zero = convlmonoCorr->pk[j][0];

        break;
      }
    }

    pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);

    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);

    if(fieldFlag == 1) fieldArea = W1area;
    if(fieldFlag == 4) fieldArea = W4area;

    for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] -= cnvldpk_zero*(fieldArea/TotalW1W4area)*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] -= cnvldpk_zero*(fieldArea/TotalW1W4area)*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    /*    
    // AP correction. 
    // xi_mu(mono_config_ap);
    // xi_mu(quad_config_ap);

    // varCalc(mono_config_ap, &variance);
    
    // u0 calc. in the absence of the value of \delta_0.
    // u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);
    // u0 = appliedClippingThreshold/sqrt(2.*variance);
    
    // printf("\nu0: %e", u0);

    // Currently evaluated at second order. 
    // clipmono(clipmono_config, mono_config, quad_config, zero_config, u0, variance);
    // clipquad(clipquad_config, mono_config, quad_config, zero_config, u0, variance);
    
    // could include hex terms of correlation fn. 
    // cnvldmonoCorr(convlmonoCorr, clipmono_config, clipquad_config, zero_config, zero_config, zero_config, zero_config);
    // cnvldquadCorr(convlquadCorr, clipmono_config, clipquad_config, zero_config, zero_config, zero_config, zero_config);
    
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    
    // pk_mu(clipmono_config);
    // pk_mu(clipquad_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    
    for(j=0; j<mono_config->N; j++){
      if(convlmonoCorr->krvals[j][0] > 0.001){
	cnvldpk_zero = convlmonoCorr->pk[j][0];

	break;
      }
    }

    //for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = convlmonoCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    //for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = convlquadCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);
    */

    for(j=0; j<mono_order; j++)  xtheory[j]                   =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    for(j=0; j<mono_order; j++)  xtheory[j + mono_order]      =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
    
    // new decorrelated, unit variance, zero mean variables. 
    for(j=0; j<order; j++){
        ytheory[j] = 0.0;
        
        gsl_matrix_get_col(col, evec, j);
        
        for(k=0; k<order; k++)  ytheory[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];
    }
    
    for(j=0; j<order; j++)  ChiSq += pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j); 
    
    // for(j=0; j<order; j++)  ChiSq += pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.); 
    
    return ChiSq;
}


int fprintfBestfit(){
    fsigma8         =     minChiSq_fsigma8;
    bsigma8         =     minChiSq_bsigma8;
    velDispersion   =       minChiSq_sigma;
    alpha_pad       =   minChiSq_alpha_pad;
    epsilon_pad     = minChiSq_epsilon_pad;
    
    // A11Sq        = 1.;
    // A11Sq        = minChiSq_A11Sq;

    printf("\n\nBest fit parameters, fsigma8: %.4e, vel. disp.: %.4e, bsigma8: %.4e, A11: %.4e", fsigma8, velDispersion, bsigma8, A11Sq);

    // Alocock-Paczynski, epsilon dependent correction terms.
    // double mono_epsilonCorrection_1, mono_epsilonCorrection_2;
    // double quad_epsilonCorrection_1, quad_epsilonCorrection_2;

    // Alcock - Paczynski corrected multipoles, stripped of the epsilon 'corrections'.
    // FFTLog_updateInput_monoAP(mono_config, fsigma8/bsigma8, velDispersion, alpha_pad, &mono_epsilonCorrection_1, &mono_epsilonCorrection_2);
    // FFTLog_updateInput_quadAP(quad_config, fsigma8/bsigma8, velDispersion, alpha_pad, &quad_epsilonCorrection_1, &quad_epsilonCorrection_2);
          
    // Alcock - Paczynski corrected multipoles, add the epsilon 'corrections'.
    // FFTLog_updateInput_monoAP_epsiloncorrections(mono_config_ap, mono_config, epsilon_pad, mono_epsilonCorrection_1, mono_epsilonCorrection_2);
    // FFTLog_updateInput_quadAP_epsiloncorrections(quad_config_ap, quad_config, epsilon_pad, quad_epsilonCorrection_1, quad_epsilonCorrection_2);

    ChiSqEval();

    // Best fit correlated variables. 
    // sprintf(filepath, "%s/Data/500s/clipped_bestfit_correlatedMonopole_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    // sprintf(filepath, "%s/W1_Spectro_V5_0/bestfit_Monopole_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/chisq_bestfit_multipoles/mocks_W1_v8.0_500_194_Nagoya_v6_Samhain_specweight_fkpweight_clipped_bestfit_mono_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    sprintf(filepath, "%s/W1_Spectro_V7_1/chisq_bestfit_multipoles_mono_kmax_%.2lf.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_order; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j], 
                                                         xtheory[j], 
                                                         1./gsl_matrix_get(sigma_norm, j, j),  
                                                         pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.));
    
    fclose(output); 
    
    // sprintf(filepath, "%s/Data/500s/clipped_bestfit_correlatedQuadrupole_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    // sprintf(filepath, "%s/W1_Spectro_V5_0/bestfit_Quadrupole_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/chisq_bestfit_multipoles/mocks_W1_v8.0_500_194_Nagoya_v6_Samhain_specweight_fkpweight_clipped_bestfit_quad_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    sprintf(filepath, "%s/W1_Spectro_V7_1/chisq_bestfit_multipoles_quad_kmax_%.2lf.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_order; j++) 
        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                         xdata[j + mono_order], 
                                                         xtheory[j + mono_order], 
                                                         1./gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), 
                                                         pow(xdata[j+mono_order] - xtheory[j+mono_order], 2.)*pow(gsl_matrix_get(sigma_norm, j+mono_order, j+mono_order), 2.));
    
    fclose(output);
    
    /******************/
    /*
    fsigma8         =                  0.5;
    bsigma8         =                  0.8;
    velDispersion   =                  5.0;
    alpha_pad       =   minChiSq_alpha_pad;
    epsilon_pad     = minChiSq_epsilon_pad;

    ChiSqEval();
    
    sprintf(filepath, "%s/W1_Spectro_V7_1/theory_bestfit_multipoles_hex.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_order; j++){ 
        fprintf(output, "%e \t %e \t %e \t %e \n", mono_config->krvals[fftlog_indices[j]][0], 
                                                       mono_config->pk[fftlog_indices[j]][0], 
                                                       quad_config->pk[fftlog_indices[j]][0], 
                                                       hex_config->pk[fftlog_indices[j]][0]);
    }
    
    fclose(output);
    */
    return 0;
}


int fprintf_model(){
    double cnvldpk_zero = 0.0; 
  
    fsigma8             = 0.9*0.79;
    bsigma8             = 1.0*0.79;
    velDispersion       =       3.;
    
    A11Sq               = 1.;

    FFTlogRes           = 768;

    pt2maskMultipoles   = &splint_VIPERS_maskMultipoles;

    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);
    
    FFTLog_setInput(mono_config, fsigma8/bsigma8, velDispersion);
    FFTLog_setInput(quad_config, fsigma8/bsigma8, velDispersion);
    FFTLog_setInput( hex_config, fsigma8/bsigma8, velDispersion);

    xi_mu(mono_config);                                                                                                                                                                      
    xi_mu(quad_config);                                                                                                                                                                      
    xi_mu( hex_config);                                                                                                                                                                      

    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, zero_config, zero_config, zero_config, zero_config);                                                                              
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, zero_config, zero_config, zero_config, zero_config);                                                                                                                                                                                                                                                      
    pk_mu(convlmonoCorr);                                                                                                                                                                    
    pk_mu(convlquadCorr);                                                                      

    for(j=0; j<mono_config->N; j++){                                                                                                                                                        
      if(convlmonoCorr->krvals[j][0] > 0.001){                                                                                                                                             
        cnvldpk_zero = convlmonoCorr->pk[j][0];                                                                                                                                              
                                                                                                                                                                                             
        break;                                                                                                                                                                               
      }                                                                                                                                                                                      
    }                                                                                                                                                                                       
                                                                                                                                                                                             
    for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] -= cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);                              
    for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] -= cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    printf("\n\nCnvld. P(k) zero point: %le", cnvldpk_zero);

    // Best fit correlated variables. 
    sprintf(filepath, "%s/W1_Spectro_V7_2/model_multipoles_quad_kmax_%.2lf.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(j=0; j<FFTlogRes; j++){  
      if((mono_config->krvals[j][0] > 0.01) && (mono_config->krvals[j][0] < 1.)){
        fprintf(output, "%e \t %e \t %e \n", convlmonoCorr->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
      }
    }
    
    fclose(output); 
    
    return 0;
}
