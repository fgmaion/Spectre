int model_compute(int aa, int bb, int cc, int dd, int ee, int print){
  double ChiSq        = 0.0;
  double fieldArea    = 0.0;
  double cnvldpk_zero = 0.0;
 
  FFTlog_updatepk(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion);
  // apmultipoles(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);
  /*
  for(j=0; j<mono_config->N; j++){
    if((0.01 < mono_config->krvals[j][0]) && (mono_config->krvals[j][0] < 1.0)){
      printf("\n%le \t %le \t %le", mono_config->krvals[j][0], mono_config->pk[j][0], quad_config->pk[j][0]);
    }
  }
  */
  xi_mu(mono_config);  // Transform to correlation function. 
  xi_mu(quad_config);  
  xi_mu( hex_config);
  
  varCalc(mono_config, &variance);
  /*
  if(d0 < 1000){  // In fact, sets mono_config by clipmono_config. 
    clip_p0p2(clipmono_config, clipquad_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
  }
  */
  // last arg. for shot noise contribution to P_*(0) as given by shot_config 
  cnvldmonoCorr_joint(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, 0);
  
  pk_mu(convlmonoCorr);
  
  cnvldpk_zero = convlmonoCorr->pk[cnvldpk_zeropoint_index][0];  // zero point determined for P'(0) for the joint fields. 
  
  cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  
  pk_mu(convlmonoCorr);
  pk_mu(convlquadCorr);
  
  // for(j=0; j<mono_config->N; j++){
  //  if((0.02 < mono_config->krvals[j][0]) && (mono_config->krvals[j][0] < 0.8)){
  //    printf("\n%lf \t %lf \t %lf", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
  //  }
  //}
  
  // double kmask_norm = get_kMask_norm();  
  // printf("\n\nKMASK NORM RATIO: %.9lf", kmask_norm*9.599*pow(10., -2.)/(cnvldpk_zero*fracArea));
  
  for(j=0; j<mono_config->N; j++){
    // convlmonoCorr->pk[j][0]  = mono_config->pk[j][0];
    // convlquadCorr->pk[j][0]  = quad_config->pk[j][0];

    // Single-field
    // convlmonoCorr->pk[j][0] -= cnvldpk_zero*FFTlog_Wk0[j];
    // convlquadCorr->pk[j][0] -= cnvldpk_zero*FFTlog_Wk2[j];

    // Joint-field 
    convlmonoCorr->pk[j][0]    -= cnvldpk_zero*fracArea*FFTlog_Wk0[j];
    convlquadCorr->pk[j][0]    -= cnvldpk_zero*fracArea*FFTlog_Wk2[j];
  }
  
  if(print == 1){
    print_basemodel();

    printf("\n\nConvolved P(k) zero point: %.4lf", cnvldpk_zero);  
  }
  
  return 0;
}

int print_basemodel(){
  printf("\n\nConvolved power spectra: \n\n");

  sprintf(filepath, "%s/models/defaultparamas_model_intcor_cnvld_W%d_zlim_%.1lf_%.1lf_d0_%d_kmax_%.1lf.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, d0, ChiSq_kmax);

  output = fopen(filepath, "w");

  for(j=0; j<mono_config->N; j++){
    if((0.01 < mono_config->krvals[j][0]) && (mono_config->krvals[j][0] < 1.0)){
      printf("%le \t %le \t %le \n", convlmonoCorr->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);

      fprintf(output, "%le \t %le \t %le \n", convlmonoCorr->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
    }
  }

  fclose(output);

  return 0;
}

int ytheory_compute(int aa, int bb, int cc, int dd, int ee){
  for(j=0; j<order; j++){
    ytheory[aa][bb][cc][dd][ee][j] = 0.0;

    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++){
      ytheory[aa][bb][cc][dd][ee][j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[aa][bb][cc][dd][ee][k];
    }
  }

  return 0;
}

int benchmark_ChiSq_eval(){
  double begin  = getRealTime();
  double end    = getRealTime();

  printf("\n\nChi Sq. benchmark (1000 runs): %.6lf", (end - begin));
  
  return 0;
}
