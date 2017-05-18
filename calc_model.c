int model_compute(int aa, int bb, int cc, int dd, int ee){
  double ChiSq        = 0.0;
  double fieldArea    = 0.0;
  double cnvldpk_zero = 0.0;
 
  FFTlog_updatepk(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion);
  // apmultipoles(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);

  for(j=0; j<mono_config->N; j++){
    if((0.01 < mono_config->krvals[j][0]) && (mono_config->krvals[j][0] < 1.0)){
      printf("\n%.4lf \t %.4lf \t %.4lf", mono_config->krvals[j][0], mono_config->pk[j][0], quad_config->pk[j][0]);
    }
  }
  
  xi_mu(mono_config);  // Transform to correlation function. 
  xi_mu(quad_config);  
  xi_mu( hex_config);
  
  varCalc(mono_config, &variance);
  
  if(d0 < 1000){
    clip_p0p2(clipmono_config, clipquad_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
  }
  
  cnvldmonoCorr_joint(convlmonoCorr, mono_config, quad_config, hex_config);
  
  pk_mu(convlmonoCorr);
  
  cnvldpk_zero = convlmonoCorr->pk[cnvldpk_zeropoint_index][0];  // zero point determined for P'(0) for the joint fields. 
  
  cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  
  pk_mu(convlmonoCorr);
  pk_mu(convlquadCorr);
  
  // double kmask_norm = get_kMask_norm();  
  // printf("\n\nKMASK NORM RATIO: %.9lf", kmask_norm*9.599*pow(10., -2.)/(cnvldpk_zero*fracArea));
  
  for(j=0; j<mono_config->N; j++){
    convlmonoCorr->pk[j][0] -= cnvldpk_zero*fracArea*FFTlog_Wk0[j];
    convlquadCorr->pk[j][0] -= cnvldpk_zero*fracArea*FFTlog_Wk2[j];

    // printf("\n%.6lf \t %.6lf", convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
  }
  
  printf("\n\nModel eval.");
  
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
