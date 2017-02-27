int model_compute ChiSqEval(int aa, int bb, int cc, int dd, int ee){
  double ChiSq        = 0.0;
  double fieldArea    = 0.0;
  double cnvldpk_zero = 0.0;

  FFTlog_updatepk(mono_config, fsigma8/bsigma8, velDispersion); // zero'd imag. of p(k).
  FFTlog_updatepk(quad_config, fsigma8/bsigma8, velDispersion);
  FFTlog_updatepk( hex_config, fsigma8/bsigma8, velDispersion);

  xi_mu(mono_config);  // Transform to correlation function. 
  xi_mu(quad_config);  
  xi_mu( hex_config);

  varCalc(mono_config, &variance);

  // Can openmp loops. mono_config updated by clip_mono.
  clip_p0p2(clipmono_config, clipquad_config, mono_config, quad_config, zero_config, zero_config, u0, variance);

  cnvldmonoCorr_joint(convlmonoCorr, mono_config, quad_config, hex_config);

  pk_mu(convlmonoCorr);

  cnvldpk_zero = convlmonoCorr->pk[cnvldpk_zeropoint_index][0];  // zero point determined for P'(0) for the joint fields. 
                                                 
  // Can openmp loops. 
  cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
  cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);

  pk_mu(convlmonoCorr);
  pk_mu(convlquadCorr);

  for(j=0; j<mono_config->N; j++){
    convlmonoCorr->pk[j][0] -= cnvldpk_zero*fracArea*FFTlog_Wk0[j];
    convlquadCorr->pk[j][0] -= cnvldpk_zero*fracArea*FFTlog_Wk2[j];
  }

  for(j=0; j<mono_order; j++){
    xtheory[aa][bb][cc][dd][ee][j]              =  convlmonoCorr->pk[fftlog_indices[j]][0]; // clipmono_config->pk[fftlog_indices[j]][0];
    xtheory[aa][bb][cc][dd][ee][j + mono_order] =  convlquadCorr->pk[fftlog_indices[j]][0]; // clipquad_config->pk[fftlog_indices[j]][0];
  }
  
  for(j=0; j<order; j++){
    ytheory[j] = 0.0;

    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++)  ytheory[aa][bb][cc][dd][ee][j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];

  }

  return ChiSq;
}


int benchmark_ChiSq_eval(){
  double begin  = getRealTime();

  for(jj=0; jj<1000; jj++)  ChiSqEval();
  
  double end    = getRealTime();

  printf("\n\nChi Sq. benchmark (1000 runs): %.6lf", (end - begin));
  
  return 0;
}


double ChiSqEval_ap(){
  // Chi sq. evaluation, accounting for AP rescaling between model and fiducial cosmology.
  double ChiSq        = 0.0;
  double dummy        = 0.0;
  double fieldArea    = 0.0;
  double cnvldpk_zero = 0.0;

  FFTLog_updateInput_multipolesAP(mono_config, quad_config, hex_config, fsigma8/bsigma8, velDispersion, alpha_pad, epsilon_pad);

  xi_mu(mono_config);
  xi_mu(quad_config);
  xi_mu( hex_config);

  pt2maskMultipoles = &splint_VIPERS_jmaskMultipoles;  // window multipoles set to joint field by FFTLog_setInput.

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

  for(j=0; j<mono_order; j++)  xtheory[j]                   =  convlmonoCorr->pk[fftlog_indices[j]][0];
  for(j=0; j<mono_order; j++)  xtheory[j + mono_order]      =  convlquadCorr->pk[fftlog_indices[j]][0];

  for(j=0; j<order; j++){
    ytheory[j] = 0.0;

    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++)  ytheory[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xtheory[k];
  }

  for(j=0; j<order; j++)  ChiSq += pow(ydata[j] - ytheory[j], 2.)/gsl_vector_get(eval, j);
  // for(j=0; j<order; j++)  ChiSq += pow(xdata[j] - xtheory[j], 2.)*pow(gsl_matrix_get(sigma_norm, j, j), 2.);

  return ChiSq;
}
