double C_n(double x, int n){                                                // (n+1)! = Gamma (n+2)
  return pow(HermitePolynomial(x, n-1), 2.)*exp(-2.*x*x)/(pi*pow(2., n)*gsl_sf_gamma(n + 2));
}


int FFTLog_zeroInput(FFTLog_config *fc, double beta, double velDispersion){
  double transformOrder;

  double logrmin  = log(fc->min);
  double logrmax  = log(fc->max);

  double dlogr    = (logrmax - logrmin)/(double)fc->N;
  double logrc    = (logrmax + logrmin)/2.0;

  double nc       = (double)(fc->N + 1)/2.0 -1;

  double logkc    = log(fc->kr)- logrc;

  transformOrder  =  fc->mu - 0.5;

  // write initial signal
  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]    = exp(logkc + ((double)i-nc)*dlogr);

    fc->krvals[i][1]    = exp(logrc + ((double)i-nc)*dlogr);

    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.
    fc->pk[i][0]        = 0.0;
    fc->pk[i][1]        = 0.0;

    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.
    fc->xi[i][0]        = 0.0;
    fc->xi[i][1]        = 0.0;
  }

  return 0;
}


int FFTLog_initialise(FFTLog_config *fc, double beta, double velDispersion){
  double transformOrder;

  double logrmin  = log(fc->min);
  double logrmax  = log(fc->max);

  double dlogr    = (logrmax - logrmin)/(double)fc->N;
  double logrc    = (logrmax + logrmin)/2.0;

  double nc       = (double)(fc->N + 1)/2.0 -1;

  double logkc    = log(fc->kr) - logrc;

  transformOrder  =  fc->mu - 0.5;

  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]    = exp(logkc + ((double)i-nc)*dlogr);
    fc->krvals[i][1]    = exp(logrc + ((double)i-nc)*dlogr);

    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.
    fc->pk[i][0]        = (*pt2Pk)(fc->krvals[i][0])*pow(bsigma8,  2.)*kaiserLorentz_multipole(fc->krvals[i][0]*velDispersion, beta, (int) transformOrder);
    fc->pk[i][1]        = 0.0;

    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.
    fc->xi[i][0]        = 0.0; // splint_VIPERS_maskMultipoles(fc->krvals[i][1], transformOrder); // (*pt2Xi)(fc->krvals[i][1]);
    fc->xi[i][1]        = 0.0;
  }

  return 0;
}


int FFTLog_initialise_mask(FFTLog_config *fc){
  double transformOrder;

  double logrmin  = log(fc->min);
  double logrmax  = log(fc->max);

  double dlogr    = (logrmax - logrmin)/(double)fc->N;
  double logrc    = (logrmax + logrmin)/2.0;

  double nc       = (double)(fc->N + 1)/2.0 -1;

  double logkc    = log(fc->kr) - logrc;

  transformOrder  = fc->mu - 0.5;

  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]  = exp(logkc + ((double)i-nc)*dlogr);
    fc->krvals[i][1]  = exp(logrc + ((double)i-nc)*dlogr);

    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.
    fc->pk[i][0]      = 0.0;
    fc->pk[i][1]      = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.
    fc->xi[i][0]      = splint_VIPERS_maskMultipoles(fc->krvals[i][1], (int) transformOrder);
    fc->xi[i][1]      = 0.0;
  }

  return 0;
}

int FFTLog_initialise_shot(FFTLog_config *fc){
  double transformOrder;

  double logrmin  = log(fc->min);
  double logrmax  = log(fc->max);

  double dlogr    = (logrmax - logrmin)/(double)fc->N;
  double logrc    = (logrmax + logrmin)/2.0;

  double nc       = (double)(fc->N + 1)/2.0 -1;

  double logkc    = log(fc->kr) - logrc;

  transformOrder  = fc->mu - 0.5;

  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]  = exp(logkc + ((double)i-nc)*dlogr);
    fc->krvals[i][1]  = exp(logrc + ((double)i-nc)*dlogr);

    if((0.01 < fc->krvals[i][0]) && (fc->krvals[i][0] < 1.0)){
      fc->pk[i][0]    = mean_shot;
    }
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.
    fc->pk[i][1]      = 0.0;

    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.
    fc->xi[i][0]      = 0.0;
    fc->xi[i][1]      = 0.0;
  }

  return 0;
}


int set_FFTlog(int FFTlogRes, double loval, double hival, double beta, double velDispersion){
  zero_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  shot_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  mono_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  quad_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);
  hex_config       = FFTLog_init(FFTlogRes, loval, hival, 0.0, 4 + 0.5);
  oct_config       = FFTLog_init(FFTlogRes, loval, hival, 0.0, 6 + 0.5);

  clipmono_config  = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  clipquad_config  = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);

  convlmonoCorr    = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  convlquadCorr    = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);
  convlhexCorr     = FFTLog_init(FFTlogRes, loval, hival, 0.0, 4 + 0.5);

  
  FFTLog_zeroInput(zero_config, 0.0, 0.0);

  FFTLog_initialise_shot(shot_config); // S

  FFTLog_initialise(mono_config,       beta, velDispersion);
  FFTLog_initialise(quad_config,       beta, velDispersion);
  FFTLog_initialise( hex_config,       beta, velDispersion);
  FFTLog_initialise( oct_config,       beta, velDispersion);

  FFTLog_initialise(clipmono_config,   beta, velDispersion);
  FFTLog_initialise(clipquad_config,   beta, velDispersion);

  FFTLog_initialise(convlmonoCorr,     beta, velDispersion);
  FFTLog_initialise(convlquadCorr,     beta, velDispersion);
  FFTLog_initialise(convlhexCorr,      beta, velDispersion);

  for(j=0; j<FFTlogRes; j++){
    xi_mu_prefactor[j]   = sqrt(pow(mono_config->krvals[j][0], 3.)/(8.*pow(pi, 3.))); // initialise arrays for xi_mu calc.
    xi_mu_postfactor[j]  = pow(mono_config->krvals[j][1], -1.5);

    pk_mu_prefactor[j]   = pow(2.*pi*mono_config->krvals[j][1], 3./2.);  // initialise arrays for pk_mu calc.
    pk_mu_postfactor[j]  = pow(mono_config->krvals[j][0], -1.5);
  }
  
  return 0;
}


int precompute_vipers_clipping_model(int FFTlogRes){
  // prep_kaiserLorentSpline();

  for(j=0; j<FFTlogRes; j++){    
    FFTlog_Pk[j]         = (*pt2Pk)(mono_config->krvals[j][0]);  // normalised to sigma_8 = 1.

    FFTlog_W0[j]         = splint_VIPERS_maskMultipoles(mono_config->krvals[j][1], 0); // splint W_0(r).
    FFTlog_W2[j]         = splint_VIPERS_maskMultipoles(mono_config->krvals[j][1], 2);
    FFTlog_W4[j]         = splint_VIPERS_maskMultipoles(mono_config->krvals[j][1], 4);
    FFTlog_W6[j]         = splint_VIPERS_maskMultipoles(mono_config->krvals[j][1], 6);

    FFTlog_Wk0[j]        = splint_VIPERS_kSpaceMono(mono_config->krvals[j][0]);  // splint \tilde W_0(k)
    FFTlog_Wk2[j]        = splint_VIPERS_kSpaceQuad(mono_config->krvals[j][0]);

    FFTlog_W0_joint[j]   = splint_VIPERS_jmaskMultipoles(mono_config->krvals[j][1], 0); // splint joint-field W_0(r).
    FFTlog_W2_joint[j]   = splint_VIPERS_jmaskMultipoles(mono_config->krvals[j][1], 2);
    FFTlog_W4_joint[j]   = splint_VIPERS_jmaskMultipoles(mono_config->krvals[j][1], 4);
    FFTlog_W6_joint[j]   = splint_VIPERS_jmaskMultipoles(mono_config->krvals[j][1], 6);

    // if(j % 10 == 0)  printf("\n%.6lf %.6lf %.6lf %.6lf %.6lf", FFTlog_Pk[j], FFTlog_W0[j], FFTlog_W2[j], FFTlog_Wk0[j], FFTlog_Wk2[j]);
  }

  for(i=0; i<mono_config->N;   i++){
    if((mono_config->krvals[i][1]) >= pow(10., -2.)){
      varcalc_index = i; // index where r>0.01; assume this is xi(0) for variance.

      break;
    }
  }

  for(i=0; i<mono_config->N;   i++){
    if((mono_config->krvals[i][0]) >= pow(10., -3.)){
      cnvldpk_zeropoint_index = i; // index where k>0.001; assume this is pk(0).

      break;
    }
  }

  xi_mu(shot_config); // set shot noise xi. 
  
  return 0;
}


int xi_mu(FFTLog_config* fc){
  // Last argument of FFTLog_init is the order of the Hankel transform.
  double transformOrder;

  transformOrder  =  fc->mu - 0.5;

  for(i=0; i<fc->N; i++){
    fc->input[i][0]   = xi_mu_prefactor[i]*fc->pk[i][0];
    fc->input[i][1]   = 0.0;  // Assumed zero'd.
  }

  // power law bias may be required. in particular, for Zel'dovich P(k) -> see FFT_log_zeldovich.c
  // for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][0], -fc->q);

  FFTLog(fc, fc->forwardplan, fc->backwardplan);

  // and debias.
  // for(i=0; i<fc->N; i++)  fc->output[i][0] *= pow(fc->krvals[i][1], -fc->q);

  for(i=0; i<fc->N;   i++) fc->xi[i][0]      = pow(-1., transformOrder/2)*fc->output[i][0]*xi_mu_postfactor[i];

  return 0;
}


int pk_mu(FFTLog_config* fc){
  // Last argument of FFTLog_init is the order of the Hankel transform.
  double transformOrder;

  transformOrder  =  fc->mu - 0.5;

  for(i=0; i<fc->N; i++){
    fc->input[i][0] = pk_mu_prefactor[i]*fc->xi[i][0];
    fc->input[i][1] = 0.0; // Assumes zero'd.
  }

  // power law bias may be required.
  // for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][1], -fc->q);

  FFTLog(fc, fc->forwardplan, fc->backwardplan);

  // and debias.
  // for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][0], -fc->q);

  for(i=0; i<fc->N;   i++) fc->pk[i][0] = pow(-1., transformOrder/2)*fc->output[i][0]*pk_mu_postfactor[i];

  return 0;
}


int varCalc(FFTLog_config* fc, double* sigmaSq){
  *sigmaSq = fc->xi[varcalc_index][0];

  // Assumes only the monopole contributes to the variance. Certainly quad and hex appear not to.  corr fn. seems to have converged/unaliased between [10**-3, 10**-2.]. Could 'trust'
  // anywhere in this interval.  Variance appears to asymptote as r->0. take variance as value of monopole at ~10**-2.
  // Strongly dependent on the limits over which the corr fn. is computed. Best not change them.

  return 0;
}


int clip_p0p2(FFTLog_config* clip_p0, FFTLog_config* clip_p2, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, double u0, double var){
  //#pragma omp parallel for private(i)
  for(i=0; i<clip_p0->N; i++){
    clip_p0->xi[i][0]  = clipmono_amp*mono->xi[i][0];
    clip_p0->xi[i][0] += (clip_distcoeff/var)*(pow(mono->xi[i][0], 2.) + pow(quad->xi[i][0], 2.)/5. + pow(hex->xi[i][0], 2.)/9. + pow(oct->xi[i][0], 2.)/13.);

    clip_p2->xi[i][0]  = clipmono_amp*quad->xi[i][0];
    clip_p2->xi[i][0] += (clip_distcoeff/var)*(2.*mono->xi[i][0]*quad->xi[i][0] + (2./7.)*pow(quad->xi[i][0], 2.) + (4./7.)*quad->xi[i][0]*hex->xi[i][0] + (100./693.)*pow(hex->xi[i][0], 2.) + (50./143.)*hex->xi[i][0]*oct->xi[i][0] + (14./143.)*pow(oct->xi[i][0], 2.));

    mono->xi[i][0]     = clip_p0->xi[i][0]; // what happens if no clipping.
    quad->xi[i][0]     = clip_p2->xi[i][0];
  }

  return 0;
}


int cnvldmonoCorr_joint(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* shot){
  // Convolved monopole calculation to 2nd order. updated via bailey.c. Shot noise adds to large scale power. 
  for(i=0; i<mono->N; i++){
    cnvld->xi[i][0]   =         (mono->xi[i][0] + shot->xi[i][0])*FFTlog_W0_joint[i]; // calculation to 2nd order; updated via bailey.c

    cnvld->xi[i][0]  +=  (1./5.)*quad->xi[i][0]*FFTlog_W2_joint[i];
    cnvld->xi[i][0]  +=   (1./9.)*hex->xi[i][0]*FFTlog_W4_joint[i];
  }

  return 0;
}


int cnvldmonoCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, FFTLog_config* dec, FFTLog_config* dodeca){
  // Convolved monopole calculation to 2nd order. updated via bailey.c
  for(i=0; i<mono->N; i++){
    cnvld->xi[i][0]  =          mono->xi[i][0]*FFTlog_W0[i]; // calculation to 2nd order; updated via bailey.c
    cnvld->xi[i][0] +=  (1./5.)*quad->xi[i][0]*FFTlog_W2[i];
    cnvld->xi[i][0] +=   (1./9.)*hex->xi[i][0]*FFTlog_W4[i];

    // cnvld->xi[i][0]  +=  (1./13.)*oct->xi[i][0]*(*pt2maskMultipoles)(mono->krvals[i][1], 6);
    // cnvld->xi[i][0]  +=  (1./17.)*dec->xi[i][0]*(*pt2maskMultipoles)(mono->krvals[i][1], 8);
    // cnvld->xi[i][0]  +=  (1./21.)*dodeca->xi[i][0]*(*pt2maskMultipoles)(mono->krvals[i][1], 10);
  }

  return 0;
}


int cnvldquadCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, FFTLog_config* dec, FFTLog_config* dodeca){
  // Convolved quadrupole calculation to 2nd order. updated via bailey.c
  for(i=0; i<mono->N; i++){
    cnvld->xi[i][0]  = mono->xi[i][0]*FFTlog_W2[i];
    cnvld->xi[i][0] += quad->xi[i][0]*(FFTlog_W0[i] + (2./7.)*FFTlog_W2[i] + (2./7.)*FFTlog_W4[i]);
    cnvld->xi[i][0] +=  hex->xi[i][0]*((2./7.)*FFTlog_W2[i] + (100./693.)*FFTlog_W4[i] + (25./143.)*FFTlog_W6[i]);
    
    // cnvld->xi[i][0] += oct->xi[i][0]*((25./143.)*(*pt2maskMultipoles)(mono->krvals[i][1], 4)+(14./143.)*(*pt2maskMultipoles)(mono->krvals[i][1], 6)+(28./221.)*(*pt2maskMultipoles)(mono->krvals[i][1], 8));

// cnvld->xi[i][0]+= dec->xi[i][0]*((28./221.)*(*pt2maskMultipoles)(mono->krvals[i][1], 6)+(24./323.)*(*pt2maskMultipoles)(mono->krvals[i][1], 8)+(225./2261.)*(*pt2maskMultipoles)(mono->krvals[i][1], 10));
}

return 0;
}

int addshotnoise(FFTLog_config *mono, double shotlevel){
  for(i=0; i<mono->N; i++){
    if((0.01 < mono_config->krvals[i][0]) && (mono_config->krvals[i][0] < 1.0)){
      mono->pk[i][0] += shotlevel;
    }
  }
  
  return 0;
}

int FFTlog_updatepk(FFTLog_config *mono, FFTLog_config *quad, FFTLog_config *hex, double beta, double velDispersion){
  int klo = 0; // kaiser Lorbsigma8entz D^2 assumes 1000 element array.

  double                           ks;
  double mu_0, mu_2, mu_4, mu_6, mu_8;
  
  for(i=0; i<mono->N; i++){
    ks   = mono->krvals[i][0]*velDispersion;
    
    mu_0 = muOrderZero(ks);  
    mu_2 = muOrderTwo(ks);   
    mu_4 = muOrderFour(ks);
    mu_6 = muOrderSix(ks);   
    mu_8 = muOrderEight(ks); 
    /*
    mu_0 = seqsp_kLmu(ks, 0, &klo);
    mu_2 = seqsp_kLmu(ks, 2, &klo);
    mu_4 = seqsp_kLmu(ks, 4, &klo);
    mu_6 = seqsp_kLmu(ks, 6, &klo);
    mu_8 = seqsp_kLmu(ks, 8, &klo);  
    */
    mono->pk[i][0] = FFTlog_Pk[i]*pow(bsigma8, 2.)*(mu_0 + 2.*beta*mu_2 + beta*beta*mu_4);
    quad->pk[i][0] = FFTlog_Pk[i]*pow(bsigma8, 2.)*((5./2.)*(-mu_0 + mu_2*(3. - 2.*beta) + mu_4*(-beta*beta + 6.*beta) + 3.*beta*beta*mu_6));
     hex->pk[i][0] = FFTlog_Pk[i]*pow(bsigma8, 2.)*((9./8.)*(35.*beta*beta*mu_8  + 10.*beta*(7. -3.*beta)*mu_6 + (35. - 60.*beta + 3.*beta*beta)*mu_4 + 6.*(beta - 5.)*mu_2 + 3.*mu_0));

    mono->pk[i][1] = 0.0;
    quad->pk[i][1] = 0.0;
     hex->pk[i][1] = 0.0;
  }

  return 0;
}


int FFTlog_checkupdatepk(FFTLog_config *mono, FFTLog_config *quad, FFTLog_config *hex, double beta, double velDispersion){
  int klo = 0; // kaiser Lorentz D^2 assumes 1000 element array.

  double   ks;
  double mu_0, mu_2, mu_4, mu_6, mu_8;

  sprintf(filepath, "%s/W1_Spectro_V7_4/kL_factors.dat", root_dir);

  output = fopen(filepath, "w");

  for(i=0; i<mono->N; i++){
    ks   = mono->krvals[i][0]*velDispersion;

    mu_0 = seqsp_kLmu(ks, 0, &klo);
    mu_2 = seqsp_kLmu(ks, 2, &klo);
    mu_4 = seqsp_kLmu(ks, 4, &klo);
    mu_6 = seqsp_kLmu(ks, 6, &klo);
    mu_8 = seqsp_kLmu(ks, 8, &klo);

    fprintf(output, "%.4le \t %.4le \t", kaiserLorentz_multipole(ks, beta, 0), (mu_0 + 2.*beta*mu_2 + beta*beta*mu_4));
    fprintf(output, "%.4le \t %.4le \t", kaiserLorentz_multipole(ks, beta, 2), (5./2.)*( -1.*mu_0 + mu_2*(3. - 2.*beta) + mu_4*(-beta*beta + 6.*beta) + 3.*beta*beta*mu_6));
    fprintf(output, "%.4le \t %.4le \n", kaiserLorentz_multipole(ks, beta, 4), ((9./8.)*(35.*beta*beta*mu_8  + 10.*beta*(7. -3.*beta)*mu_6 + (35. - 60.*beta + 3.*beta*beta)*mu_4 + 6.*(beta-5.)*mu_2 + 3.*mu_0)));
  }

  fclose(output);

  return 0;
}
