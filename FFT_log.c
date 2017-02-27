double C_n(double x, int n){                                                // (n+1)! = Gamma (n+2)
  return pow(HermitePolynomial(x, n-1), 2.)*exp(-2.*x*x)/(pi*pow(2., n)*gsl_sf_gamma(n + 2));
}


int FFTlog_memory(int FFTlogRes, double loval, double hival, double beta, double velDispersion){
    zero_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
    mono_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
    quad_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);
     hex_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 4 + 0.5);
     oct_config      = FFTLog_init(FFTlogRes, loval, hival, 0.0, 6 + 0.5);
                   
    clipmono_config  = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);  
    clipquad_config  = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);  
    
    convlmonoCorr    = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);  
    convlquadCorr    = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);  
    convlhexCorr     = FFTLog_init(FFTlogRes, loval, hival, 0.0, 4 + 0.5);  
        
    FFTLog_zeroInput(zero_config, 0.0, 0.0);

    FFTLog_initialise( mono_config,      beta, velDispersion); 
    FFTLog_initialise( quad_config,      beta, velDispersion);
    FFTLog_initialise(  hex_config,      beta, velDispersion);
    FFTLog_initialise(  oct_config,      beta, velDispersion);
    
    FFTLog_initialise(clipmono_config,   beta, velDispersion);
    FFTLog_initialise(clipquad_config,   beta, velDispersion);
    
    FFTLog_initialise(convlmonoCorr,     beta, velDispersion);
    FFTLog_initialise(convlquadCorr,     beta, velDispersion);
    FFTLog_initialise(convlhexCorr,      beta, velDispersion);

    prep_kaiserLorentSpline();
    
    for(j=0; j<FFTlogRes; j++){
      FFTlog_Pk[j]         = (*pt2Pk)(mono_config->krvals[j][0])*pow(app_sigma8, -2.);  // normalise to sigma_8 = 1. 
      
      xi_mu_prefactor[j]   = sqrt(pow(mono_config->krvals[j][0], 3.)/(8.*pow(pi, 3.))); // initialise arrays for xi_mu calc.  
      xi_mu_postfactor[j]  = pow(mono_config->krvals[j][1], -1.5);

      pk_mu_prefactor[j]   = pow(2.*pi*mono_config->krvals[j][1], 3./2.);  // initialise arrays for pk_mu calc.
      pk_mu_postfactor[j]  = pow(mono_config->krvals[j][0], -1.5);

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
    
    return 0; 
}


int FFTlog_lnnormmemory(int FFTlogRes, double loval, double hival, double beta, double velDispersion){
  lnnorm_mono_config   = FFTLog_init(FFTlogRes, loval, hival, 0.0, 0 + 0.5);
  lnnorm_quad_config   = FFTLog_init(FFTlogRes, loval, hival, 0.0, 2 + 0.5);
  
  FFTLog_initialise(lnnorm_mono_config,    beta, velDispersion);
  FFTLog_initialise(lnnorm_quad_config,    beta, velDispersion);
  
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
    fc->pk[i][0]        = 
                          //** MUST TURN ON NO RSD, OTHERWISE QUAD is SET TO MONO ETC. **
                          (*pt2Pk)(fc->krvals[i][0])
                          
                          *pow(app_sigma8, -2.)
                          *pow(   bsigma8,  2.) 
                          *kaiserLorentz_multipole(fc->krvals[i][0]*velDispersion, beta, (int) transformOrder);
                    
    fc->pk[i][1]        = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
    fc->xi[i][0]        = splint_VIPERS_maskMultipoles(fc->krvals[i][1], transformOrder); // (*pt2Xi)(fc->krvals[i][1]);
    fc->xi[i][1]        = 0.0;
  }
  
  return 0;
}


int FFTlog_updatepk(FFTLog_config *mono, FFTLog_config *quad, FFTLog_config *hex, double beta, double velDispersion){
  int klo = 0; // kaiser Lorentz D^2 assumes 1000 element array.

  double   ks;
  double mu_0, mu_2, mu_4, mu_6, mu_8;
    
  for(i=0; i<mono->N; i++){
    ks   = mono->krvals[i][0]*velDispersion;

    mu_0 = seqsp_kLmu(ks, 0, &klo);
    mu_2 = seqsp_kLmu(ks, 2, &klo);
    mu_4 = seqsp_kLmu(ks, 4, &klo);
    mu_6 = seqsp_kLmu(ks, 6, &klo);
    mu_8 = seqsp_kLmu(ks, 8, &klo);
    
    mono->pk[i][0] = FFTlog_Pk[i]*pow(bsigma8, 2.)*(mu_0 + 2.*beta*mu_2 + beta*beta*mu_4); //  *kaiserLorentz_multipole(ks, beta, 0)
    quad->pk[i][0] = FFTlog_Pk[i]*pow(bsigma8, 2.)*((5./2.)*(-1.*mu_0 + mu_2*(3. - 2.*beta) + mu_4*(-beta*beta + 6.*beta) + 3.*beta*beta*mu_6));
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
    fprintf(output, "%.4le \t %.4le \n", kaiserLorentz_multipole(ks, beta, 4), ((9./8.)*(35.*beta*beta*mu_8  + 10.*beta*(7. -3.*beta)*mu_6 + (35. - 60.*beta + 3.*beta*beta)*mu_4 + 6.*(beta - 5.)*mu_2 + 3.*mu_0)));
  }

  fclose(output);
   
  return 0;
}


/*
int FastLegendre_setInput(FFTLog_config *fc_mono, FFTLog_config *fc_quad, FFTLog_config *fc_hex, FFTLog_config *fc_oct, double beta, double velDispersion){  
  parameters[2]   =          beta;
  parameters[3]   = velDispersion;
  
  double kmod;
  
  // write initial signal
  for(ii=0; ii<fc_mono->N; ii++){
    // kmod sent to FastLegendre function.c
    kmod                     = fc_mono->krvals[ii][0];
    parameters[4]            = kmod;
    
    if(ggb_coeff_fft(fLegendre_func, ggb, 256, 0.5, 5000, 2, 1e-12) == 1)  error (" Error computing the Gegenbauer coefficients");
    
    // printf("\n%d \t %e \t %e", ii, fc_mono->krvals[ii][0], fc_mono->pk[ii][0]);
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.     
    fc_mono->pk[ii][0]       = (*pt2Pk)(fc_mono->krvals[ii][0])*ggb[0];
    fc_mono->pk[ii][1]       = 0.0;
    
    // printf("\n%d \t %e \t %e \t %e \n", ij, fc_mono->krvals[ij][0], fc_mono->pk[ij][0], mono_config->pk[ij][0]);
    
    fc_quad->pk[ii][0]       = (*pt2Pk)(fc_quad->krvals[ii][0])*ggb[2];
    fc_quad->pk[ii][1]       = 0.0;
    
    fc_hex->pk[ii][0]        = (*pt2Pk)( fc_hex->krvals[ii][0])*ggb[4];
    fc_hex->pk[ii][1]        = 0.0;
    
    fc_oct->pk[ii][0]        = (*pt2Pk)( fc_oct->krvals[ii][0])*ggb[6];
    fc_oct->pk[ii][1]        = 0.0;
  }
  
  return 0;
}
*/

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
  clip_p0->xi[i][0] = clipmono_amp*mono->xi[i][0];
  clip_p0->xi[i][0]+= (clip_distcoeff/var)*(pow(mono->xi[i][0], 2.) + pow(quad->xi[i][0], 2.)/5. + pow(hex->xi[i][0], 2.)/9. + pow(oct->xi[i][0], 2.)/13.);

  clip_p2->xi[i][0] = clipmono_amp*quad->xi[i][0];
  clip_p2->xi[i][0]+= (clip_distcoeff/var)*(2.*mono->xi[i][0]*quad->xi[i][0] + (2./7.)*pow(quad->xi[i][0], 2.) + (4./7.)*quad->xi[i][0]*hex->xi[i][0] + (100./693.)*pow(hex->xi[i][0], 2.) + (50./143.)*hex->xi[i][0]*oct->xi[i][0] + (14./143.)*pow(oct->xi[i][0], 2.));

  mono->xi[i][0] = clip_p0->xi[i][0]; // what happens if no clipping. 
  quad->xi[i][0] = clip_p2->xi[i][0];
 }

 return 0;
}


int clipmono(FFTLog_config* clip, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, double u0, double var){
  // Clipped monopole calculation to 2nd order.
  for(i=0; i<clip->N; i++){
    clip->xi[i][0]   = clipmono_amp*mono->xi[i][0];


    clip->xi[i][0]  += (C_n(u0, 1)/var)*(pow(mono->xi[i][0], 2.) + pow(quad->xi[i][0], 2.)/5. + pow(hex->xi[i][0], 2.)/9. + pow(oct->xi[i][0], 2.)/13.);
  }

  return 0;
}


int clipquad(FFTLog_config* clip, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, double u0, double var){
  // Clipped quadrupole calculation to 2nd order.
  for(i=0; i<clip->N; i++){
    clip->xi[i][0]   = 0.25*pow(1. + gsl_sf_erf(u0), 2.)*quad->xi[i][0];

    clip->xi[i][0]  += (C_n(u0, 1)/var)*(2.*mono->xi[i][0]*quad->xi[i][0] + (2./7.)*pow(quad->xi[i][0], 2.) + (4./7.)*quad->xi[i][0]*hex->xi[i][0] + (100./693.)*pow(hex->xi[i][0], 2.) + (50./143.)*hex->xi[i][0]*oct->xi[i][0] + (14./143.)*pow(oct->xi[i][0], 2.));
  }

  return 0;
}


int cliphex(FFTLog_config* clip, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, double u0, double var){
  // Clipped hexadecapole calculation to 2nd order.                                                                                                                                                                                         
  for(i=0; i<clip->N; i++){
    clip->xi[i][0]   = 0.25*pow(1. + gsl_sf_erf(u0), 2.)*hex->xi[i][0];

    clip->xi[i][0]  += (C_n(u0, 1)/var)*( (18./35.)*pow(quad->xi[i][0], 2.) + 2.*mono->xi[i][0]*hex->xi[i][0] + (40./77.)*quad->xi[i][0]*hex->xi[i][0] + (162./1001.)*pow(hex->xi[i][0], 2.) + (90./143.)*quad->xi[i][0]*oct->xi[i][0]                                                                                 + (40./143.)*hex->xi[i][0]*oct->xi[i][0] +(252./2431.)*pow(oct->xi[i][0], 2.));
  }

  return 0;
}


int lnnorm(FFTLog_config* lnnorm, FFTLog_config* mono){
    // Correlation fn. in the log-normal model, when the Gaussian field is **isotropic**   
    for(i=0; i<lnnorm->N; i++){  
        if((mono->krvals[i][1] > 0.0001) && (mono->krvals[i][1] < 10000.)){
            lnnorm->xi[i][0]   = exp(mono->xi[i][0]) - 1.;
        }
    }
    
    return 0;
}


int invert_lnnorm(FFTLog_config* mono, FFTLog_config* mono2){
    // inversion of the correlation fn. in the log-normal model, when the Gaussian field is **isotropic**   
    for(i=0; i<mono->N; i++){  
        if((mono->krvals[i][1] > 0.0001) && (mono->krvals[i][1] < 10000.)){ 
              mono->xi[i][0]   = log(1. + mono2->xi[i][0]);
        }
    }
    
    return 0;
}

/*
int clip_lnnorm(FFTLog_config* clip_lnnorm, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, double u0, double var, int order){
    double sigma = sqrt(var);
    
    double Jn[15];
    
    printf("\n\nclipped ln norm, sigma: %.3lf, u0: %.3lf", sigma, u0);
    
    //J_0
    Jn[0]  =               J0_calc(u0, sigma);
    
    Jn[1]  = recursion(Jn[0],   1, u0, sigma);
    Jn[2]  = recursion(Jn[1],   2, u0, sigma);
    Jn[3]  = recursion(Jn[2],   3, u0, sigma);
    Jn[4]  = recursion(Jn[3],   4, u0, sigma);
    Jn[5]  = recursion(Jn[4],   5, u0, sigma);
    Jn[6]  = recursion(Jn[5],   6, u0, sigma);
    Jn[7]  = recursion(Jn[6],   7, u0, sigma);
    Jn[8]  = recursion(Jn[7],   8, u0, sigma);
    Jn[9]  = recursion(Jn[8],   9, u0, sigma);
    Jn[10] = recursion(Jn[9],  10, u0, sigma);
    Jn[11] = recursion(Jn[10], 11, u0, sigma);
    Jn[12] = recursion(Jn[11], 12, u0, sigma);
    Jn[13] = recursion(Jn[12], 13, u0, sigma);
    Jn[14] = recursion(Jn[13], 14, u0, sigma);
            
    printf("\n\nJn*Jn/(n+1)!: %.5lf \t %.5lf \t %.5lf \t %.5lf \t %.5lf", Jn[0]*Jn[0]/gsl_sf_fact(1), Jn[1]*Jn[1]/gsl_sf_fact(2), Jn[2]*Jn[2]/gsl_sf_fact(3), 
                                                                          Jn[3]*Jn[3]/gsl_sf_fact(4), Jn[4]*Jn[4]/gsl_sf_fact(5));
            
    // Clipped ln norm. assumes isotropic correlation fn. 
    for(i=0; i<clip_lnnorm->N; i++){
      clip_lnnorm->xi[i][0]   = mono->xi[i][0]; 
      
      if((mono->krvals[i][1] > 0.001) && (mono->krvals[i][1] < 2000.)){
        clip_lnnorm->xi[i][0] = 0.0;
      
        for(k=0; k<order; k++)  clip_lnnorm->xi[i][0]  += series_term(mono->xi[i][0], Jn[k], k);
      }
    }
    
    // printf("\n\nu0: %.3lf, 2.*J_0^2/J_1^2: %.3lf", u0, 2.*Jn[0]*Jn[0]/Jn[1]/Jn[1]);

    return 0;
}
*/

int add_fogRSD(FFTLog_config *multipole, int Order){
    for(j=0; j<multipole->N; j++){
        if((multipole->krvals[j][0] > 0.0001) && (multipole->krvals[j][0] < 2000.)){
            multipole->pk[j][0] *= 1.; //Lorentz_multipole(multipole->krvals[j][0]*velDispersion, 0);
        }
    }

    return 0;
}


int cnvldmonoCorr_joint(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex){
  // Convolved monopole calculation to 2nd order. updated via bailey.c
  for(i=0; i<mono->N; i++){
    cnvld->xi[i][0]   =          mono->xi[i][0]*FFTlog_W0_joint[i]; // calculation to 2nd order; updated via bailey.c
    cnvld->xi[i][0]  +=  (1./5.)*quad->xi[i][0]*FFTlog_W2_joint[i];
    cnvld->xi[i][0]  +=  (1./9.)*hex->xi[i][0]*FFTlog_W4_joint[i];
  }

  return 0;
}


int cnvldmonoCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, FFTLog_config* oct, FFTLog_config* dec, FFTLog_config* dodeca){
 // Convolved monopole calculation to 2nd order. updated via bailey.c   
 for(i=0; i<mono->N; i++){
   cnvld->xi[i][0]   =          mono->xi[i][0]*FFTlog_W0[i]; // calculation to 2nd order; updated via bailey.c    
      
   cnvld->xi[i][0]  +=  (1./5.)*quad->xi[i][0]*FFTlog_W2[i];       
    
   cnvld->xi[i][0]  +=  (1./9.)*hex->xi[i][0]*FFTlog_W4[i];
    
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
      
  cnvld->xi[i][0] += quad->xi[i][0]*(FFTlog_W0[i]+(2./7.)*FFTlog_W2[i]+(2./7.)*FFTlog_W4[i]);  
      
  cnvld->xi[i][0] += hex->xi[i][0]*((2./7.)*FFTlog_W2[i]+(100./693.)*FFTlog_W4[i]+(25./143.)*FFTlog_W6[i]);
    
  // cnvld->xi[i][0] += oct->xi[i][0]*((25./143.)*(*pt2maskMultipoles)(mono->krvals[i][1], 4)+(14./143.)*(*pt2maskMultipoles)(mono->krvals[i][1], 6)+(28./221.)*(*pt2maskMultipoles)(mono->krvals[i][1], 8));
    
  // cnvld->xi[i][0]+= dec->xi[i][0]*((28./221.)*(*pt2maskMultipoles)(mono->krvals[i][1], 6)+(24./323.)*(*pt2maskMultipoles)(mono->krvals[i][1], 8)+(225./2261.)*(*pt2maskMultipoles)(mono->krvals[i][1], 10));
 }
    
 return 0;
}


/*
int clipped_lnnorm_pkCalc(){
    FFTlogRes     =               4096;
    
    fsigma8       =                0.0;  
    bsigma8       =        0.557*1.495;
    velDispersion =               0.00;
    
    
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), fsigma8/bsigma8, velDispersion);
     
    xi_mu(mono_config);
    
    // non-linear halofit p(k) to p(k) of the underlying Gaussian field in the lognormal model.  
    invert_lnnorm(mono_config, mono_config);

    // variance of the underlying Gaussian. 
    varCalc(mono_config, &variance);
    
    
    double f = 0.66;
    
    // u0 calc. in the absence of the value of \delta_0.
    u0 = inverse_erf(2.*f - 1.)  + sqrt(variance)/sqrt(2.);
      
    printf("\n\nvariance: %.3lf, \t u0: %.3lf", variance, u0);
        
    // lnnorm(lnnorm_mono_config, mono_config);    

    // correlation fn. in the clipped ln normal model. 
    clip_lnnorm(clipmono_config, mono_config, quad_config, zero_config, u0, variance, 14);
    
    
    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/fftlog_multipoles_clipped_lnnormal_hodcube_%.2lf.dat", root_dir, f);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){
      if((mono_config->krvals[j][1] > 1.) && (mono_config->krvals[j][1] < 250.)){
        fprintf(output, "%e \t %e \t %e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0], lnnorm_mono_config->xi[j][0], clipmono_config->xi[j][0]);
      }
    }
    
    fclose(output); 
    
    
    pk_mu(mono_config);
    
    // pk_mu(lnnorm_mono_config);
    
    pk_mu(clipmono_config);
    
    
    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/fftlog_multipoles_clipped_lnnormal_hodcube_f_%.2lf_pk.dat", root_dir, f);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){
      if((mono_config->krvals[j][0] > 0.01) && (mono_config->krvals[j][0] < 1.)){
        fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0], clipmono_config->pk[j][0]);
      }
    }
    
    fclose(output); 
    
    return 0;
}
*/

int fog_xiCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), beta, velDispersion);
    
    xi_mu(mono_config);
    
    xi_mu(quad_config);
    
    // pk_mu(lnnorm_mono_config);
    
    // xi_mu(lnnorm_mono_config);
    
    
    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_clustered_fog_500_ximultipoles/fftlog_multipoles_xi_Gaussian_fogRSD.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][1] > 5.) && (mono_config->krvals[j][1] < 300.)){
            fprintf(output, "%e \t %e \t %e\n", mono_config->krvals[j][1], mono_config->xi[j][0], quad_config->xi[j][0]);
        }
    }
    
    fclose(output); 

    return 0;
}

/*
int spherical_randDistribution(){
    int cnvldpknorm;

    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    pt2Pk = &spherical_tophat_pk;
    
    
    FFTLog_setInput(W2_config, beta, velDispersion);
    
    xi_mu(W2_config);
    
    
    pt2Pk = &linearPk_Gaussian;
    
    FFTLog_setInput(mono_config,      beta, velDispersion);
    
    xi_mu(mono_config);
    
    
    for(j=0; j<mono_config->N; j++)   convlmonoCorr->xi[j][0] = mono_config->xi[j][0]*W2_config->xi[j][0];
    
    
    pk_mu(mono_config);
    
    pk_mu(convlmonoCorr);
    
    
    cnvldpknorm = cnvldpk_norm(mono_config);
    
    printf("\n %e", mono_config->pk[cnvldpknorm][0]);
    
    
    for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] *= mono_config->pk[cnvldpknorm][0]/convlmonoCorr->pk[cnvldpknorm][0];
    
    
    sprintf(filepath, "%s/Data/stacpolly/spherical_RandDistribution_monoxi.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<W2_config->N; j++){  
        if((W2_config->krvals[j][1] > 1.) && (W2_config->krvals[j][1] < 1000.)){
            fprintf(output, "%e \t %e \n", W2_config->krvals[j][1], W2_config->xi[j][0]);
        }
    }
    
    fclose(output); 
    
    
    sprintf(filepath, "%s/Data/stacpolly/spherical_RandDistribution_cnvldpk.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 10.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0], convlmonoCorr->pk[j][0]);
        }
    }
    
    fclose(output); 

    return 0; 
}
*/

int VIPERS_mask_clippedpk(){
    FFTlogRes     =     4096;
    // best fit unclipped, kmax 0.7. fsig8 set to truth
    // fsigma8       =    0.491;  
    // bsigma8       =    0.825;
    // velDispersion =      5.6;
    
    fsigma8       =        0.35;  
    bsigma8       =       0.825;
    velDispersion =         3.0;
    
    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);
        
    FFTLog_initialise(mono_config, fsigma8/bsigma8, velDispersion);
    FFTLog_initialise(quad_config, fsigma8/bsigma8, velDispersion);
    

    // Transform to correlation functions, mono and quad. 
    xi_mu(mono_config);
    xi_mu(quad_config);
    // xi_mu( hex_config);
    
    varCalc(mono_config, &variance);
    
    // u0 calc. in the absence of the value of \delta_0.
    u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);
    // u0  = appliedClippingThreshold/sqrt(2.*variance);
    
    printf("\nu0: %e, variance: %e", u0, variance);
        
    // cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    // cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    
    // Currently evaluated at second order. 
    clipmono(clipmono_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
    clipquad(clipquad_config, mono_config, quad_config, zero_config, zero_config, u0, variance);
    
    // could include hex terms of correlation fn. 
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);

    
    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/cnvldpk/clipped_d0_6_BA_cnvldpk_W1_Nagoya_v6_Samhain_incmock_specweight_nbar_xi_0.6_0.9_hex_linear_matterpk.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 10.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
        }
    }
    
    fclose(output); 
    
    return 0; 
}


int VIPERS_mask_cnvldpk(){
    FFTlogRes     =       4096;

    // z_eff: 0.75.  sigma_8(z_eff): 0.5930. f(z_eff): 0.8275. f(z_eff)*sigma_8(z_eff): 0.4907.
    // z_eff: 1.05.  sigma_8(z_eff): 0.5180. f(z_eff): 0.8823. f(z_eff)*sigma_8(z_eff): 0.4570.

    if(lo_zlim < 0.8){  
      fsigma8       =     0.4907;  
      // fsigma8       =       0.25; 
      bsigma8       =       0.70;
      velDispersion =        5.5;
    }

    if(lo_zlim > 0.8){
      fsigma8       =        0.4570;
      bsigma8       =        0.8;
      velDispersion =        4.0;
    }

    pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), fsigma8/bsigma8, velDispersion);
    
    // printf("\n\nStarting fft log clock.");
    // clock_t start = clock(), diff;
    // for(ii=0; ii<1000; ii++){
    
    double fieldArea    = 0.0;
    double cnvldpk_zero = 0.0;

    pt2maskMultipoles = &splint_VIPERS_jmaskMultipoles;

    FFTLog_initialise(mono_config, fsigma8/bsigma8, velDispersion);
    FFTLog_initialise(quad_config, fsigma8/bsigma8, velDispersion);
    FFTLog_initialise( hex_config, fsigma8/bsigma8, velDispersion);

    // Transform to correlation functions, mono and quad.                                                                                                                                                         
    xi_mu(mono_config);
    xi_mu(quad_config);
    xi_mu( hex_config);

    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, zero_config, zero_config, zero_config);

    pk_mu(convlmonoCorr);

    // Zero point determined for P'(0) for the joint fields.                                                                                                                                                      
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
    
    // VIPERS window systematics plot. 
    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = mono_config->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = quad_config->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);
    
    // diff = clock() - start;
    // int msec = diff*1000/CLOCKS_PER_SEC;
    // printf("\n\nFFT log approach to convolution, for 1000 models: %d seconds %d milliseconds", msec/1000, msec%1000);
    
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/cnvldpk/cnvldpk_W1_Nagoya_v6_Samhain_incmock_specweight_nbar_xi_0.6_0.9_hex_linearRSD.dat", root_dir);
    // sprintf(filepath, "%s/W1_Spectro_V7_1/cnvldpk_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_xi_0.6_0.9_hex.dat", root_dir, fieldFlag);
    sprintf(filepath, "%s/W1_Spectro_V7_2/cnvldpk_models/cnvldpk_nonlinear_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_xi_%.1lf_%.1lf_intcor_hex_fsig8_%.2lf.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, fsigma8);

    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 10.)){
            fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0], mono_config->pk[j][0], quad_config->pk[j][0]);
        }
    }
    
    fclose(output); 
    
    return 0; 
}


int printf_kMask_multipoles(){
    //  Hankel transform pair counts of the window. 
    int cnvldpknorm;

    FFTlogRes = 4096;
    
    pt2maskMultipoles = &splint_VIPERS_maskMultipoles;
    
    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), beta, velDispersion);
    
    // correlation_fns already assigned to the window in FFTlog_memory. 
    pk_mu(mono_config);    
    pk_mu(quad_config);
        
    double norm = 0.0;
    // norm = 4.700981*1.823239*pow(10., 6.);
        
    for(j=0; j<mono_config->N; j++){
        if(mono_config->krvals[j][0] > 0.0001){
            norm = mono_config->pk[j][0];
        
            break;
        }
    }
       
    
    sprintf(filepath, "%s/W1_Spectro_V7_4/Qmultipoles/mask_kmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf.dat", root_dir, fieldFlag, 
                                                                                                                                                  lo_zlim, hi_zlim);    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 100.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0]/norm, quad_config->pk[j][0]/norm);
        }
    }
    
    fclose(output);
    
    return 0; 
}


int Edinburgh_PDRA(){
    FFTlogRes     =      4068; 
    
    beta          =       0.5;
    velDispersion =       0.0000001;

    FFTlog_memory(FFTlogRes, pow(10., -6.), pow(10., 6.), beta, velDispersion);

    xi_mu(mono_config);
    xi_mu(quad_config);


    sprintf(filepath, "%s/EdinburghPDRA_ximultipoles_beta_%.2lf.dat", root_dir, beta);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][1] > 1.00) && (mono_config->krvals[j][1] < 100.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0], quad_config->xi[j][0]);
        }
    }
    
    fclose(output); 

    return 0;
}

/*
int kaiserLorentz_convergence(){
    // low, but sufficient resolution. 
    // FFTlogRes  =       768;
    FFTlogRes     =      4068; 
    
    beta          =       0.5;
    velDispersion =       5.0;
    
    // pt2maskMultipoles = &splint_unit_maskMultipoles;
    pt2maskMultipoles = &splint_VIPERS_maskMultipoles;
    
    // Currently beta is set by hand in the input P(k) function. 10**-4. to 10.**4. range is sufficient. 
    FFTlog_memory(FFTlogRes, pow(10., -6.), pow(10., 6.), beta, velDispersion);
    
    FastLegendre_setInput(mono_config, quad_config, hex_config, oct_config, 0.5, 2.0);
    
    // correlation_fns already assigned to the window in FFTlog_memory. 
    xi_mu(mono_config);
    xi_mu(quad_config);
    xi_mu( hex_config);
    xi_mu( oct_config);
        
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config, oct_config, zero_config, zero_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config, oct_config, zero_config, zero_config);

    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    
    sprintf(filepath, "%s/Scripts/FastLegendre/kaiserLorentz_vipers_window_cnvldpk_beta_0.5_sigma_5.0_inc_intcnstrnt.dat", root_dir, velDispersion);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){      
      if((mono_config->krvals[j][0] > 0.01) && (mono_config->krvals[j][0] < 1.)){      
        fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \n",mono_config->krvals[j][0],mono_config->pk[j][0],quad_config->pk[j][0],hex_config->pk[j][0], (*pt2Pk)(mono_config->krvals[j][0])*kaiserLorentz_multipole(mono_config->krvals[j][0]*velDispersion, 0.5,0), (*pt2Pk)(mono_config->krvals[j][0])*kaiserLorentz_multipole(mono_config->krvals[j][0]*velDispersion, 0.5,2), (*pt2Pk)(mono_config->krvals[j][0])*kaiserLorentz_multipole(mono_config->krvals[j][0]*velDispersion, 0.5,4));
        // fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0], quad_config->xi[j][0], hex_config->xi[j][0], oct_config->xi[j][0]);
      }
    }
    
    fclose(output);
    
    return 0; 
}
*/


/* FFTLog */                                                        
void FFTLog(FFTLog_config *fc, fftw_plan p_forward, fftw_plan p_backward){
  /* Hamilton 2000. http://casa.colorado.edu/~ajsh/FFTLog/

     The FFTLog algorithm for taking the discrete Hankel transform, equation (22), 
     of a sequence an of N logarithmically spaced points is:

     * FFT an to obtain the Fourier coefficients cm, equation (13);
     * multiply by um given by equations (18) and (19) to obtain cm um;
     * FFT cm um back to obtain the discrete Hankel transform ãn, equation (21). 
  */
  
  //cm's: FFT forward
  fftw_execute(p_forward);
    
  // result of FFT placed in fc->cm, as is clear from the FFTw plan. 
  // um*cm, real and imaginary parts. 
  fc->cmum[0][0] = fc->cm[0][0]*fc->um[0][0] - fc->cm[0][1]*fc->um[0][1];
  fc->cmum[0][1] = fc->cm[0][0]*fc->um[0][1] + fc->cm[0][1]*fc->um[0][0];

  for(i=1; i<fc->N/2+1; i++){
    fc->cmum[i][0] = fc->cm[i][0]*fc->um[i][0] - fc->cm[i][1]*fc->um[i][1];
    fc->cmum[i][1] = fc->cm[i][0]*fc->um[i][1] + fc->cm[i][1]*fc->um[i][0];
    
    //Hermitian symetry (i.e. to get a real signal after FFT back)
    fc->cmum[fc->N-i][0] =  fc->cmum[i][0];
    fc->cmum[fc->N-i][1] = -fc->cmum[i][1];
  }
  
  /*
  // appropriate for the inverse transform, with (kr)^q biasing. 
  fc->cmum[0][0] = (fc->cm[0][0]*fc->um[0][0] + fc->cm[0][1]*fc->um[0][1])*pow(fc->um[0][0]*fc->um[0][0] + fc->um[0][1]*fc->um[0][1], -1.);
  fc->cmum[0][1] = (fc->cm[0][1]*fc->um[0][0] - fc->cm[0][0]*fc->um[0][1])*pow(fc->um[0][0]*fc->um[0][0] + fc->um[0][1]*fc->um[0][1], -1.);

  for(i=1; i<fc->N/2+1; i++){
    fc->cmum[i][0] = (fc->cm[i][0]*fc->um[i][0] + fc->cm[i][1]*fc->um[i][1])*pow(fc->um[i][0]*fc->um[i][0] + fc->um[i][1]*fc->um[i][1], -1.);
    fc->cmum[i][1] = (fc->cm[i][1]*fc->um[i][0] - fc->cm[i][0]*fc->um[i][1])*pow(fc->um[i][0]*fc->um[i][0] + fc->um[i][1]*fc->um[i][1], -1.);
    
    //Hermitian symetry (i.e. to get a real signal after FFT back)
    fc->cmum[fc->N-i][0] =  fc->cmum[i][0];
    fc->cmum[fc->N-i][1] = -fc->cmum[i][1];
  }
  */
  
  //xi's: FFT backward
  fftw_execute(p_backward);
  
  for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->output[i][0], fc->output[fc->N-i-1][0]);
  
  for(i=0; i<fc->N;   i++) fc->output[i][0] = (fc->output[i][0]/(double)fc->N);
  
  return;
}


FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu){
  /*Initializes what FFTLog needs.*/
  
  FFTLog_config *fc = (FFTLog_config*)  malloc(sizeof(FFTLog_config));

  //FFTW3 Initialization
  fc->min        = min;
  fc->max        = max;
  fc->q          = q;
  fc->mu         = mu;
  fc->N          = N;
  
  fc->pk         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fc->xi         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  
  fc->cm         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fc->um         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fc->cmum       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
 
  fc->input      = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fc->output     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  
  fc->krvals     = malloc(sizeof(double*)*N);
  
  for(j=0; j<N; j++)  fc->krvals[j]  = malloc(2*sizeof(double));

  fftw_import_wisdom_from_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");
  
  fc->forwardplan  = fftw_plan_dft_1d(N, fc->input,    fc->cm,  FFTW_FORWARD,  FFTW_MEASURE);
  fc->backwardplan = fftw_plan_dft_1d(N, fc->cmum, fc->output, FFTW_BACKWARD,  FFTW_MEASURE);

  fftw_export_wisdom_to_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");
  
  //um's
  FFTLog_complex z, result;
  
  double L = log(max)-log(min);
  
  fc->kr   = 1.0;
  
  int i;
  
  for(i=0;i<fc->N/2+1;i++){
    z.re   = fc->q;
    z.im   = 2.0*M_PI*(double)i/L;
    
    result = FFTLog_U_mu(mu,z);
    
    //Multiply by (kr)^-2PIim/L
    result.amp *= 1.0;
    result.arg += -2.0*M_PI*(double)i*log(fc->kr)/L;
    
    fc->um[i][0] = result.amp*cos(result.arg);
    fc->um[i][1] = result.amp*sin(result.arg);
  }
  
  //If N even, mutiply by real part only
  if(PARITY(fc->N) == EVEN) fc->um[fc->N/2][1] = 0.0;
  
  return fc;
}


FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z){
  /*Computes 2^z Gamma[(mu + 1 - z)/2]/Gamma[(mu + 1 - z)/2]
              1                2                 3
  */
  
  double amp1, arg1;
  
  gsl_sf_result lnamp2, arg2, lnamp3, arg3;
  
  FFTLog_complex result;
  
  //2^z
  amp1 = exp(z.re*log(2.0));
  
  arg1 = z.im*log(2.0);
  
  //Gamma 1
  FFTLog_complex zplus;
  
  zplus.re = (mu + 1.0 + z.re)/2.0;
  
  zplus.im = z.im/2.0;
  
  gsl_sf_lngamma_complex_e(zplus.re, zplus.im, &lnamp2, &arg2);
  
  //Gamma 2
  FFTLog_complex zminus;
  
  zminus.re = (mu + 1.0 - z.re)/2.0;
  
  zminus.im = - z.im/2.0;
  
  gsl_sf_lngamma_complex_e(zminus.re,zminus.im,&lnamp3,&arg3);

  //Result
  result.amp = amp1*exp(lnamp2.val)*exp(-lnamp3.val);
  
  result.arg = arg1 + arg2.val - arg3.val;
  
  result.re = result.amp*cos(result.arg);
  
  result.im = result.amp*sin(result.arg);
  
  return result;
}


int FFTLog_free(FFTLog_config *fc){
    fftw_destroy_plan(fc->forwardplan);
    fftw_destroy_plan(fc->backwardplan);
  
    fftw_free(fc->pk);
    fftw_free(fc->xi);
   
    fftw_free(fc->input);
    fftw_free(fc->output);
  
    free(fc->krvals);
    
    fftw_free(fc->cm);
    fftw_free(fc->um);  
    fftw_free(fc->cmum);
  
    free(fc);
  
    return 0;
}
