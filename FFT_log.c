#include <assert.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>


int FFTlog_memory(int FFTlogRes, double beta, double velDispersion){
    mono_config          = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);
    quad_config          = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);
     
     W2_config           = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5); 
     
     hex_config          = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 4 + 0.5);

    zero_config          = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);
     
    clipmono_config      = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    clipquad_config      = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    
    convlmonoCorr        = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    convlquadCorr        = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    convlhexCorr         = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 4 + 0.5);  
    
    lnnorm_mono_config   = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    lnnorm_quad_config   = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    
    
    FFTLog_setInput(W2_config,        beta, velDispersion);
    
    
    FFTLog_setInput(mono_config,      beta, velDispersion);
    
    FFTLog_setInput(quad_config,      beta, velDispersion);
    
    FFTLog_setInput( hex_config,      beta, velDispersion);
    
    // FFTLog_setInput(octupole_config, fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(clipmono_config,  beta, velDispersion);
    
    FFTLog_setInput(clipquad_config,  beta, velDispersion);
    
    
    FFTLog_setInput(convlmonoCorr,    beta, velDispersion);
    FFTLog_setInput(convlquadCorr,    beta, velDispersion);
    FFTLog_setInput(convlhexCorr,     beta, velDispersion);
    
    
    FFTLog_setInput(lnnorm_mono_config,    beta, velDispersion);
    FFTLog_setInput(lnnorm_quad_config,    beta, velDispersion);
    
    return 0;
}


int FFTLog_setInput(FFTLog_config *fc, double beta, double velDispersion){  
  double transformOrder;

  double logrmin  = log(fc->min);
  
  double logrmax  = log(fc->max);
  
  double dlogr    = (logrmax - logrmin)/(double)fc->N;
  
  double logrc    = (logrmax + logrmin)/2.0;
  
  double nc       = (double)(fc->N + 1)/2.0 -1;
  
  double logkc    = log(fc->kr)- logrc;
  
  transformOrder  =  fc->mu - 0.5;
  
  // printf("\n\n%d", (int) transformOrder);
  
  // write initial signal
  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]    = exp(logkc + ((double)i-nc)*dlogr);
    
    fc->krvals[i][1]    = exp(logrc + ((double)i-nc)*dlogr);
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.     
    fc->pk[i][0]        = 
                          // spherical_tophat_pk(fc->krvals[i][0]);
                          // haloModel_pk(fc->krvals[i][0], 0.7, transformOrder);
                          (*pt2Pk)(fc->krvals[i][0])*toyRSD_OnePlusOneHalfMuSq(transformOrder);
                          // kaiserLorentz_multipole(fc->krvals[i][0]*velDispersion, beta, (int) transformOrder); 
                          // (*pt2RSD_k)(fc->krvals[i][0]*velDispersion, beta, transformOrder);
    fc->pk[i][1]        = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
    fc->xi[i][0]        = splint_VIPERS_maskMultipoles(fc->krvals[i][1], transformOrder); // (*pt2Xi)(fc->krvals[i][1]);
    fc->xi[i][1]        = 0.0;
  }
  
  return 0;
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
 

int xi_mu(FFTLog_config* fc){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    for(i=0; i<fc->N; i++)  fc->input[i][0] = sqrt(pow(fc->krvals[i][0], 3.)/(8.*pow(pi, 3.)))*fc->pk[i][0];
    for(i=0; i<fc->N; i++)  fc->input[i][1] = 0.0;

    // power law bias may be required. in particular, for Zel'dovich P(k) -> see FFT_log_zeldovich.c
    // for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][0], -fc->q);
    
    FFTLog(fc, fc->forwardplan, fc->backwardplan); 

    // and debias.                                                                                                                        
    // for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][1], -fc->q);

    for(i=0; i<fc->N;   i++) fc->xi[i][0] = pow(-1., transformOrder/2)*fc->output[i][0]*pow(fc->krvals[i][1], -1.5);
    
    return 0;
}  


int pk_mu(FFTLog_config* fc){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    for(i=0; i<fc->N; i++)  fc->input[i][0] = pow(2.*pi*fc->krvals[i][1], 3./2.)*fc->xi[i][0];
    for(i=0; i<fc->N; i++)  fc->input[i][1] = 0.0;
    
    // power law bias may be required.                                                                                                                        
    // for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][1], -fc->q);  

    FFTLog(fc, fc->forwardplan, fc->backwardplan); 

    // and debias.                                                                                                                                           
    // for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][0], -fc->q);

    for(i=0; i<fc->N;   i++) fc->pk[i][0] = pow(-1., transformOrder/2)*fc->output[i][0]*pow(fc->krvals[i][0], -1.5);

    return 0;
}


int varCalc(FFTLog_config* fc, double* sigmaSq, double* u0){
    // Assumes only the monopole contributes to the variance. Certainly quad and hex appear not to.  corr fn. seems to have converged/unaliased between [10**-3, 10**-2.]. Could 'trust'
    // anywhere in this interval.  Variance appears to asymptote as r->0. take variance as value of monopole at ~10**-2.
    
    // Strongly dependent on the limits over which the corr fn. is computed. Best not change them.
    
    for(i=0; i<fc->N;   i++){ 
        if((fc->krvals[i][1]) >= pow(10., -2.)){
            *sigmaSq = fc->xi[i][0];
            
            *u0      = appliedClippingThreshold/sqrt(2.**sigmaSq);
            
            printf("\n\nfft log method: variance: %e, u0:  %e, suppression factor: %e", *sigmaSq, *u0, 0.25*pow(1.0 + gsl_sf_erf(*u0), 2.));
            
            break;
        }   
    }

    return 0;
}


int clipmono(FFTLog_config* clip, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, double u0, double var){
    // Clipped monopole calculation to 2nd order.   
    for(i=0; i<clip->N; i++){
      clip->xi[i][0]   = 0.25*pow(1. + gsl_sf_erf(u0), 2.)*mono->xi[i][0];
             
      clip->xi[i][0]  += (C_n(u0, 1)/var)*(pow(8.*mono->xi[i][0] - 4.*quad->xi[i][0] + 3.*hex->xi[i][0], 2.)/64. - pow(quad->xi[i][0], 2.)/20. + mono->xi[i][0]*(quad->xi[i][0] - 3.*hex->xi[i][0]/4.) + 3.*quad->xi[i][0]*hex->xi[i][0]/8. - (17./576.)*hex->xi[i][0]*hex->xi[i][0]);
    }

    return 0;
}


int clipquad(FFTLog_config* clip, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex, double u0, double var){
    // Clipped monopole calculation to 2nd order.   
    for(i=0; i<clip->N; i++){
      clip->xi[i][0]   = 0.25*pow(1. + gsl_sf_erf(u0), 2.)*quad->xi[i][0];
             
      clip->xi[i][0]  += (C_n(u0, 1)/var)*(2./693.)*(693.*mono->xi[i][0]*quad->xi[i][0] + 99.*pow(quad->xi[i][0], 2.) + 198.*quad->xi[i][0]*hex->xi[i][0] + 50.*pow(hex->xi[i][0], 2.));
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


int add_fogRSD(FFTLog_config *multipole, int Order){
    for(j=0; j<multipole->N; j++){
        if((multipole->krvals[j][0] > 0.0001) && (multipole->krvals[j][0] < 2000.)){
            multipole->pk[j][0] *= 1.; //Lorentz_multipole(multipole->krvals[j][0]*velDispersion, 0);
        }
    }

    return 0;
}


int cnvldmonoCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex){
    // Convolved monopole calculation to 2nd order.   
    for(i=0; i<mono->N; i++){
      cnvld->xi[i][0]   =     mono->xi[i][0]*splint_VIPERS_maskMono(mono->krvals[i][1]);     // splint_windfn_rSpaceMono(mono->krvals[i][1]);
      
      cnvld->xi[i][0]  += 0.2*quad->xi[i][0]*splint_VIPERS_maskQuad(mono->krvals[i][1]);     //splint_windfn_rSpaceQuad(mono->krvals[i][1]);    
    
      // cnvld->xi[i][0]  += 0.2*hex->xi[i][0]*splint_VIPERS_maskHex;
    }

    return 0;
}


int cnvldquadCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex){
    // Convolved quadrupole calculation to 2nd order.   
    for(i=0; i<mono->N; i++){      
      cnvld->xi[i][0]   = mono->xi[i][0]*splint_VIPERS_maskQuad(mono->krvals[i][1]); //splint_windfn_rSpaceQuad(mono->krvals[i][1]);
      
      cnvld->xi[i][0]  += quad->xi[i][0]*(splint_VIPERS_maskMono(mono->krvals[i][1])  + (2./7.)*splint_VIPERS_maskQuad(mono->krvals[i][1]) + (2./7.)*splint_VIPERS_maskHex(mono->krvals[i][1]));  
                       // quad->xi[i][0]*(splint_windfn_rSpaceMono(mono->krvals[i][1]) + (2./7.)*splint_windfn_rSpaceQuad(mono->krvals[i][1]) + (2./7.)*splint_windfn_rSpaceHex(mono->krvals[i][1]));    
    }

    return 0;
}


int cnvldpk_norm(FFTLog_config* fc){    
    for(i=0; i<fc->N;   i++){ 
        if((fc->krvals[i][0]) >= 0.25){
            return i;
     
            // probably never evaluated.        
            break;
        }   
    }
}


int convolvedPkCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
        
    xi_mu(mono_config);

    // xi_mu(quad_config);
    
    // Currently evaluated at second order. 
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config);
    
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config);
    
    // Power spectra given input correlation functions. 
    pk_mu(mono_config);
    
    pk_mu(convlmonoCorr);
    
    pk_mu(convlquadCorr);

    return 0;
}


int clippedPkCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);

    xi_mu(mono_config);

    // xi_mu(quad_config);

    varCalc(mono_config, &variance, &u0);
        
    // Currently evaluated at second order. 
    // clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    
    clipmono(clipmono_config, mono_config, zero_config, zero_config, u0, variance);
    
    // clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // pk_mu(mono_config);
    
    pk_mu(clipmono_config);
    
    // pk_mu(clipquad_config);
    
    // pk_mu(clipquad_config);
    
    // pk_mu(lnnorm_config);

    return 0;
}


int lnnormPkCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    xi_mu(mono_config);
    
    xi_mu(quad_config);
        
    lnnorm(lnnorm_mono_config, mono_config);
    
    pk_mu(lnnorm_mono_config);
    
    add_fogRSD(lnnorm_mono_config, 0);
    
    xi_mu(lnnorm_mono_config);
    
    sprintf(filepath, "%s/Data/lnnorm_anisotropic/fftlog_multipoles_xi_fogRSD.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++)  fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][1], lnnorm_mono_config->xi[j][0], lnnorm_quad_config->xi[j][0]);
    
    fclose(output); 

    return 0;
}


int fog_xiCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
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


int hodmodel_xi(){
    FFTlogRes        = 4096;
    
    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    xi_mu(mono_config);

    xi_mu(quad_config);

    
    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_monoxi.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][1] > 1.) && (mono_config->krvals[j][1] < 300.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0], quad_config->xi[j][0]);
        }
    }
    
    fclose(output); 

    return 0; 
}


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


int VIPERS_mask_cnvldpk(){
    int cnvldpknorm;

    FFTlogRes = 4096;
    
    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    // printf("\n\nStarting fft log clock.");
    // clock_t start = clock(), diff;
    // for(ii=0; ii<1000; ii++){
    
    FFTLog_setInput(mono_config, beta, velDispersion);
    FFTLog_setInput(quad_config, beta, velDispersion);
    
    xi_mu(mono_config);
    xi_mu(quad_config);
    
    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->xi[j][0] = mono_config->xi[j][0]; //*splint_VIPERS_maskMono(mono_config->krvals[j][1]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->xi[j][0] = mono_config->xi[j][0]; // *splint_VIPERS_maskQuad(mono_config->krvals[j][1]);
    // for(j=0; j<mono_config->N; j++)   convlhexCorr->xi[j][0]  = mono_config->xi[j][0]*splint_VIPERS_maskHex(mono_config->krvals[j][1]);
    
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config);
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config);
    
    pk_mu(convlmonoCorr);
    pk_mu(convlquadCorr);
    // pk_mu(convlhexCorr);

    double cnvldpk_zero;

    // convolved pk zero point. 
    for(i=0; i<convlmonoCorr->N;   i++){ 
        if((convlmonoCorr->krvals[i][0]) >= pow(10., -2.)){
            cnvldpk_zero = convlmonoCorr->pk[i][0];
            
            printf("\n\nConvolved P(k) zero point: %e", cnvldpk_zero);
            
            break;
        }   
    }

    // for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = convlmonoCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    // for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = convlquadCorr->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    for(j=0; j<mono_config->N; j++)   convlmonoCorr->pk[j][0] = mono_config->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceMono(convlmonoCorr->krvals[j][0]);
    for(j=0; j<mono_config->N; j++)   convlquadCorr->pk[j][0] = quad_config->pk[j][0] - cnvldpk_zero*splint_VIPERS_kSpaceQuad(convlmonoCorr->krvals[j][0]);

    /*
    sprintf(filepath, "%s/Data/VIPERS_window2/rand_VIPERS_W1_xi_500_mask_0.7_0.8_gridded_hex_normed_jointmultipoles_hihiRes.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][1] > 0.01) && (mono_config->krvals[j][1] < 400.)){
            fprintf(output, "%e \t %e \t %e \t %e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0],  splint_VIPERS_maskMono(mono_config->krvals[j][1]), splint_VIPERS_maskQuad(mono_config->krvals[j][1]), splint_VIPERS_maskHex(mono_config->krvals[j][1]));
        }
    }

    fclose(output);
    */
    // }
    
    // diff = clock() - start;
    // int msec = diff*1000/CLOCKS_PER_SEC;
    // printf("\n\nFFT log approach to convolution, for 1000 models: %d seconds %d milliseconds", msec/1000, msec%1000);
    
    
    sprintf(filepath, "%s/Data/VIPERS_window2/VIPERS_window_500_cnvldpk_hex_aniso_intcor_truth_intcor.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 100.)){
            fprintf(output, "%e \t %e \t %e \t %e \n", mono_config->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0], (*pt2Pk)(mono_config->krvals[j][0]));
        }
    }
    
    fclose(output); 
    
    return 0; 
}


int VIPERS_mask_intCnsrt(){
    //  Hankel transform pair counts of the window. 

    int cnvldpknorm;

    FFTlogRes = 4096;
    
    // Currently beta is set by hand in the input P(k) function. 
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    pk_mu(mono_config);
    
    pk_mu(quad_config);
        
        
    sprintf(filepath, "%s/Data/VIPERS_window2/VIPERS_window_kMultipoles_fftLog.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 100.)){
            fprintf(output, "%e \t %e \t %e \n", mono_config->krvals[j][0], mono_config->pk[j][0]/(1.823239*pow(10., 6.)), quad_config->pk[j][0]/(1.823239*pow(10., 6.)));
        }
    }
    
    fclose(output);
    
    return 0; 
}


/*----------------------------------------------------------------*
 *FFTLog                                                          *
 *----------------------------------------------------------------*/
 
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
  
  for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->output[i][0],fc->output[fc->N-i-1][0]);
  
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
  
  for(j=0; j<N; j++){  
    fc->krvals[j]  = malloc(2*sizeof(double));
  }
  
  fc->forwardplan  = fftw_plan_dft_1d(N, fc->input,    fc->cm,  FFTW_FORWARD,  FFTW_ESTIMATE);
  fc->backwardplan = fftw_plan_dft_1d(N, fc->cmum, fc->output, FFTW_BACKWARD,  FFTW_ESTIMATE);
  
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
