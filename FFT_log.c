#include <assert.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define CC 299792.458

int FFTlog_memory(int FFTlogRes, double beta, double velDispersion){
    mono_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);
    quad_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);
     hex_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 4 + 0.5);

    // octupole_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 6 + 0.5);
     
    clipmono_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    clipquad_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    
    convlmonoCorr   = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    convlquadCorr   = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    
    
    FFTLog_setInput(mono_config,      beta, velDispersion);
    
    FFTLog_setInput(quad_config,      beta, velDispersion);
    
    FFTLog_setInput( hex_config,      beta, velDispersion);
    
    // FFTLog_setInput(octupole_config, fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(clipmono_config,  beta, velDispersion);
    
    FFTLog_setInput(clipquad_config,  beta, velDispersion);
    
    
    FFTLog_setInput(convlmonoCorr,    beta, velDispersion);
    FFTLog_setInput(convlquadCorr,    beta, velDispersion);
    
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
  
  // write initial signal
  for(i=0; i<fc->N; i++){
    fc->krvals[i][0]    = exp(logkc + ((double)i-nc)*dlogr);
    
    fc->krvals[i][1]    = exp(logrc + ((double)i-nc)*dlogr);
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.     
    fc->pk[i][0]        = (*pt2Pk)(fc->krvals[i][0])*toyRSD_OnePlusOneHalfMuSq(transformOrder); // (*pt2RSD_k)(fc->krvals[i][0]*velDispersion, beta, transformOrder);
    fc->pk[i][1]        = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
    fc->xi[i][0]        = (*pt2Xi)(fc->krvals[i][1]);
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
    
    FFTLog(fc, fc->forwardplan, fc->backwardplan); 

    for(i=0; i<fc->N;   i++) fc->xi[i][0] = pow(-1., transformOrder/2)*fc->output[i][0]*pow(fc->krvals[i][1], -1.5);
    
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


int pk_mu(FFTLog_config* fc){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    for(i=0; i<fc->N; i++)  fc->input[i][0] = pow(2.*pi*fc->krvals[i][1], 3./2.)*fc->xi[i][0];
    for(i=0; i<fc->N; i++)  fc->input[i][1] = 0.0;
    
    FFTLog(fc, fc->forwardplan, fc->backwardplan); 

    for(i=0; i<fc->N;   i++) fc->pk[i][0] = pow(-1., transformOrder/2)*fc->output[i][0]*pow(fc->krvals[i][0], -1.5);

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
  
  // printf("\n\nAnother round.");
  // for(j=0; j<fc->N/2 + 1; j++)  printf("\n%e \t %e", fc->cm[j][0], fc->cm[j][1]);
  
  fc->cmum[0][0] = fc->cm[0][0]*fc->um[0][0] - fc->cm[0][1]*fc->um[0][1];
  fc->cmum[0][1] = fc->cm[0][0]*fc->um[0][1] + fc->cm[0][1]*fc->um[0][0];

  for(i=1; i<fc->N/2+1; i++){
    fc->cmum[i][0] = fc->cm[i][0]*fc->um[i][0] - fc->cm[i][1]*fc->um[i][1];
    fc->cmum[i][1] = fc->cm[i][0]*fc->um[i][1] + fc->cm[i][1]*fc->um[i][0];
    
    //Hermitian symetry (i.e. to get a real signal after FFT back)
    fc->cmum[fc->N-i][0] =  fc->cmum[i][0];
    fc->cmum[fc->N-i][1] = -fc->cmum[i][1];
  }
  
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

/*
FFTLog_config *limrange_create(FFTLog_config* fc, int korr, double krlo, double krhi){
    // korr = 0 if an FFTlog object corresponding to a mid range in k is wanted, 1 if a mid range in r.   
    // The idea here is to isolate a mid range over which either p(k) or xi(r) is unaffected by aliasing/end effects before an inverse transformation. 

    int firstgood  = fc->N; 
    int  lastgood  =     0;
    int goodpoints =     0;
    
    krlo           =    pow(10., krlo);
    krhi           =    pow(10., krhi);
    
    for(i=0; i<fc->N; i++){
        if((fc->krvals[i][korr] >= krlo)  && (i <  firstgood))                firstgood   = i;
    
        if((fc->krvals[i][korr] >= krlo)  && (fc->krvals[i][korr] <= krhi))   goodpoints += 1; 
        
        if((fc->krvals[i][korr] <= krhi)  && (i >   lastgood))                lastgood    = i;
    } 
    
    if(goodpoints%2 != 0)  goodpoints += 1;
  
    FFTLog_config *midfc = (FFTLog_config*)  malloc(sizeof(FFTLog_config));
  
    //FFTW3 Initialization
    midfc->min          = fc->krvals[firstgood][korr];
    midfc->max          = fc->krvals[lastgood][korr];
  
    midfc->q            = fc->q;
    midfc->mu           = fc->mu;
    midfc->N            = goodpoints; // array indexing 0 to (goodpoints-1).
  
    midfc->pk           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
    midfc->xi           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);

    midfc->cm           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
    midfc->um           = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
    midfc->cmum         = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
 
    midfc->input        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
    midfc->output       = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*goodpoints);
  
    midfc->krvals       = malloc(sizeof(double*)*goodpoints);
  
    for(j=0; j<goodpoints; j++)  midfc->krvals[j]  = malloc(2*sizeof(double));

    midfc->forwardplan  = fftw_plan_dft_1d(goodpoints, midfc->input,    midfc->cm,  FFTW_FORWARD,  FFTW_ESTIMATE);
    midfc->backwardplan = fftw_plan_dft_1d(goodpoints, midfc->cmum, midfc->output, FFTW_BACKWARD,  FFTW_ESTIMATE);
  
    // copy selected signal over required range.
    for(i=firstgood; i<=lastgood; i++){
      midfc->krvals[i - firstgood][0]    = fc->krvals[i][0];
    
      midfc->krvals[i - firstgood][1]    = fc->krvals[i][1];
    
      // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform.     
      midfc->pk[i - firstgood][0]        = fc->pk[i][0]; //*(*pt2RSD_k)(fc->kvals[i]*velDispersion, beta, transformOrder);
      midfc->pk[i - firstgood][1]        = 0.0;
    
      // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
      midfc->xi[i - firstgood][0]        = fc->xi[i][0];
      midfc->xi[i - firstgood][1]        = 0.0;
    }
  
    //um's. Requires arguments int N, double min, double max, double q, double mu. all of which can be provided for the midfc object. 
    FFTLog_complex z, result;
  
    double L = log(midfc->max)-log(midfc->min);
  
    midfc->kr   = 1.0;
  
    int i;
  
    for(i=0;i<midfc->N/2+1;i++){
      z.re   = midfc->q;
      z.im   = 2.0*M_PI*(double)i/L;
    
      result = FFTLog_U_mu(midfc->mu,z);
    
      //Multiply by (kr)^-2PIim/L
      result.amp *= 1.0;                                                      // This seems HIGHLY suspicious? //
      result.arg += -2.0*M_PI*(double)i*log(midfc->kr)/L;
    
      midfc->um[i][0] = result.amp*cos(result.arg);
      midfc->um[i][1] = result.amp*sin(result.arg);
    }
  
    //If N even, mutiply by real part only
    if(PARITY(midfc->N) == EVEN) midfc->um[midfc->N/2][1] = 0.0;

    return midfc;
}
*/

int print_fftlogdetails(FFTLog_config* fc){
    printf("\nmin:     %e", fc->min);
    printf("\nmin:     %e", fc->max);

    printf("\nbias, q: %d", fc->q);
    printf("\nmu:      %e", fc->mu);
  
    printf("\nRes:     %d", fc->N);

    // copy selected signal over required range.
    for(i=0; i< fc->N; i++)  printf("\n%e \t %e \t %e \t %e \t %e \t %e", fc->krvals[i][0], fc->krvals[i][1], fc->pk[i][0], fc->xi[i][0], fc->um[i][0], fc->um[i][1]);

  return 0;
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


int clippedPkCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    // Spectral distortion calculation. 
    // Correlation functions given input P(k).
    xi_mu(mono_config);

    xi_mu(quad_config);
    
    // xi_mu( hex_config);
    
    // limmono_config   = limrange_create(mono_config, 1, -4., 3.);
    
    // print_fftlogdetails(mono_config);
    
    // print_fftlogdetails(limmono_config);
    
    // varCalc(mono_config, &variance, &u0);
        
    // Currently evaluated at second order. 
    // clipmono(clipmono_config, mono_config, quad_config, hex_config, u0, variance);
    
    // clipquad(clipquad_config, mono_config, quad_config, hex_config, u0, variance);
    
    // Power spectra given input correlation functions. 
    // pk_mu(mono_config);
    
    // pk_mu(clipmono_config);
    
    // pk_mu(clipquad_config);
    
    printf("\nClipped mono evaluated.");
    
    print_xi();
    
    // printClippedPk();

    return 0;
}


int cnvldmonoCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex){
    // Clipped monopole calculation to 2nd order.   
    for(i=0; i<mono->N; i++){
      // printf("%e \t %e \n", cnvld->krvals[i][1], mono->krvals[i][1]);
      
      cnvld->xi[i][0]   = mono->xi[i][0]*splint_windfn_rSpaceMono(mono->krvals[i][1]);
      
      cnvld->xi[i][0]  += quad->xi[i][0]*splint_windfn_rSpaceQuad(mono->krvals[i][1]);    
    }

    return 0;
}


int cnvldquadCorr(FFTLog_config* cnvld, FFTLog_config* mono, FFTLog_config* quad, FFTLog_config* hex){
    // Clipped monopole calculation to 2nd order.   
    for(i=0; i<mono->N; i++){
      // printf("%e \t %e \n", cnvld->krvals[i][1], mono->krvals[i][1]);
      
      cnvld->xi[i][0]   = mono->xi[i][0]*splint_windfn_rSpaceQuad(mono->krvals[i][1]);
      
      cnvld->xi[i][0]  += quad->xi[i][0]*(0.2*splint_windfn_rSpaceMono(mono->krvals[i][1]) + (2./7.)*splint_windfn_rSpaceQuad(mono->krvals[i][1]));    
    }

    return 0;
}


int convolvedPkCalc(){
    FFTlogRes        = 4096;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    // Spectral distortion calculation. 
    // Correlation functions given input P(k).
    xi_mu(mono_config);

    xi_mu(quad_config);
    
    xi_mu( hex_config);

    // Currently evaluated at second order. 
    cnvldmonoCorr(convlmonoCorr, mono_config, quad_config, hex_config);
    
    cnvldquadCorr(convlquadCorr, mono_config, quad_config, hex_config);
    
    // Power spectra given input correlation functions. 
    pk_mu(mono_config);
    
    pk_mu(convlmonoCorr);
    
    pk_mu(convlquadCorr);
        
    printf("\nConvolved mono + quad evaluated.");
    
    print_xi();
    
    printCnvldPk();

    return 0;
}



int print_xi(){
    sprintf(filepath, "%s/Data/SpectralDistortion/fftlog_hodpk_NoRSD_%d_%.2e_%.2e_02Sep_xi.dat", root_dir, FFTlogRes, mono_config->min, mono_config->max);

    output = fopen(filepath, "w");

    for(i=0; i<mono_config->N; i++) fprintf(output, "%e \t %e \t %e \t %e \n", mono_config->krvals[i][1], mono_config->xi[i][0], quad_config->xi[i][0], mono_config->xi[i][0]*splint_windfn_rSpaceMono(mono_config->krvals[i][1]));
    
    fclose(output);

    return 0;
}


int printClippedPk(){
    sprintf(filepath, "%s/Data/SpectralDistortion/fftlog_hodpk_NoRSD_%d_%.2e_%.2e_08Aug_clip_pk.dat", root_dir, FFTlogRes, clipmono_config->min, clipmono_config->max);

    output = fopen(filepath, "w");

    for(i=0; i<clipmono_config->N; i++) fprintf(output, "%e \t %e \t %e \t %e \t %e\n", clipmono_config->krvals[i][0], clipmono_config->pk[i][0], clipquad_config->pk[i][0], 0.25*pow(1. + gsl_sf_erf(u0), 2.)*mono_config->pk[i][0], 0.25*pow(1. + gsl_sf_erf(u0), 2.)*quad_config->pk[i][0]);
    
    fclose(output);

    return 0;
}


int printCnvldPk(){
    sprintf(filepath, "%s/Data/SpectralDistortion/fftlog_hodpk_NoRSD_%d_%.2e_%.2e_04Sep_cnvld_firstOrder_pk.dat", root_dir, FFTlogRes, mono_config->min, mono_config->max);

    output = fopen(filepath, "w");

    for(i=0; i<mono_config->N; i++){ 
        fprintf(output, "%e \t %e \t %e \t %e \n", mono_config->krvals[i][0], mono_config->pk[i][0], convlmonoCorr->pk[i][0], convlquadCorr->pk[i][0]);
    }
    
    fclose(output);

    return 0;
}
