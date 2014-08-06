#include <assert.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define CC 299792.458

int FFTlog_memory(int FFTlogRes, double beta, double velDispersion){
    fftlogr         = malloc(FFTlogRes*sizeof(double)); 
    fftlogk         = malloc(FFTlogRes*sizeof(double)); 

    mono_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);
    quad_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);
     hex_config     = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 4 + 0.5);

    // octupole_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 6 + 0.5);
     
    clipmono_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 0 + 0.5);  
    clipquad_config = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), 0.0, 2 + 0.5);  
    
    FFTLog_setInput(hex_config,      fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(mono_config,     fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(quad_config,     fftlogk, fftlogr, beta, velDispersion);
    
    // FFTLog_setInput(octupole_config, fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(clipmono_config, fftlogk, fftlogr, beta, velDispersion);
    
    FFTLog_setInput(clipquad_config, fftlogk, fftlogr, beta, velDispersion);
    
    return 0;
}


int FFTLog_setInput(FFTLog_config *fc, double *k, double *r, double beta, double velDispersion){
  // Set k and r arrays. 
  
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
    k[i] = exp(logkc + ((double)i-nc)*dlogr);
    
    r[i] = exp(logrc + ((double)i-nc)*dlogr);
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform. 
    //fc->pk[i][0]  = sqrt(pow(k[i], 3.)/(8.*pow(pi, 3.)))*Pk_powerlaw(k[i], 5., 1.8)*kaiser_multipole(k[i], beta, (int) transformOrder);
    
    fc->pk[i][0]    = sqrt(pow(k[i], 3.)/(8.*pow(pi, 3.)))*(*pt2Pk)(k[i]); //*(*pt2RSD_k)(k[i]*velDispersion, beta, transformOrder);
    fc->pk[i][1]    = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
    fc->xi[i][0]    = pow(2.*pi*r[i], 3./2.)*(*pt2Xi)(r[i]);
    fc->xi[i][1]    = 0.0;
  }

  return 0;
}
 
 
int xi_mu(FFTLog_config* fc, int mu, double* fftlogr){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    FFTLog(fc, fc->pk_forwardplan, fc->pk_backwardplan); 
    
    // reverse array
    for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->xi[i][0],fc->xi[fc->N-i-1][0]);
  
    for(i=0; i<fc->N;   i++) fc->xi[i][0] = (fc->xi[i][0]/(double)fc->N);

    for(i=0; i<fc->N;   i++) fc->xi[i][0] = pow(-1., mu/2)*fc->xi[i][0]*pow(fftlogr[i], -1.5);

    // FFTLog_free(fc);    
    
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

int pk_mu(FFTLog_config* fc, int mu, double* fftlogk){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    FFTLog(fc, fc->xi_forwardplan, fc->xi_backwardplan); 
    
    // reverse array
    for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->pk[i][0],fc->pk[fc->N-i-1][0]);
  
    for(i=0; i<fc->N;   i++) fc->pk[i][0] = (fc->pk[i][0]/(double)fc->N);

    for(i=0; i<fc->N;   i++) fc->pk[i][0] = pow(-1., mu/2)*fc->pk[i][0]*pow(fftlogk[i], -1.5);

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
  
  //xi's: FFT backward
  fftw_execute(p_backward);
  
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
 
  fc->pk_forwardplan  = fftw_plan_dft_1d(N, fc->pk,   fc->cm, FFTW_FORWARD,  FFTW_ESTIMATE);
  fc->pk_backwardplan = fftw_plan_dft_1d(N, fc->cmum, fc->xi, FFTW_BACKWARD, FFTW_ESTIMATE);
  
  fc->xi_forwardplan  = fftw_plan_dft_1d(N, fc->xi,   fc->cm, FFTW_FORWARD,  FFTW_ESTIMATE);
  fc->xi_backwardplan = fftw_plan_dft_1d(N, fc->cmum, fc->pk, FFTW_BACKWARD, FFTW_ESTIMATE);
  
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


void FFTLog_free(FFTLog_config *fc){
  fftw_destroy_plan(fc->pk_forwardplan);
  fftw_destroy_plan(fc->pk_backwardplan);
  
  fftw_destroy_plan(fc->xi_forwardplan);
  fftw_destroy_plan(fc->xi_backwardplan);
  
  fftw_free(fc->pk);
  fftw_free(fc->xi);
  
  fftw_free(fc->cm);
  fftw_free(fc->um);  
  fftw_free(fc->cmum);
  
  free(fc);
  return;
}
