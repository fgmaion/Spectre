#include <assert.h>
#include <string.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#define CC 299792.458

// Pk model, power law with high k truncation.
double Pk_powerlaw_truncated(double k, double k0, double kc, int n){
    double dn = (double) n;

    return 2.*pow(pi, 2.)*pow(k0, -3.)*pow(k/k0, dn)*exp(-1.*k/kc);
}
 
double Pk_powerlaw_truncated_xi(double r, double k0, double kc, int n){
    double dn = (double) n;
    double y  = kc*r;
    
    return pow(kc/k0, dn + 3.)*sin((2.+dn)*atan(y))*tgamma(2.+dn)*pow(y*pow(1.+y*y, 1.+dn/2.), -1.);
} 
 
double Pk_powerlaw(double k, double r0, double gamma){
    return 4.*pi*pow(k, -3.)*pow(k*r0, gamma)*tgamma(2.-gamma)*sin((pi/2.)*(2.-gamma)); 
} 

double Pk_powerlaw_xi(double r, double r0, double gamma){
    return pow(r0/r, gamma);
}

double kaiser_multipole_xifactors(double gamma, int monoQuad){
    switch(monoQuad){
        case 0:
            return 1.0;
        case 2:
            return -gamma*pow(3.-gamma, -1.);
        case 4:
            return gamma*(2.+gamma)*pow(-5.+gamma, -1.)*pow(-3.+gamma, -1.);
    }
}

double pk_Setup(double k, double transformOrder){
    return sqrt(pow(k, 3.)/(8.*pow(pi, 3.)))*(*pt2Pk)(k)*kaiserLorentz_multipole(k*velDispersion, beta, (int) transformOrder);
}

double xi_Setup(double r, double transformOrder){
    return pow(2.*pi*r, 3./2.)*Pk_powerlaw_xi(r, 5., 1.8)*kaiser_multipole(r, beta, (int) transformOrder)*kaiser_multipole_xifactors(1.8, transformOrder);
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
  
  printf("\n\nTransform order: %.2e\n", transformOrder);
  
  // write signal
  for(i=0; i<fc->N; i++){
    k[i] = exp(logkc + ((double)i-nc)*dlogr);
    
    r[i] = exp(logrc + ((double)i-nc)*dlogr);
    
    // purely real.  Set P(k) to obtain xi(r) by inverse Hankel transform. 
    //fc->pk[i][0]  = sqrt(pow(k[i], 3.)/(8.*pow(pi, 3.)))*Pk_powerlaw(k[i], 5., 1.8)*kaiser_multipole(k[i], beta, (int) transformOrder);
    
    fc->pk[i][0]    = sqrt(pow(k[i], 3.)/(8.*pow(pi, 3.)))*Pk_powerlaw(k[i], 5., 1.8)*kaiser_multipole(k[i], beta, (int) transformOrder);
    fc->pk[i][1]    = 0.0;
    
    // purely real.  Set xi(r) to obtain P(k) by Hankel transform.  
    fc->xi[i][0]    = pow(2.*pi*r[i], 3./2.)*Pk_powerlaw_xi(r[i], 5., 1.8);
    fc->xi[i][1]    = 0.0;
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


int print_fcprop(FFTLog_config *fc){
    printf("\n\nP(k) object properties:");
    printf("\nminimum: %e", fc->min);
    printf("\nmaximum: %e", fc->max);
    printf("\nbias   : %d", fc->q);
    printf("\norder:   %d",   fc->mu);
    printf("\nlength: %d",  fc->N);
    
    return 0;
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
