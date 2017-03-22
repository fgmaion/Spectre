void FFTLog(FFTLog_config *fc, fftw_plan p_forward, fftw_plan p_backward){
  /* Hamilton 2000. http://casa.colorado.edu/~ajsh/FFTLog/                                                                                                                                  
     The FFTLog algorithm for taking the discrete Hankel transform, equation (22),                                                                                                          
     of a sequence an of N logarithmically spaced points is:                                                                                                                                
     -- FFT an to obtain the Fourier coefficients cm, equation (13);                                                                                                                         
     -- multiply by um given by equations (18) and (19) to obtain cm um;                                                                                                                    
     -- FFT cm um back to obtain the discrete Hankel transform ãn, equation (21).                                                                                                           
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
  for(i=1; i<fc->N/2+1; i++){                                                                                                                                                                    fc->cmum[i][0] = (fc->cm[i][0]*fc->um[i][0] + fc->cm[i][1]*fc->um[i][1])*pow(fc->um[i][0]*fc->um[i][0] + fc->um[i][1]*fc->um[i][1], -1.);                                                
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

