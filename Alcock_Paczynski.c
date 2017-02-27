int ln_PR(double ln_k){
  // Given ln_k, return ln(P_R(k)), where P_R is the real space power spectrum. 
  double k;

  k = exp(ln_k);

  // Removed dependence on pow(bsigma8,  2.), for MCMC the scaling of dP2_dlnk with bsigma8 should be accounted for. 
  return log( (*pt2Pk)(k)*pow(app_sigma8, -2.));
}


// Numerical evaluation of d ln(P_R) / d ln k with gsl.
int prep_dlnPR_dlnk(){
  gsl_function F;

  // Sanity check, d ln(P_R) / d ln k = 0.96 for large scale power law from inflation.
  double result, abserr, k, ln_k;

  spline_lnk      = malloc(2000*sizeof(*spline_lnk));
  
  dlnPR_dlnk      = malloc(2000*sizeof(*dlnPR_dlnk));
  dlnPR_dlnk_2D   = malloc(2000*sizeof(*dlnPR_dlnk_2D));

  F.function      = &ln_PR;
  F.params        =      0;

  for(j=0; j<2000; j++){
    k             = pow(10., -3. + 4.*j/2000.);

    ln_k          = log(k);

    spline_lnk[j] = ln_k;

    gsl_deriv_central(&F, ln_k, 1e-8, &result, &abserr);
  
    dlnPR_dlnk[j] = result; 

    // printf ("%.10lf \t %.10lf +/- %.10lf\n", k, dlnPR_dlnk[j], abserr);
  }

  spline(spline_lnk, dlnPR_dlnk, 2000, 1.0e31, 1.0e31, dlnPR_dlnk_2D);

  return 0;
}


double splint_dlnPR_dlnk(double k){
  double Interim, lnk;
  
  lnk = log(k);

  splint(spline_lnk, dlnPR_dlnk, dlnPR_dlnk_2D, 2000, lnk, &Interim);

  // stripped of bsigma8 dependence. which is killed by the dlnk derivative anyway. 
  return Interim;
}


// assumes Kaiser-Lorentz model. 
double dP2_dlnk(double k, double beta, double sigma){
  double result;

  double ks, ks2, ks4, PR, P0, P2;
  double M0, M2, M4, M6;

  // kappa in Cole et al. notation.                                                                          
  ks     = k*sigma;
  ks2    =   ks*ks;
  ks4    = ks2*ks2;

  PR     = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0     = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2     = PR*kaiserLorentz_multipole(ks, beta, 2);

  // M0     = muOrderZero(ks);
  // M2     =  muOrderTwo(ks);
  // M4     = muOrderFour(ks);
  // M6     =  muOrderSix(ks);

  if(ks<0.03){
    M0 = 1.0 - 0.166667*ks2 + 0.05*ks4 - 0.0178571*ks4*ks2;
    M2 = 1./3. - 0.1*ks2 + 0.0357143*ks4 - 0.0138889*ks4*ks2;
    M4 = 0.2 - 0.0714286*ks2 + 0.0277778*ks4 - 0.0113636*ks4*ks2;
    M6 = 0.142857 - 0.0555556*ks2 + 0.0227273*ks4 - 0.00961538*ks4*ks2;
  }

  else{
    M0     = muOrderZero(ks);
    M2     = 2.*(1. - M0)/ks2;
    M4     = 2.*(-2. + ks2/3. + 2*M0)/ks4;  
    M6     = 2.*(4. - 2.*ks2/3. + ks4/5. - 4.*M0)/(ks4*ks2);
  }

  result = P2*splint_dlnPR_dlnk(k) + (5./2.)*PR*(  M0 - 3.*(3.- 2.*beta)*M2 - 5.*beta*(6.-beta)*M4 -21.*beta*beta*M6 + 2.*pow(1.+ beta, 2.)/(1.+ks*ks/2.)  );

  // n=4 powerlaw 
  // printf("%.4le \t %.4le \n", result, 4.*P2);

  return result;
}


// assumes Kaiser-Lorentz model.                                                                                                                                                                                 
double dP4_dlnk(double k, double beta, double sigma){
  double result;

  double ks, ks2, ks4, PR, P0, P2, P4;
  double M0, M2, M4, M6, M8;

  // kappa in Cole et al. notation.                                                                                                                                                                               
  ks     = k*sigma;
  ks2    =   ks*ks;
  ks4    = ks2*ks2;

  PR     = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0     = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2     = PR*kaiserLorentz_multipole(ks, beta, 2);
  P4     = PR*kaiserLorentz_multipole(ks, beta, 4);

  // M0     =  muOrderZero(ks);
  // M2     =   muOrderTwo(ks);
  // M4     =  muOrderFour(ks);
  // M6     =   muOrderSix(ks);
  // M8     = muOrderEight(ks);

  if(ks<0.03){
    M0 = 1.0 - 0.166667*ks2 + 0.05*ks4 - 0.0178571*ks4*ks2;
    M2 = 1./3. - 0.1*ks2 + 0.0357143*ks4 - 0.0138889*ks4*ks2;
    M4 = 0.2 - 0.0714286*ks2 + 0.0277778*ks4 - 0.0113636*ks4*ks2;
    M6 = 0.142857 - 0.0555556*ks2 + 0.0227273*ks4 - 0.00961538*ks4*ks2;
  }

  else{
    M0     = muOrderZero(ks);
    M2     = 2.*(1. - M0)/ks2;
    M4     = 2.*(-2. + ks2/3. + 2*M0)/ks4;
    M6     = 2.*(4. - 2.*ks2/3. + ks4/5. - 4.*M0)/(ks4*ks2);
  }

  if(ks>0.18){
    M8     = -16./(ks4*ks4) + (8./3.)/(ks4*ks2) - 0.8/ks4 + (2./7.)/ks2 + 16.*M0/(ks4*ks4);
  }

  else{
    M8     = 1./9. - 0.0454545*ks2 + 0.0192308*ks4 - 0.00833333*ks4*ks2;
  }

  result = P4*splint_dlnPR_dlnk(k) + (9./8.)*PR*( -3.*M0 - 18.*(beta - 5.)*M2 - 5.*(3*beta*beta -60.*beta + 35.)*M4 - 70.*beta*(7.-3.*beta)*M6 - 315.*beta*beta*M8 + 8.*pow(1.+ beta, 2.)/(1.+ks*ks/2.)    );

  // n=4 powerlaw                                                                                                                                                                                                 
  // printf("%.4le \t %.4le \t %.4le \n", result, 4.*P4, M8);                                                                                                                                                              
  return result;
}

// assumes Kaiser-Lorentz model.                                                                                                                                                                      
double dP6_dlnk(double k, double beta, double sigma){
  double result, L;

  double ks, ks2, ks4, PR, P0, P2, P4,  P6;
  double M0, M2, M4, M6, M8, M10;
  double A0, A2, A4, A6, A8, A10;

  // kappa in Cole et al. notation.                                                                                                                                                                               
  ks     = k*sigma;
  ks2    =   ks*ks;
  ks4    = ks2*ks2;

  PR     = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  L      = 1./(1. + 0.5*ks*ks);

  A0     = -5.;
  A2     =  5.*(21. - 2.*beta);
  A4     = -5.*beta*beta + 210.*beta -315.;
  A6     = 231. -630.*beta + 105.*beta*beta;
  A8     = beta*(462. - 315.*beta); 
  A10    = 231.*beta*beta;

  P0     = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2     = PR*kaiserLorentz_multipole(ks, beta, 2);
  P4     = PR*kaiserLorentz_multipole(ks, beta, 4);
  P6     = PR*kaiserLorentz_multipole(ks, beta, 6);

  // M0     =  muOrderZero(ks);
  // M2     =   muOrderTwo(ks);
  // M4     =  muOrderFour(ks);
  // M6     =   muOrderSix(ks);
  // M8     = muOrderEight(ks);
  // M10    =   muOrderTen(ks);

  if(ks<0.03){
    M0 = 1.0 - 0.166667*ks2 + 0.05*ks4 - 0.0178571*ks4*ks2;
    M2 = 1./3. - 0.1*ks2 + 0.0357143*ks4 - 0.0138889*ks4*ks2;
    M4 = 0.2 - 0.0714286*ks2 + 0.0277778*ks4 - 0.0113636*ks4*ks2;
    M6 = 0.142857 - 0.0555556*ks2 + 0.0227273*ks4 - 0.00961538*ks4*ks2;
  }

  else{
    M0     = muOrderZero(ks);
    M2     = 2.*(1. - M0)/ks2;
    M4     = 2.*(-2. + ks2/3. + 2*M0)/ks4;
    M6     = 2.*(4. - 2.*ks2/3. + ks4/5. - 4.*M0)/(ks4*ks2);
  }

  if(ks>0.18){
    M8     = -16./(ks4*ks4) + (8./3.)/(ks4*ks2) - 0.8/ks4 + (2./7.)/ks2 + 16.*M0/(ks4*ks4);
  }

  else{
    M8     = 1./9. - 0.0454545*ks2 + 0.0192308*ks4 - 0.00833333*ks4*ks2;
  }

  if(ks>0.14){
    M10    = (2./315.)*(35./ks2 -90./ks4 +252./pow(ks,6.) -840./pow(ks,8.) +5040./pow(ks,10.) -5040.*M0/pow(ks,10.));
  }

  else{
    M10    = 1./11. - ks2/26. + ks4/60. - ks2*ks4/136. + ks4*ks4/304.;
  }

  result = P6*splint_dlnPR_dlnk(k) + (13./16.)*PR*(-A0*M0 -3.*M2*A2 - 5.*M4*A4 -7.*M6*A6 -9.*M8*A8 -11.*M10*A10 + 16.*pow(1.+beta, 2.)*L);
  
  // n=4 powerlaw                                                                                                                                                                                     
  // printf("%.4le \t %.4le \n", result, 4.*P6); 

  return result;
}


double dPn_dlnk(double k, double beta, double sigma, double PR, double P0, double P2, double P4, double P6, double* dP0_dlnk, double* dP2_dlnk, double* dP4_dlnk, double* dP6_dlnk){
  double result, L, NL, dlnPR_dlnk;

  double ks, ks2, ks4, beta2;
  double M0, M2, M4, M6, M8, M10;
  double A0, A2, A4, A6, A8, A10;

  // kappa in Cole et al. notation.                                                                                                                                                                                                          
  ks     = k*sigma;
  ks2    =   ks*ks;
  ks4    = ks2*ks2;

  dlnPR_dlnk = splint_dlnPR_dlnk(k);

  // PR      = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  L      =  1./(1. + 0.5*ks*ks);
  NL     = pow(1. + beta, 2.)*L;

  A0     = -5.;
  A2     =  5.*(21. - 2.*beta);
  A4     = -5.*beta*beta + 210.*beta -315.;
  A6     = 231. -630.*beta + 105.*beta*beta;
  A8     = beta*(462. - 315.*beta);
  A10    = 231.*beta*beta;

  // P0     = PR*kaiserLorentz_multipole(ks, beta, 0);
  // P2     = PR*kaiserLorentz_multipole(ks, beta, 2);
  // P4     = PR*kaiserLorentz_multipole(ks, beta, 4);
  // P6     = PR*kaiserLorentz_multipole(ks, beta, 6);

  // M0     =  muOrderZero(ks);                                                                                                                                                                                                             
  // M2     =   muOrderTwo(ks);                                                                                                                                                                                                             
  // M4     =  muOrderFour(ks);                                                                                                                                                                                                               // M6     =   muOrderSix(ks);                                                                                                                                                                                                             
  // M8     = muOrderEight(ks);                                                                                                                                                                                                             
  // M10    =   muOrderTen(ks);                                                                                                                                                                                                              
  if(ks<0.03){
    M0 = 1.0 - 0.166667*ks2 + 0.05*ks4 - 0.0178571*ks4*ks2;
    M2 = 1./3. - 0.1*ks2 + 0.0357143*ks4 - 0.0138889*ks4*ks2;
    M4 = 0.2 - 0.0714286*ks2 + 0.0277778*ks4 - 0.0113636*ks4*ks2;
    M6 = 0.142857 - 0.0555556*ks2 + 0.0227273*ks4 - 0.00961538*ks4*ks2;
  }

  else{
    M0     = muOrderZero(ks);
    M2     = 2.*(1. - M0)/ks2;
    M4     = 2.*(-2. + ks2/3. + 2*M0)/ks4;
    M6     = 2.*(4. - 2.*ks2/3. + ks4/5. - 4.*M0)/(ks4*ks2);
  }

  if(ks>0.18){
    M8     = -16./(ks4*ks4) + (8./3.)/(ks4*ks2) - 0.8/ks4 + (2./7.)/ks2 + 16.*M0/(ks4*ks4);
  }

  else{
    M8     = 1./9. - 0.0454545*ks2 + 0.0192308*ks4 - 0.00833333*ks4*ks2;
  }

  if(ks>0.14){
    M10    = (2./315.)*(35./ks2 -90./ks4 +252./(ks4*ks2) -840./(ks4*ks4) +5040./pow(ks,10.) -5040.*M0/pow(ks,10.));
  }

  else{
    M10    = 1./11. - ks2/26. + ks4/60. - ks2*ks4/136. + ks4*ks4/304.;
  }

  *dP0_dlnk = P0*dlnPR_dlnk + PR*(NL - M0 - 6.*beta*M2 - 5.*beta*beta*M4);
  *dP2_dlnk = P2*dlnPR_dlnk + (5./2.)*PR*(  M0 - 3.*(3.- 2.*beta)*M2 - 5.*beta*(6.-beta)*M4 -21.*beta*beta*M6 + 2.*NL  ); 
  *dP4_dlnk = P4*dlnPR_dlnk + (9./8.)*PR*( -3.*M0 - 18.*(beta - 5.)*M2 - 5.*(3*beta*beta -60.*beta + 35.)*M4 - 70.*beta*(7.-3.*beta)*M6 - 315.*beta*beta*M8 + 8.*NL);
  *dP6_dlnk = P6*dlnPR_dlnk + (13./16.)*PR*(-A0*M0 -3.*M2*A2 - 5.*M4*A4 -7.*M6*A6 -9.*M8*A8 -11.*M10*A10 + 16.*NL);

  return result;
}


double dP0_dlnk(double k, double beta, double sigma){
  double result;

  double ks, ks2, ks4, PR, P0, P2;
  double M0, M2, M4, M6;

  // kappa in Cole et al. notation.                                                                          
  ks     = k*sigma;
  ks2    =   ks*ks;
  ks4    = ks2*ks2;


  PR     = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0     = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2     = PR*kaiserLorentz_multipole(ks, beta, 2);

  // M0     = muOrderZero(ks);
  // M2     =  muOrderTwo(ks);
  // M4     = muOrderFour(ks);

  if(ks<0.03){
    M0 = 1.0 - 0.166667*ks2 + 0.05*ks4 - 0.0178571*ks4*ks2;
    M2 = 1./3. - 0.1*ks2 + 0.0357143*ks4 - 0.0138889*ks4*ks2;
    M4 = 0.2 - 0.0714286*ks2 + 0.0277778*ks4 - 0.0113636*ks4*ks2;
  }

  else{
    M0     = muOrderZero(ks);
    M2     = 2.*(1. - M0)/ks2;
    M4     = 2.*(-2. + ks2/3. + 2*M0)/ks4;
  }

  result = P0*splint_dlnPR_dlnk(k) + PR*(pow(1. + beta, 2.)/(1.+ ks*ks/2.) - M0 - 6.*beta*M2 - 5.*beta*beta*M4); 

  // n=4 powerlaw
  // printf("%.4le \t %.4le \n", result, 4.*P0);

  return result;
}

// FFTlog implementation for AP corrected multipoles, see mike.pdf
int FFTLog_updateInput_monoAP(FFTLog_config *fc, double beta, double velDispersion, double alpha_pad, double epsilon_pad){    
  double  kprime;        // k = k'/alpha_pad.
  
  for(i=0; i<fc->N; i++){
    kprime              = fc->krvals[i][0];

    // Now P'(k')                                                                                                                                                                                                
    fc->pk[i][0]        = AP_P0(kprime, beta, velDispersion, alpha_pad, epsilon_pad);
  }
  
  return 0;
}

// FFTlog implementation for AP corrected multipoles, see mike.pdf
int FFTLog_updateInput_quadAP(FFTLog_config *fc, double beta, double velDispersion, double alpha_pad, double epsilon_pad){    
  double  kprime;
  
  for(i=0; i<fc->N; i++){
    // By matching fc->krvals[i][0] to xdata[0] with kvals_matchup(), fc->krvals[i][0] represent k' values.
    kprime              = fc->krvals[i][0];
    
    // Now P'(k')
    fc->pk[i][0]        = AP_P2(kprime, beta, velDispersion, alpha_pad, epsilon_pad);
  }
  
  return 0;
}


// FFTlog implementation for AP corrected multipoles, see mike.pdf                                                                                                                                                
int FFTLog_updateInput_hexAP(FFTLog_config *fc, double beta, double velDispersion, double alpha_pad, double epsilon_pad){
  double  kprime;

  for(i=0; i<fc->N; i++){
    // By matching fc->krvals[i][0] to xdata[0] with kvals_matchup(), fc->krvals[i][0] represent k' values.                                                                                                      
    kprime              = fc->krvals[i][0];

    // Now P'(k')                                                                                                                                                                                                
    fc->pk[i][0]        = AP_P4(kprime, beta, velDispersion, alpha_pad, epsilon_pad);
  }

  return 0;
}


double AP_PN(double kprime, double beta, double sigma, double local_alpha, double epsilon, double* AP_P0, double* AP_P2, double* AP_P4){
  double k, ks, Jacobian, PR, P0, P2, P4, P6, result;

  double dP0_dlnk, dP2_dlnk, dP4_dlnk, dP6_dlnk;

  k        = kprime/local_alpha;
  ks       = k*sigma;

  Jacobian = pow(local_alpha, 3.);

  PR       = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0       = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2       = PR*kaiserLorentz_multipole(ks, beta, 2);
  P4       = PR*kaiserLorentz_multipole(ks, beta, 4);
  P6       = PR*kaiserLorentz_multipole(ks, beta, 6);

  dPn_dlnk(k, beta, sigma, PR, P0, P2, P4, P6, &dP0_dlnk, &dP2_dlnk, &dP4_dlnk, &dP6_dlnk);

  *AP_P0   = (P0 - (2./5.)*epsilon*dP2_dlnk - (6./5.)*epsilon*P2)/Jacobian; 
  *AP_P2   = ((1.- 6.*epsilon/7.)*P2 - 2.*epsilon*dP0_dlnk - (4./7.)*epsilon*dP2_dlnk - (20./7.)*epsilon*P4 - 4.*epsilon*dP4_dlnk/7.)/Jacobian;
  *AP_P4   = (P4 - 3.*epsilon*( (-24./35.)*P2 + (20./77.)*P4 + (210./143.)*P6) - 2.*epsilon*( (18./35.)*dP2_dlnk + (20./77.)*dP4_dlnk + (45./143.)*dP6_dlnk))/Jacobian;

  return 0;
}

int FFTLog_updateInput_multipolesAP(FFTLog_config *mono, FFTLog_config *quad, FFTLog_config *hex, double beta, double velDispersion, double alpha_pad, double epsilon_pad){
  double  kprime;

  double ap_P0, ap_P2, ap_P4;

  for(i=0; i<mono->N; i++){
    // By matching fc->krvals[i][0] to xdata[0] with kvals_matchup(), fc->krvals[i][0] represent k' values.                                                                                                                               
    kprime              = mono->krvals[i][0];

    AP_PN(kprime, beta, velDispersion, alpha_pad, epsilon_pad, &ap_P0, &ap_P2, &ap_P4);

    // Now P'(k')                                                                                                                                                                                                                           
    mono->pk[i][0]      = ap_P0;
    quad->pk[i][0]      = ap_P2;
     hex->pk[i][0]      = ap_P4;
  }

  return 0;
}


double AP_P0(double kprime, double beta, double sigma, double local_alpha, double epsilon){
  // Alcock-Paczynski correction to the monopole, epsilon is the "warping" parameter.
  // See Padmanabhan and White. 

  double k, Jacobian, PR, P0, P2, P4, result;

  k        = kprime/local_alpha;

  Jacobian = pow(local_alpha, 3.);

  PR       = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0       = PR*kaiserLorentz_multipole(k*sigma, beta, 0);
  P2       = PR*kaiserLorentz_multipole(k*sigma, beta, 2);
  P4       = PR*kaiserLorentz_multipole(k*sigma, beta, 4);

  // + O(epsilon^2)
  result   = P0 - (2./5.)*epsilon*dP2_dlnk(k, beta, sigma) - (6./5.)*epsilon*P2; 

  // n=4 powerlaw
  // printf("\n%e \t %e", result, P0 - 14.*epsilon*P2/5.);

  return result/Jacobian;
}


double AP_P2(double kprime, double beta, double sigma, double local_alpha, double epsilon){
  double k, Jacobian, PR, P0, P2, P4, result;
  
  k         = kprime/local_alpha;

  Jacobian  = pow(local_alpha, 3.);

  PR        = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0        = PR*kaiserLorentz_multipole(k*sigma, beta, 0); // sigma, beta are unprimed. 
  P2        = PR*kaiserLorentz_multipole(k*sigma, beta, 2);
  P4        = PR*kaiserLorentz_multipole(k*sigma, beta, 4);

  // + O(epsilon^2)
  result    = (1.- 6.*epsilon/7.)*P2 - 2.*epsilon*dP0_dlnk(k, beta, sigma) - (4./7.)*epsilon*dP2_dlnk(k, beta, sigma) - (20./7.)*epsilon*P4 - 4.*epsilon*dP4_dlnk(k, beta, sigma)/7.; 

  // n=4 powerlaw
  // printf("\n%e \t %e \t %e \t %e \t %e", P0, P2, P4, result, -8.*epsilon*P0 + (1.-22.*epsilon/7.)*P2 - 20.*epsilon*P4/7. - (4.*epsilon/7.)*4.*P4);

  return result/Jacobian;
}


double AP_P4(double kprime, double beta, double sigma, double local_alpha, double epsilon){
  double k, Jacobian, PR, P0, P2, P4, P6, result;

  k         = kprime/local_alpha;

  Jacobian  = pow(local_alpha, 3.);

  PR        = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

  P0        = PR*kaiserLorentz_multipole(k*sigma, beta, 0); // sigma, beta are unprimed.                                                                                                               
  P2        = PR*kaiserLorentz_multipole(k*sigma, beta, 2);
  P4        = PR*kaiserLorentz_multipole(k*sigma, beta, 4);
  P6        = PR*kaiserLorentz_multipole(k*sigma, beta, 6);

  // + O(epsilon^2)                                                                                                                                                                                   
  result    = P4 - 3.*epsilon*( (-24./35.)*P2 + (20./77.)*P4 + (210./143.)*P6) - 2.*epsilon*( (18./35.)*dP2_dlnk(k, beta, sigma) + (20./77.)*dP4_dlnk(k, beta, sigma) + (45./143.)*dP6_dlnk(k, beta, sigma)  );
  
  // n=4 powerlaw, sigma=0                                                                                                                                                                                      
  // printf("\n%e \t %e \t %e \t %e \t %e", P0, P2, P4, result, P4 - epsilon*(216.*P2/35. + 295.*P4/77.  + 90.*P6/13.));

  return result/Jacobian;
}


double n4_beff_multipoles(double beta, double F, double f_perp, double kprime, double sigma, int MonoQuad){
  // n=4 limit of eqn. (A9) of Ballinger, reduces to an effective Kaiser factor -> can use Kaiser Lorentz                   
  // multipoles.                                                                                                             
  double beta_eff;

  beta_eff = (beta + 1.)/(F*F) - 1.;

  // printf("\n%e \t %e \t %e", beta, F, f_perp);

  // printf("\n%d \t %e \t %e", MonoQuad, beta_eff, kaiser_multipole(kprime, beta_eff, MonoQuad));

  return kaiser_multipole(kprime, beta_eff, MonoQuad)*(*pt2Pk)(kprime)*pow(app_sigma8, -2.)*pow(bsigma8,  2.)/(F*pow(f_perp, 3. + 4.));
}


int pw_ap_corrected_multipoles(){
  double f_perp, f_para, beta, k, F;
   
  velDispersion  = 6.00;  // both unprimed 
  beta           = 0.50;

  f_para         = 1.00;  
  f_perp         = 1.00;
  
  // constraint for alpha = 1.
  // f_perp = pow(f_para, -0.5);

  // Flattening factor.
  F             = f_para/f_perp;

  printf("\n\nFlattening factor: %.4lf", F);

  alpha_pad     = pow(pow(f_perp, 2.)*f_para, 1./3.);
  epsilon_pad   = pow(F, 1./3.) - 1.;
  
  alpha_pad     = 1.03;
  epsilon_pad   = 0.03;

  printf("\n\nPadmanabhan & White parameters, alpha: %.4lf \t epsilon: %.4lf\n\n", alpha_pad, epsilon_pad);

  // Simple n=4 power law, corresponds to only an effecive Kaiser factor. See (A9) of Ballinger.  Otherwise, feeds 
  // off usual inputLinearPk() etc.                                                                                              
  // pt2Pk   = &powerlaw;

  prep_dlnPR_dlnk();

  
  int sampling = 200;

  double PR, P0, P2, P4, P6, ap_P0, ap_P2, ap_P4;
  
  sprintf(filepath, "%s/W1_Spectro_V7_2/pw_ap_corrected_multipoles_hex_%.2lf_%.2lf_%.2lf_%.2lf_%.2lf_repeatTest2.dat", root_dir, alpha_pad, epsilon_pad, beta, bsigma8, velDispersion);

  output = fopen(filepath, "w");

  for(jj=0; jj<sampling; jj++){
    // kprime 
    k = pow(10., -2. + 3.*jj/sampling);    

    if(k>0.8)  break;

    AP_PN(k, beta, velDispersion, alpha_pad, epsilon_pad, &ap_P0, &ap_P2, &ap_P4);

    fprintf(output, "%.4le \t %.4le \t %.4le \t %.4le \n", k, ap_P0, ap_P2, ap_P4);
  
    // printf("%.4le \t %.4le \n", k, KaiserLorentz_apMultipoles(beta,F,f_perp,k,sigma,2), AP_P2(k,beta,sigma,epsilon,local_alpha));
  }

  /*
  sprintf(filepath, "%s/W1_Spectro_V7_2/pw_no_ap_multipoles_hex_%.2lf_%.2lf_%.2lf_%.2lf_%.2lf.dat", root_dir, 1.00, 0.00, beta, bsigma8, velDispersion);

  output = fopen(filepath, "w");

  for(jj=0; jj<sampling; jj++){
    // kprime                                                                                                                                                                                                    
    k = pow(10., -2. + 3.*jj/sampling);

    if(k>0.8)  break;

    PR       = (*pt2Pk)(k)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

    P0       = PR*kaiserLorentz_multipole(k*velDispersion, beta, 0);
    P2       = PR*kaiserLorentz_multipole(k*velDispersion, beta, 2);
    P4       = PR*kaiserLorentz_multipole(k*velDispersion, beta, 4);
    P6       = PR*kaiserLorentz_multipole(k*velDispersion, beta, 6);

    fprintf(output, "%.4le \t %.4le \t %.4le \t %.4le \t %.4le \n", k, P0, P2, P4, P6);
  }
  
  fclose(output);
  
  test_mulitpole_calc(f_perp, f_para, beta, velDispersion);
  */
  return 0;
}


double KaiserLorentz_AP(double f_perp, double f_para, double kmod, double mu, double beta, double sigma){
  // eqn. (A8) of Ballinger. Returns P'(k') with appropriate mapping from (kmod, mu) - mode in assumed cosmology,
  // to that in the true cosmology.                             

  double F, J, arg, ap_aniso, eff_kaiser, damping;

  // Flattening factor.                                                                                          
  F = f_para/f_perp;

  // Jacobian = alpha (Xu et al., 2013)                                                                                       
  J = f_perp*f_perp*f_para;

  // Power-spectrum eval. kmod & mu are primed quantities.                                                                                         
  arg = (kmod/f_perp)*sqrt(1. + mu*mu*(pow(F, -2.) -1.));

  // AP-anisotropy factor. No RSD dependence.                                                                    
  ap_aniso = pow(1. + mu*mu*(pow(F, -2.) -1.), -2.);

  // Effective Kaiser factor, RSD dependence.                                                                    
  eff_kaiser = pow(1. + mu*mu*((1.+beta)*pow(F, -2.) - 1.), 2.);

  // Effective damping.                                                                                          
  damping = 1./(1. + 0.5*pow(kmod*mu*sigma/f_para, 2.));

  return (*pt2Pk)(arg)*pow(app_sigma8, -2.)*pow(bsigma8,  2.)*ap_aniso*eff_kaiser*damping/J;
}


int test_mulitpole_calc(double f_perp, double f_para, double beta, double sigma){
  // Alcock-Paczynski effect yields a new P'(') which differs in amplitude, by (1/f_perp^2)*(1/f_para) factor, 
  // eqn (A5) of Ballinger and otherwise is the old power spectrum evaluated at a different mode, given by factors f_perp 
  // and f_para.

  int m0, m1, m2;
  
  double pk;

  // Ballinger amplitude term for P(k).
  double Jacobian;

  Jacobian = f_perp*f_perp*f_para; 

  // Unprimed variables - corresponding to those values in the true cosmology.  
  double  noprime_kx, noprime_ky, noprime_kz, noprime_kSq, noprime_kmod, noprime_mu;

  polar_pkcount = 0;

  // prep_pkRegression(-2., log10(modkMax), kbin_no);

  // assign2DPkMemory();

  for(k=0; k<n0; k++){    
    for(j=0; j<n1; j++){
      for(i=0; i<n2; i++){
	    m0 = k;
	    m1 = j;
	    m2 = i;

       	    if(m2>n2/2)  m2                   -= n2;
	    if(m1>n1/2)  m1                   -= n1;
	    if(m0>n0/2)  m0                   -= n0;

       	    // Assumed geometry, correspond to primed quantities in Ballinger.
	    k_x                                = kIntervalx*m2;
	    k_y                                = kIntervaly*m1;
	    k_z                                = kIntervalz*m0;

    	    Index                              = k*n1*n2 + j*n2 + i;

    	    // Primed, 
	    kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

    	    kmodulus                           = pow(kSq, 0.5);
    
    	    mu                                 = k_z/kmodulus;
	    if(kmodulus < 0.000001)       mu   = 0.0;

    	    // For mode assumed cosmology k_x, calculate corresponding mode in the true cosmology. 
	    // true                            // fiducial
	    noprime_kx                         = k_x/f_perp;
	    noprime_ky                         = k_y/f_perp;
            noprime_kz                         = k_z/f_para;

    	    noprime_kSq                        = pow(noprime_kx, 2.) + pow(noprime_ky, 2.) + pow(noprime_kz, 2.);

      	    noprime_kmod                       = pow(noprime_kSq, 0.5);

	    noprime_mu                         = noprime_kz/noprime_kmod;

	    // Observed pk in the presence of the AP effect. 
	    pk                                 = (*pt2Pk)(noprime_kmod)*pow(1. + beta*noprime_mu*noprime_mu, 2.)*pow(app_sigma8, -2.)*pow(bsigma8,  2.);

	    pk                                /= 1.+ 0.5*pow(noprime_kmod*noprime_mu*sigma, 2.);

	    // eqn. (A8) of Ballinger. Analytic solution for Kaiser-Lorentz model p(k) in presence of AP effect.
	    // pk                                 = KaiserLorentz_AP(f_perp, f_para, kmodulus, mu, beta, sigma);

	    // Negative quadrupole test.
	    // pk                              = (2. - 3.*0.5*(3*mu*mu - 1.));
  
	    // Multipoles calculated for primed (observed) variables
            polar_pk[polar_pkcount][0]         =           kmodulus;
	    polar_pk[polar_pkcount][1]         =           fabs(mu);
	    polar_pk[polar_pkcount][2]         =        pk/Jacobian;

	    // twodim_pk[polar_pkcount][0]        = fabs(k_z);
	    // twodim_pk[polar_pkcount][1]        = pow(k_y*k_y + k_x*k_x, 0.5);
	    // twodim_pk[polar_pkcount][2]        = pk;

	    polar_pkcount                     += 1;
      }
    }
  }

  sprintf(filepath, "%s/W1_Spectro_V7_2/AP_Multipoles_Ballinger_pk_prediction_hex_%.2lf_%.2lf_%.2lf_%.2lf_%.2lf.dat", root_dir, alpha_pad, epsilon_pad, beta, bsigma8, velDispersion);

  // MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 1);

  HexadecapoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, kHexadecapole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 1);

  // Cartesian2Dpk(polar_pkcount);
  
  return 0;
}

