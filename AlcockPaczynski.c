double ln_PR(double ln_k, void* p){
  double k;
  
  k = exp(ln_k);

  return log((*pt2Pk)(k));  // Given ln_k, return ln(P_R(k)), where P_R is the real space power spectrum.
}

int prep_dlnPR_dlnk(){
  gsl_function F;
  double       result, abserr, k, ln_k;  // Sanity check, d ln(P_R) / d ln k = 0.96 for large scale power law from inflation.
  
  dlnPR_dlnk      = malloc(FFTlogRes*sizeof(*dlnPR_dlnk));

  F.function      = &ln_PR;
  F.params        =      0;

  for(j=0; j<FFTlogRes; j++){
    gsl_deriv_central(&F, log(mono_config->krvals[j][0]), 1e-8, &result, &abserr);
  
    dlnPR_dlnk[j] = result; 

    // printf ("%.10lf \t %.10lf +/- %.10lf\n", ln_k, dlnPR_dlnk[j], abserr);
  }

  return 0;
}

int apmultipoles(FFTLog_config *mono, FFTLog_config *quad, FFTLog_config *hex, double beta, double velDispersion, double alpha, double epsilon){
  if(epsilon == 0.){
    // no ap geometric distortion. 
    FFTlog_updatepk(mono, quad, hex, beta, velDispersion);
  }

  else{
    int i;

    double kprime;

    double A0, A2, A4, A6, A8, A10;
    double M0, M2, M4, M6, M8, M10;
    double k, ks, ks2, ks4, beta2, L, NL, Jacobian, PR, P0, P2, P4, P6, dP0_dlnk, dP2_dlnk, dP4_dlnk, dP6_dlnk, result;
  
    beta2    = beta*beta;
  
    A0       = -5.;
    A2       =  5.*(21. - 2.*beta);
    A4       = -5.*beta2 + 210.*beta -315.;
    A6       = 231. -630.*beta + 105.*beta2;
    A8       = beta*(462. - 315.*beta);
    A10      = 231.*beta2;

    Jacobian = pow(alpha, 3.);
  
    for(i=0; i<mono->N; i++){
      kprime          = mono->krvals[i][0];  // By matching fc->krvals[i][0] to xdata[0] with kvals_matchup(), fc->krvals[i][0] represent k' values. 
    
      k               = kprime/alpha;

      ks              = k*velDispersion;

      ks2             = ks*ks;
      ks4             = ks2*ks2;
    
      PR              = FFTlog_Pk[i]*pow(bsigma8,  2.);

      L               = 1./(1. + ks2/2.);
      NL              = pow(1. + beta, 2.)*L;
    
      M0              = muOrderZero(ks);
      M2              = muOrderTwo(ks);
      M4              = muOrderFour(ks);
      M6              = muOrderSix(ks);
      M8              = muOrderEight(ks);
      M10             = muOrderTen(ks);
    
      P0              = PR*(M0 + 2.*beta*M2 + beta2*M4);
      P2              = PR*((5./2.)*(-M0 + M2*(3. - 2.*beta) + M4*(-beta2 + 6.*beta) + 3.*beta2*M6));
      P4              = PR*((9./8.)*(35.*beta2*M8  + 10.*beta*(7. -3.*beta)*M6 + (35. - 60.*beta + 3.*beta2)*M4 + 6.*(beta - 5.)*M2 + 3.*M0));
      P6              = 0.0;
    
      dP0_dlnk        = P0*dlnPR_dlnk[i] + PR*(NL - M0 - 6.*beta*M2 - 5.*beta2*M4);
      dP2_dlnk        = P2*dlnPR_dlnk[i] + (5./2.)*PR*(  M0 - 3.*(3.- 2.*beta)*M2 - 5.*beta*(6.-beta)*M4 -21.*beta2*M6 + 2.*NL);
      dP4_dlnk        = P4*dlnPR_dlnk[i] + (9./8.)*PR*( -3.*M0 - 18.*(beta - 5.)*M2 - 5.*(3*beta2 -60.*beta + 35.)*M4 - 70.*beta*(7.-3.*beta)*M6 - 315.*beta2*M8 + 8.*NL);
      dP6_dlnk        = P6*dlnPR_dlnk[i] + (13./16.)*PR*(-A0*M0 -3.*M2*A2 - 5.*M4*A4 -7.*M6*A6 -9.*M8*A8 -11.*M10*A10 + 16.*NL);
      
      // Now P'(k')
      mono->pk[i][0]  = (P0 - (2./5.)*epsilon*dP2_dlnk - (6./5.)*epsilon*P2)/Jacobian;;
      quad->pk[i][0]  = ((1.- 6.*epsilon/7.)*P2 - 2.*epsilon*dP0_dlnk - (4./7.)*epsilon*dP2_dlnk - (20./7.)*epsilon*P4 - 4.*epsilon*dP4_dlnk/7.)/Jacobian;;
       hex->pk[i][0]  = (P4 - 3.*epsilon*( (-24./35.)*P2 + (20./77.)*P4 + (210./143.)*P6) - 2.*epsilon*( (18./35.)*dP2_dlnk + (20./77.)*dP4_dlnk + (45./143.)*dP6_dlnk))/Jacobian;
    }
  }
  
  return 0;
}
