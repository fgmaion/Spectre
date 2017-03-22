int set_u0(){
  if(data_mock_flag == 0){                                                                                                                                 
    // if(d0      == 1000.)  A11Sq  =      1.00; // Where is the redshift distinction; evolution of bias.                           
    // else if(d0 ==   10.)  A11Sq  =   1./1.30;                                                                                    
    // else if(d0 ==    6.)  A11Sq  =   1./1.85;                                                                                    
    // else if(d0 ==    4.)  A11Sq  =   1./3.00;                                                                                    
                                                                                                                                                          
    if(d0      == 1000.)  A11Sq  =      1.00;                                                                                       
    else if(d0 ==    4.)  A11Sq  =   1./4.50;                                                                                       
  }                                                                                                                                                        
                                                                                                                                                          
  if(data_mock_flag == 1){                                                                                                                                 
    if(fieldFlag == 1){                                                                                                                                    
      if(d0      == 1000.)  A11Sq  =      1.00;  // Correct monopole.
      else if(d0 ==   10.)  A11Sq  =   1./1.30;  // Note: as a difference estimator, covariance does not have to explicity corrected.
      else if(d0 ==    6.)  A11Sq  =   1./1.70;                                                                                      
      else if(d0 ==    4.)  A11Sq  =   1./2.70;                                                                                      
    }                                                                                                                                                      
                                                                                                                                                           
    if(fieldFlag == 4){                                                                                                                                   
      if(d0      == 1000.)  A11Sq  =      1.00;  // Correct monopole.                                                              
      else if(d0 ==   10.)  A11Sq  =   1./1.30;  // Note: as a difference estimator, covariance does not have to explicity corrected.
      else if(d0 ==    6.)  A11Sq  =   1./1.90;
      else if(d0 ==    4.)  A11Sq  =   1./2.70;
    }
  }

  u0 = inverse_erf(2.*sqrt(A11Sq) - 1.);  // Estimate u0 from amplitude suppression.

  clipmono_amp   = 0.25*pow(1. + gsl_sf_erf(u0), 2.); // assumes u0 is fixed.
  clip_distcoeff = C_n(u0, 1);
  
  printf("\n\nu0 set: %.3lf", u0);
  
  return 0;
}


int snipping_amplitudeCorrection(double clipped_pk[], int len){
  // Note I:  As a difference estimator, covariance does not have to explicity corrected for shotnoise.
  // Note II: Assumes same rescaling for each field.

  for(j=0; j<len; j++){
    if(hi_zlim == 1.20){
      if(d0 == 10.){
        clipped_pk[j]              *=   1.30;
        clipped_pk[j + mono_order] *=   1.30;
      }
    }

    else if(d0 ==  6.){
      clipped_pk[j]              *=   1.85;
      clipped_pk[j + mono_order] *=   1.85;
    }

    else if(d0 ==  4.){
      clipped_pk[j]              *=   3.00;
      clipped_pk[j + mono_order] *=   3.00;
    }

    else if(hi_zlim == 1.0){
      if(d0 ==  4.){
        clipped_pk[j]              *=   5.00;
        clipped_pk[j + mono_order] *=   5.00;
      }
    }
  }

  return 0;
}


int snipping_shotnoise(double clipped_pk[], int len){
  for(j=0; j<len; j++){
    if(fieldFlag == 1){
      if(hi_zlim == 1.20){
        if(d0      == 1000.)  clipped_pk[j] -= 277.16;  // Correct monopole.
        else if(d0 ==   10.)  clipped_pk[j] -= 243.94;  // Note: as a difference estimator, covariance does not have to explicity corrected.
        else if(d0 ==    6.)  clipped_pk[j] -= 208.99;
        else if(d0 ==    4.)  clipped_pk[j] -= 167.90;
      }
    }
    
    else if(hi_zlim == 1.00){
      if(d0      == 1000.)    clipped_pk[j] -= 564.00;
      else if(d0 ==    4.)    clipped_pk[j] -= 227.85;
    }
    
    if(fieldFlag == 4){
      if(d0      == 1000.)    clipped_pk[j] -= 289.16;  // Correct monopole.
      else if(d0 ==   10.)    clipped_pk[j] -= 249.26;  // Note: as a difference estimator, covariance does not have to be shot noise corrected.
      else if(d0 ==    6.)    clipped_pk[j] -= 208.00;
      else if(d0 ==    4.)    clipped_pk[j] -= 168.59;
    }
  }
  
  return 0;
}
