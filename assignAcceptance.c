double AcceptedMax(double a[], bool b[], int n){
  double max = -99.;

  for(j=0; j<n; j++){
    if(b[j] == true){
      max = a[j];
      break;
    }
  }

  for(j=0; j<n; j++){
    if((b[j] == true) &&  (a[j] > max)){
      max = a[j];
    }
  }

  return max;
}


double AcceptedMin(double a[], bool b[], int n){
  double min = 99.;

  for(j=0; j<n; j++){
    if(b[j] == true){
      min = a[j];
      break;
    }
  }

  for(j=0; j<n; j++){
    if((b[j]==true) && (a[j]<min)){
      min = a[j];
    }
  }

  return min;
}


int assignAcceptance_rand(){
    for(j=0; j<rand_number; j++)                        rand_accept[j]  = false;
    
    for(j=0; j<rand_number; j++){
        if((loChi<rand_chi[j]) && (rand_chi[j]<hiChi))  rand_accept[j]  = true;
    }
        
    accepted_rand = 0;
    
    for(j=0; j<rand_number; j++){
        if(rand_accept[j] == true)  accepted_rand += 1;
    }

    printf("\n\nTotal number of randoms %d, accepted %d", rand_number, accepted_rand);
    
    return 0;
}


int assignAcceptance_WX_SPECTRO_V7(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = false;
    
    int accepted_count = 0;
    
    for(j=0; j<Vipers_Num; j++){
        if((lo_zlim<=zobs[j]) && (zobs[j]<=hi_zlim)){
            if((2<=zflag[j]) && (zflag[j]<=9)){
                if(photoMask[j] == 1){
                    // Inclusive limits. 
                    Acceptanceflag[j]  = true;
                    
                    accepted_count    +=    1;
                }
            }
        }
    }
        
    printf("\n\nNumber of accepted galaxies: %d", accepted_count);


    double daccepted_gals     = 0.0;

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true)  daccepted_gals   += 1./sampling[j];
    }

    accepted_gals = (int) round(daccepted_gals);

    printf("\n\nTotal number of galaxies %d, accepted (weighted) %.2lf", Vipers_Num, daccepted_gals);

    alpha = 1.*daccepted_gals/accepted_rand;

    printf("\n\nalpha %.4f", alpha);

    printf("\n\naccepted. inverted, rotated & translated");
        
    return 0;
}


int assignAcceptance(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = false;
    
    for(j=0; j<Vipers_Num; j++){
        if((lo_zlim<gal_z[j]) && (gal_z[j]<hi_zlim)){
            if((lo_MBlim<M_B[j]) && (M_B[j]<hi_MBlim)){
                Acceptanceflag[j]  = true;
            }
        }
    }
        
    double daccepted_gals     = 0.0;
    
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true)  daccepted_gals   += 1./sampling[j];
    }

    printf("\n\nTotal number of galaxies %d, accepted (weighted) %.2lf", Vipers_Num, daccepted_gals);
    
    alpha = 1.*daccepted_gals/accepted_rand;
    
    printf("\n\nalpha %.4f", alpha);
                                                                                                                                                                              
    printf("\n\naccepted. inverted, rotated & translated");  
                     
    return 0;
}


int assignAcceptanceCube(){
    accepted_gals = 0;
    
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
            accepted_gals += 1;
        }
    }

    alpha = 1.*accepted_gals/accepted_rand;

    printf("\n\nAccepted gals:%d \t, rands: %d", accepted_gals, accepted_rand);

    printf("\n\nalpha: %.4e", alpha);
    
    return 0;
}
