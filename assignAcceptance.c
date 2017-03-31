int vollim_cutbyMB(double lim){
  for(j=0; j<Vipers_Num; j++){
    if(M_B[j] > lim){
      Acceptanceflag[j]  = false;
    }
  }

  return 0;
}


int assignAcceptance(){
  for(j=0; j<Vipers_Num; j++){
    Acceptanceflag[j]  = false;
    
    if((lo_zlim <= gal_z[j]) && (gal_z[j] <= hi_zlim)){
      if(data_mock_flag == 1){
        if(dec[j] >= -5.97)  Acceptanceflag[j]  = true;  // dec problem in mocks; cut data to have the same boundary; should specify W4 is not cut.  
      }

      else{
        Acceptanceflag[j]  = true;
      }
    }
  }
  
  return 0;
}


int assignAcceptance_rand(){
  //** randoms bool flag is obsolete. 
  //** Assumes random catalogue is sorted by dec, most negative last.  
  if(data_mock_flag == 0){ // dec problem in mocks. 
    for(j=0; j<rand_number; j++){
      // randoms satisfy chi limits by construction.
      if(rand_dec[j] >= -5.97){
	rand_number = j;  // dec problem in mocks.

	break;
	// fkp_accepted_rand += pow(1. + (*pt2nz)(rand_chi[j])*fkpPk, -1.);
      }
    }

    accepted_rand = rand_number; // Superfluous. 
  }

  else{
    // assumes randoms have same boundary as (cut) data. (beyond dec = -5.97)
    accepted_rand = rand_number;
  }

  printf("\n\nTotal number of randoms %d, accepted %d", rand_number, accepted_rand);

  return 0;
}


int accepted_gal(){
  accepted      =   0;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      accepted           += 1;
    }
  }
  
  return 0;
}


int alpha_calc(){
   accepted      =   0;
  daccepted_gals = 0.0;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      accepted           += 1;

      // daccepted_gals     += 1./sampling[j];
      daccepted_gals     += clip_galweight[j]/sampling[j]; //  13/05/16.
    }
  }

  alpha = 1.*daccepted_gals/accepted_rand;
  
  printf("\nInfo: d0=% 5d;  Total number of galaxies on input: %d, accepted: %d, accepted(weighted): % 6.2lf. (1./alpha): %.9lf; FKP norm: %.6lf", d0, Vipers_Num,
                                                                                                                   accepted, daccepted_gals, 1./alpha, sqrt(alpha)*bare_fkp_norm);
  
  return 0;
}


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


int assignAcceptance_true(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = true;

    return 0;
}


double assignAcceptance_parent(){
  int accepted = 0;

  for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = false;

  for(j=0; j<Vipers_Num; j++){
    if((lo_zlim<=gal_z[j]) && (gal_z[j]<=hi_zlim) && (ra[j] > LowerRAlimit) && (ra[j] < UpperRAlimit) && (LowerDecLimit < dec[j]) && (dec[j] < UpperDecLimit)){
      if(data_mock_flag ==0){
	// dec problem in mocks.                                                                                                                                                                                                         
	if(dec[j] >= -5.97)  Acceptanceflag[j]  = true;
      }

      else{
	Acceptanceflag[j]  = true;
      }
    }
  }

  double daccepted_gals     = 0.0;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      accepted             += 1;
    
      daccepted_gals       += 1.;
    }
  }

  printf("\n\nTotal number of galaxies on input: %d, accepted: %.2lf", Vipers_Num, daccepted_gals);

  return daccepted_gals;
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
