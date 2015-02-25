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


int assignAcceptance(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = false;
    
    for(j=0; j<Vipers_Num; j++){
        if((lo_zlim<gal_z[j]) && (gal_z[j]<hi_zlim) && (lo_MBlim<M_B[j]) && (M_B[j]<hi_MBlim)){
            Acceptanceflag[j]  = true;
        }
    }
        
    accepted_gals = 0;
    
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
            accepted_gals += 1;
        }
    }

    printf("\n\nTotal number of galaxies %d, accepted %d", Vipers_Num, accepted_gals);
    
    return 0;
}


int assignAcceptanceCube(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j] = true;

    return 0;
}
