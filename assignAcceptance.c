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

    // Cannot apply acceptance based on ra & dec cuts currently, due to the application of Acceptance during
    // Comoving number density calc., which includes both W1 and W4.

    for(j=0; j<Vipers_Num; j++){
        if((redshiftLowLimit<zUtilized[j]) && (zUtilized[j]<redshiftHiLimit) && (M_B[j]<absMagCut)){
            Acceptanceflag[j]  = true;
        }
    }

    /*
    printf("\n\nAccepted limits (in the VIPERS basis..)");
    printf("\nx max:  %lf \t x min:  %lf", AcceptedMax(xCoor, Acceptanceflag, Vipers_Num), AcceptedMin(xCoor, Acceptanceflag, Vipers_Num));
    printf("\ny max:  %lf \t y min:  %lf", AcceptedMax(yCoor, Acceptanceflag, Vipers_Num), AcceptedMin(yCoor, Acceptanceflag, Vipers_Num));
    printf("\nz max:  %lf \t z min:  %lf", AcceptedMax(zCoor, Acceptanceflag, Vipers_Num), AcceptedMin(zCoor, Acceptanceflag, Vipers_Num));
    */
    
    printf("\n\nVIPERS gals. meeting redshift acceptance.");
    printf("\nx max:  %e \t x min:  %e",     AcceptedMax(xCoor, Acceptanceflag, Vipers_Num), AcceptedMin(xCoor, Acceptanceflag, Vipers_Num));
    printf("\ny max:  %e \t y min:  %e",     AcceptedMax(yCoor, Acceptanceflag, Vipers_Num), AcceptedMin(yCoor, Acceptanceflag, Vipers_Num));
    printf("\nz max:  %e \t z min:  %e",     AcceptedMax(zCoor, Acceptanceflag, Vipers_Num), AcceptedMin(zCoor, Acceptanceflag, Vipers_Num));
    printf("\nchi max:  %e \t chi min:  %e", AcceptedMax(rDist, Acceptanceflag, Vipers_Num), AcceptedMin(rDist, Acceptanceflag, Vipers_Num));
        
    int AcceptedNumber = 0;
    
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
            AcceptedNumber += 1;
        }
    }

    printf("\nNumber of accepted VIPERS gals.: %d", AcceptedNumber);

    return 0;
}


int assignAcceptanceCube(){
    int subVolAccepted = 0;

    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j] = false;
    
    for(j=0; j<Vipers_Num; j++){
        if((AxisLimsArray[0][0]<xCoor[j]) && (xCoor[j]<AxisLimsArray[1][0])){
            if((AxisLimsArray[0][1]<yCoor[j]) && (yCoor[j]<AxisLimsArray[1][1])){
                if((AxisLimsArray[0][2]<zCoor[j]) && (zCoor[j]<AxisLimsArray[1][2])){
                    Acceptanceflag[j]  = true;
                    
                    subVolAccepted    +=    1;
                }
            }
        }
    }

    printf("\nNumber of galaxies in sub volume: %d", subVolAccepted);

    return 0;
}
