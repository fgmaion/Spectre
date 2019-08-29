// created 10/02/2017
int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][2])/xCellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/yCellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][0])/zCellSize);
    
    boxlabel                = (int)                         xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){    
    /*
    double  eff_weights_norm = 0.0;
    double* eff_weights;

    eff_weights = malloc(Vipers_Num*sizeof(double));

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){
        eff_weights[j]      = fkp_galweight[j]*clip_galweight[j]/sampling[j];

        eff_weights_norm   += eff_weights[j]; 
      }
    }
    */
    /*
    double* randomisedweights;

    randomisedweights = malloc(accepted*sizeof(double));

    int cc = 0; 

    // double blah[5] = {1, 2, 3, 4, 5};

    // shuffle(blah, 5);

    // for(j=0; j<5; j++)  printf("shuffled array: %.1lf", blah[j]);

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){
	randomisedweights[cc] = clip_galweight[j];

	cc += 1;  
      }
    }
    
    shuffle(randomisedweights, accepted);
    */
    // Clean overdensity before galaxy assignment. 
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;
    
    // cc = 0;

    // cic assign the galaxies to overdensity, weighted by sampling, fkp weights and clipping weight.
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){  
	    cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], (1./sampling[j])*fkp_galweight[j]*clip_galweight[j]);

	    // cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], (1./sampling[j])*fkp_galweight[j]*randomisedweights[cc]);

	    // cc += 1;

	    // cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], eff_weights[j]/eff_weights_norm);
      
	    // printf("\n%.2lf \t %.2lf \t %.2lf", 1./sampling[j], fkp_galweight[j], clip_galweight[j]);
      }
    }
       
    for(j=0; j<n0*n1*n2; j++)  surveyMask[j] = 0.0;
        
    for(j=0; j<rand_number; j++){
        // assign randoms to surveyMask weighted by fkp.
        if(rand_accept[j] == true)     cic_assign(0, rand_x[j], rand_y[j], rand_z[j], rand_weight[j]);    
    }
    
    // free if not pair counting window. 
    // free(rand_x);
    // free(rand_y);
    // free(rand_z); 
    
    return 0;
}


static int rand_int(int n) {
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do {
    rnd = rand();
  }

  while (rnd >= limit);

  return rnd % n;
}


void shuffle(double *array, int n) {
  int i, j;

  double tmp;

  for(i = n - 1; i > 0; i--){
    j = rand_int(i + 1);

    tmp = array[j];

    array[j] = array[i];

    array[i] = tmp;
  }
}
