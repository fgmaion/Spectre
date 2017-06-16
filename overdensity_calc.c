int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel           = (int) trunc((xCoor[rowNumber] - AxisLimsArray[0][2])/xCellSize);    
    ylabel           = (int) trunc((yCoor[rowNumber] - AxisLimsArray[0][1])/yCellSize);
    zlabel           = (int) trunc((zCoor[rowNumber] - AxisLimsArray[0][0])/zCellSize);
    
    boxlabel         = (int)                         xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){    
  walltime("Wall time at start of overdensity calc.");

  for(j=0; j<n0*n1*n2; j++)  overdensity[j] = 0.0;
    
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){
        // ngp_assign(xCoor[j], yCoor[j], zCoor[j], (fkp_galweight[j]/sqrt(alpha))*clip_galweight[j]/sampling[j]);
        cic_assign(xCoor[j], yCoor[j], zCoor[j], (fkp_galweight[j]/sqrt(alpha))*clip_galweight[j]/sampling[j]);  // cic assign, weighted by sampling, fkp wghts and clipping.
      }
    }

    for(j=0; j<rand_number; j++){
      // ngp_assign(rand_x[j], rand_y[j], rand_z[j], -sqrt(alpha)*rand_weight[j]);
      cic_assign(rand_x[j], rand_y[j], rand_z[j], -sqrt(alpha)*rand_weight[j]);  // assumes all randoms up to rand_number are accepted.    
    }
      
    walltime("Wall time after overdensity calc.");
    
    return 0;
}


static int rand_int(int n){
  int limit = RAND_MAX - RAND_MAX % n;
  int rnd;

  do{
    rnd = rand();
  }

  while (rnd >= limit);

  return rnd % n;
}


void shuffle(double *array, int n){
  int i, j;

  double tmp;

  for(i = n - 1; i > 0; i--){
    j = rand_int(i + 1);

    tmp = array[j];

    array[j] = array[i];

    array[i] = tmp;
  }
}
