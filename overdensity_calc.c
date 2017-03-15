int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][2])/xCellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/yCellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][0])/zCellSize);
    
    boxlabel                = (int)                         xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){    
    for(j=0; j<n0*n1*n2; j++){
      // overdensity[j][0] = 0.0;  // Clean before galaxy/random assignment.
      // overdensity[j][1] = 0.0;

      overdensity[j] = 0.0;
    }
    
    // Overkill: #pragma omp parallel for private(j)
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){  
        ngp_assign(xCoor[j], yCoor[j], zCoor[j], fkp_galweight[j]*clip_galweight[j]/sampling[j]);
        // cic_assign(1, xCoor[j], yCoor[j], zCoor[j], fkp_galweight[j]*clip_galweight[j]/sampling[j]);  // cic assign, weighted by sampling, fkp wghts and clipping.
      }
    }
    
    #pragma omp parallel for private(j) if(thread == 1)
    for(j=0; j<rand_number; j++){
      ngp_assign(rand_x[j], rand_y[j], rand_z[j], -alpha*rand_weight[j]);
      // cic_assign(1, rand_x[j], rand_y[j], rand_z[j], -alpha*rand_weight[j]);  // assumes all randoms up to rand_number are accepted.    
    }

    printf("\n\nOverdensity:");
    
    for(j=0; j<n0*n1*n2; j++)  printf("\n%.6lf", overdensity[j]);
    
    walltime("Wall time after cloud-in-cell assignment");
    
    return 0;
}

/*
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
*/
