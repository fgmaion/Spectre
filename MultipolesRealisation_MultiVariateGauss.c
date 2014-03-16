int multipolesRealisation(){
  // Multivariate Gaussian realisation of Multipoles.
  
  int kBinNumb = 20;
  
  double** CholDecomCov;
  double*  cholVector_y;
  
  double jWaveNumber, mu_j;

  CholDecomCov = malloc((kBinNumb-1)*sizeof(*CholDecomCov));

  for(j=0; j<(kBinNumb-1); j++){
      CholDecomCov[j] = malloc((kBinNumb-1)*sizeof(**CholDecomCov));
  }

  CholeskyDecomposition(Covariance, (kBinNumb-1), CholDecomCov);
  
  mvGauss = malloc((kBinNumb-1)*sizeof(*mvGauss));

  for(j=0; j<(kBinNumb-1); j++){
      mvGauss[j] =  gsl_ran_gaussian(gsl_ran_r, 1.0);
  }  
  
  cholVector_y       = (double *)  malloc((kBinNumb-1)*sizeof(*cholVector_y));

  for(j=0; j<(kBinNumb-1); j++) cholVector_y[j] = 0.0;

  for(j=0; j<(kBinNumb-1); j++){
      for(k=j; k<(kBinNumb-1); k++){
          // printf("\n%e \t %e", Covariance[j][k], CholDecomCov[k][j]);
      
	      cholVector_y[j] += CholDecomCov[j][k]*mvGauss[k];   
      }
  }

  // See Numerical Recipes. 
  sprintf(filepath, "%s/Data/Multipoles/Realisation.dat", root_dir);
  output = fopen(filepath, "w");

  double bestfit;

  for(j=0; j<(kBinNumb-1); j++){
      jWaveNumber = kMultipoles[j];

      mu_j        = (*pt2Pk)(jWaveNumber)*kaiserGauss_Monofactor(jWaveNumber*velDispersion, beta);

      mvGauss[j]  = mu_j + cholVector_y[j];

      // bestfit     = (*pt2Pk)(jWaveNumber)*kaiserGauss_Monofactor(jWaveNumber*1.783333, 0.133333)*(2.5/2.587);

      fprintf(output, "\n%e \t %e \t %e", kMultipoles[j], mvGauss[j], mu_j);

      mvGauss[j] /= 2.587;
  }  

  fclose(output);

  return 0;
}


int CholeskyDecomposition(double** Matrix, int size, double** cholDecomMatrix){
  int  INFO;
  char UPLO;

  double* flat;
  flat = malloc(size*size*sizeof(*flat));
  
  int count = 0;

  for(j=0; j<size; j++){
      for(k=0; k<size; k++){
          flat[count] = Matrix[j][k] ;
          count      += 1;
      }
  } 

  UPLO ='L';

  // Cholesky decomposition. Upper triangular matrix, lower triangular                                                       
  // left as original matrix.                                                                                                

  dpotrf_(&UPLO, &size, flat, &size, &INFO);

  count = 0;

  for(j=0; j<size; j++){
      for(k=0; k<size; k++){
          cholDecomMatrix[j][k] = flat[count];
          count += 1;
      }
  }

  return 0;
}
