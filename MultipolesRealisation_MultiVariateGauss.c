int multipolesRealisation(){
  // Multivariate Gaussian realisation of Multipoles.
  
  int kBinNumb = 20;
  
  double** CholDecomCov;
  double*  cholVector_y;
  
  double jWaveNumber, mu_j;

  CholDecomCov = malloc(2*(kBinNumb-1)*sizeof(*CholDecomCov));

  for(j=0; j<2*(kBinNumb-1); j++){
      CholDecomCov[j] = malloc(2*(kBinNumb-1)*sizeof(**CholDecomCov));
  }

  CholeskyDecomposition(Covariance, 2*(kBinNumb-1), CholDecomCov);
  
  /*
  for(j=0; j<2*(kBinNumb-1); j++){
    for(k=0; k<2*(kBinNumb-1); k++){
        
        printf("%d \t %d \t %e \t %e \n", j, k, Covariance[j][k], CholDecomCov[j][k]);
    }  
  }
  */
  
  mvGauss = malloc(2*(kBinNumb-1)*sizeof(*mvGauss));

  for(j=0; j<2*(kBinNumb-1); j++){
      mvGauss[j]     =  gsl_ran_gaussian(gsl_ran_r, 1.0);
  }  
  
  cholVector_y       = (double *)  malloc(2*(kBinNumb-1)*sizeof(*cholVector_y));

  for(j=0; j<2*(kBinNumb-1); j++) cholVector_y[j] = 0.0;

  for(j=0; j<2*(kBinNumb-1); j++){
      for(k=j; k<2*(kBinNumb-1); k++){
          // printf("\n%e \t %e", Covariance[j][k], CholDecomCov[k][j]);
      
	      cholVector_y[j] += CholDecomCov[j][k]*mvGauss[k];   
      }
  }

  // See Numerical Recipes. 

  double bestfit;
  double Mono_j;
  double Quad_j;
  
  sprintf(filepath, "%s/Data/Multipoles/Realisation.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<(kBinNumb-1); j++){
      jWaveNumber = kMultipoles[j];

      Mono_j      = (*pt2Pk)(jWaveNumber)*kaiserGauss_Monofactor(jWaveNumber*velDispersion, beta);

      mvGauss[j]  = Mono_j; // + cholVector_y[j];

      bestfit     = (*pt2Pk)(jWaveNumber)*kaiserGauss_Monofactor(jWaveNumber*velDispersion, beta)/A11Sq;

      mvGauss[j] /= A11Sq;

      fprintf(output, "\n%e \t %e \t %e", jWaveNumber, mvGauss[j], bestfit);
  }  
  
  for(j=0; j<(kBinNumb-1); j++){
      jWaveNumber = kMultipoles[j];

      Quad_j      = (*pt2Pk)(jWaveNumber)*kaiserGauss_Quadfactor(jWaveNumber*velDispersion, beta);

      mvGauss[(kBinNumb-1) + j]  = Quad_j; //  + cholVector_y[j];

      mvGauss[(kBinNumb-1) + j] /= A11Sq;
      
      bestfit     = (*pt2Pk)(jWaveNumber)*kaiserGauss_Quadfactor(jWaveNumber*velDispersion, beta)/A11Sq;
      
      fprintf(output, "\n%e \t %e \t %e", jWaveNumber, mvGauss[(kBinNumb-1) + j], bestfit);
  }  
 
  fclose(output);

  return 0;
}


int CholeskyDecomposition(double** Matrix, int size, double** cholDecomMatrix){
  int  INFO;
  
  // Either 'U' for upper triangular, or 'L' for lower triangular.
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

  // Return a lower triangular matrix. 
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
