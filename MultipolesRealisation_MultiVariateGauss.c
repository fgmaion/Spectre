int multipolesRealisation(int kBinNumb){
  // Multivariate Gaussian realisation of Multipoles.
  
  double** CholDecomCov;
  double*  cholVector_y;
  
  double   mu_j;
  
  double   Mono_j;
  double   Quad_j;

  mvGauss             = malloc(2*(kBinNumb-1)*sizeof(*mvGauss));

  cholVector_y        = malloc(2*(kBinNumb-1)*sizeof(*cholVector_y));

  CholDecomCov        = malloc(2*(kBinNumb-1)*sizeof(*CholDecomCov));

  for(j=0; j<2*(kBinNumb-1); j++){
      CholDecomCov[j] = malloc(2*(kBinNumb-1)*sizeof(**CholDecomCov));
      
      cholVector_y[j] = 0.0;
      
      // Gaussian distributed variable, of sigma = 1.
      mvGauss[j]      =  rand_gaussian(gsl_ran_r, 1.);  
  }

  CholeskyDecomposition(Covariance, 2*(kBinNumb-1), CholDecomCov);

  // See Numerical Recipes. 
  for(j=0; j<2*(kBinNumb-1); j++){
      for(k=j; k<2*(kBinNumb-1); k++){
	      cholVector_y[j] += mvGauss[k]*CholDecomCov[j][k];   
      }
  }

  sprintf(filepath, "%s/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_fullCube_kbin_0.010_000.dat", root_dir, k);

  inputfile = fopen(filepath, "r");

  for(i=0; i<kBinNumb-1; i++){
	fscanf(inputfile, "%*lf \t %lf \t %lf \t %*d \t %*lf \t %*lf \n", &mvGauss[i], &mvGauss[(kBinNumb-1) + i]);
  }

  fclose(inputfile);
  
  
  sprintf(filepath, "%s/Data/Multipoles/zCube_clipThreshold_1.0e+03_subVols_Realisation_beta_%.2f_sigma_%.2f.dat", root_dir, beta, velDispersion);

  output = fopen(filepath, "w");
  
  for(j=0; j<(kBinNumb-1); j++){
      Mono_j                     = (*pt2Pk)(kMultipoles[j])*kaiserLorentz_Monofactor(kMultipoles[j]*velDispersion, beta);
      
      Quad_j                     = (*pt2Pk)(kMultipoles[j])*kaiserLorentz_Quadfactor(kMultipoles[j]*velDispersion, beta);

      // mvGauss[j]                 = Mono_j + cholVector_y[j];
      
      // mvGauss[(kBinNumb-1) + j]  = Quad_j + cholVector_y[(kBinNumb-1) + j];
      
      mvGauss[j]                /= A11Sq;
      
      mvGauss[(kBinNumb-1) + j] /= A11Sq;

      fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", kMultipoles[j], mvGauss[j], mvGauss[(kBinNumb-1) + j], Mono_j, Quad_j, Mono_j/A11Sq, Quad_j/A11Sq);
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
