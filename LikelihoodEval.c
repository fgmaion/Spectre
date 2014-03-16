int minimiseChiSq(){
   Res =   200;
  dRes = 200.0;

  minChiSq = pow(10., 6.);
 
  ChiSqGrid  = (double ***) malloc(Res*sizeof(*ChiSqGrid));
  
  for(j=0; j<Res; j++)  ChiSqGrid[j] = (double **) malloc(Res*sizeof(**ChiSqGrid));
  
  for(j=0; j<Res; j++){
      for(k=0; k<Res; k++){
          ChiSqGrid[j][k]    = (double *) malloc(Res*sizeof(***ChiSqGrid));
      }
  }

  printf("\nbeta: %.2e, sigma: %.2e", beta, velDispersion);
  
  printf("\n\nBeginning Chi sq. minimisation.");
  
  for(jj=0; jj<Res; jj++){
      beta = 0.25 + 0.2*(jj/dRes);

      for(kk=0; kk<Res; kk++){
          velDispersion = 1.35 + 0.5*(kk/dRes);

          for(ii=0; ii<Res; ii++){
              A11Sq = 2.2 + 0.8*(ii/dRes);
      
              ChiSqGrid[jj][kk][ii] = LikelihoodEval();

              if(ChiSqGrid[jj][kk][ii] < minChiSq){
                  minChiSq = ChiSqGrid[jj][kk][ii];
                  
                  printf("\n%e \t %e \t %e \t %e", beta, velDispersion, A11Sq, ChiSqGrid[jj][kk][ii]);
              }
	    }
     }
  }

  printf("\n\nMinChiSq: %e", minChiSq);
  
  return 0;
}


double LikelihoodEval(){
    double ChiSq = 0.0;
    
    double mu_j;
    double mu_k;
    
    double jWaveNumber;
    double kWaveNumber;
    
    for(j=0; j<(kBinNumb-1); j++){
        for(k=0; k<(kBinNumb-1); k++){

	        jWaveNumber = kMultipoles[j];
	        kWaveNumber = kMultipoles[k];

	        mu_j = (*pt2Pk)(jWaveNumber)*kaiserGauss_Monofactor(jWaveNumber*velDispersion, beta);
	        mu_k = (*pt2Pk)(kWaveNumber)*kaiserGauss_Monofactor(kWaveNumber*velDispersion, beta);

	        ChiSq += (A11Sq*mvGauss[j] - mu_j)*invCov[j][k]*(A11Sq*mvGauss[k] - mu_k);
	    }
    }

    // Given covariance, Maximum likelihood solution requires minimisation of chi sq. 

    return ChiSq;
}
