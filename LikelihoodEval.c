int minimiseChiSq(){
    minChiSq   = pow(10., 12.);

    ChiSqGrid  = (double ***) malloc(Res*sizeof(*ChiSqGrid));
  
    for(j=0; j<Res; j++)  ChiSqGrid[j] = (double **) malloc(Res*sizeof(**ChiSqGrid));
  
    for(j=0; j<Res; j++){
        for(k=0; k<Res; k++){
            ChiSqGrid[j][k]    = (double *) malloc(Res*sizeof(***ChiSqGrid));
        }
    }

    printf("\nbeta: %.2e, sigma: %.2e, A11 sq.: %.2e", beta, velDispersion, A11Sq);
  
    printf("\n\nBeginning Chi sq. minimisation.");
  
    printf("\n\nPriors:");
    printf("\n%.2f < beta  < %.2f", min_beta, max_beta);  
    printf("\n%.2f < sigma < %.2f", min_velDisperse, max_velDisperse);
    printf("\n%.2f < A11 sq. < %.2f", min_A11Sq, max_A11Sq);
    printf("\n");
    
    printf("\nTrue chi^2: %e", LikelihoodEval());
    
    for(jj=0; jj<Res; jj++){
        beta = min_beta + (max_beta - min_beta)*(jj/dRes);

        for(kk=0; kk<Res; kk++){
            velDispersion = min_velDisperse + (max_velDisperse - min_velDisperse)*(kk/dRes);

            for(ii=0; ii<Res; ii++){
                A11Sq = min_A11Sq + (max_A11Sq - min_A11Sq)*(ii/dRes);
      
                ChiSqGrid[jj][kk][ii] = LikelihoodEval();

                if(ChiSqGrid[jj][kk][ii] < minChiSq){
                    printf("\n%e \t %e \t %e \t %e", beta, velDispersion, A11Sq, ChiSqGrid[jj][kk][ii]);
                
                    minChiSq = ChiSqGrid[jj][kk][ii];
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
    
    int    MonoCount;
    int    QuadCount;
    
    int kBinNumb = 20;
    
    for(MonoCount=1; MonoCount<2; MonoCount++){
        for(QuadCount=1; QuadCount<2; QuadCount++){
            for(j=0; j<(kBinNumb-1); j++){
                for(k=0; k<(kBinNumb-1); k++){
        	        jWaveNumber = kMultipoles[j];
	                kWaveNumber = kMultipoles[k];
    
	                mu_j   = (*pt2Pk)(jWaveNumber)*kaiserGauss_multipole(jWaveNumber*velDispersion, beta, MonoCount);
	                mu_k   = (*pt2Pk)(kWaveNumber)*kaiserGauss_multipole(kWaveNumber*velDispersion, beta, QuadCount);

	                ChiSq += (A11Sq*mvGauss[MonoCount*(kBinNumb-1) + j] - mu_j)*invCov[MonoCount*(kBinNumb-1) + j][QuadCount*(kBinNumb-1) + k]*(A11Sq*mvGauss[QuadCount*(kBinNumb-1) + k] - mu_k);
	            }
	       }
	    }
    }

    // Given covariance, Maximum likelihood solution requires minimisation of chi sq. 

    return ChiSq;
}
