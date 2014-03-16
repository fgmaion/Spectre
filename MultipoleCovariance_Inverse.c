int CovarianceInverse(){
    // Inverse of C_ij using LAPACK.
    
    int mockNumber = 64;
    int kBinNumb   = 20;
    
    // double Covariance[3*(kBinNumb-1)][3*(kBinNumb-1)];

    double* flatCov;
    flatCov     = malloc((kBinNumb-1)*(kBinNumb-1)*sizeof(*flatCov));

    invCov      = (double **) malloc((kBinNumb-1)*sizeof(*invCov));

    for(j=0; j<(kBinNumb-1); j++){
      invCov[j] = malloc((kBinNumb-1)*sizeof(**invCov));
    }

    int count =0;

    for(j=0; j<(kBinNumb-1); j++) flatCov[j] = 0.0;
 
    // Assuming diagonal covariance. 
    for(j=0; j<(kBinNumb-1); j++){
        for(k=0; k<(kBinNumb-1); k++){
            flatCov[count] = Covariance[j][k];
		  
		  count += 1;
	    }

    }
    
    int    order, order2;
    
    order  = (kBinNumb-1);
    order2 = order*order;

    // length of N+1
    double pivotArray[order + 1];
    
    int    errorHandler; 

    double lapackWorkspace[order2];
    
    dgetrf_(&order, &order, flatCov, &order, pivotArray, &errorHandler);
    
    dgetri_(&order, flatCov, &order, pivotArray, lapackWorkspace, &order2, &errorHandler);
    
    count = 0;

    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
            invCov[j][k] = flatCov[count];
            
           	count += 1; 
	        
	        // printf("%e \t", invCov[j][k]);
        }
    
        // printf("\n");
    }
    
    return 0;
}
