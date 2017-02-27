int CovarianceInverse(){
    // Inverse of C_ij using LAPACK.
    
    int mockNumber = 64;
    int kBinNumb   = 20;
    
    // double Covariance[3*(kBinNumb-1)][3*(kBinNumb-1)];

    double* flatCov;
    flatCov     = malloc(2*(kBinNumb-1)*2*(kBinNumb-1)*sizeof(*flatCov));

    invCov      = (double **) malloc(2*(kBinNumb-1)*sizeof(*invCov));

    for(j=0; j<2*(kBinNumb-1); j++){
      invCov[j] = malloc(2*(kBinNumb-1)*sizeof(**invCov));
    }

    int count   = 0;

    for(j=0; j<4*(kBinNumb-1); j++) flatCov[j] = 0.0;
 
    for(j=0; j<2*(kBinNumb-1); j++){
        for(k=0; k<2*(kBinNumb-1); k++){
            flatCov[count] = Covariance[j][k];
		  
   		    count += 1;
	    }
    }

    int    order, order2;
    
    order  = 2*(kBinNumb-1);
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
        }
    }
    
    sprintf(filepath, "%s/Data/Covariance/Clipped_HODCubeColumns_InverseCovariance.dat", root_dir);

    output = fopen(filepath, "w"); 

    for(j=0; j<2*(kBinNumb-1); j++){
        for(k=0; k<2*(kBinNumb-1); k++){
            fprintf(output, "%e \t", fabs(invCov[j][k]));                  
        }
    
        fprintf(output, "\n");
    }

    fclose(output);
    
    return 0;
}
