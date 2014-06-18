int CovarianceEigenVecs(){
    // Numerical recipes routines start at 1. 
    order     = hiMultipoleOrder*chiSq_kmaxIndex + 1;
    order2    = order*order;
    
    U         = malloc(order*sizeof(*U));
    eigenVecs = malloc(order*sizeof(*eigenVecs)); 
    sigmaNorm = malloc(order*sizeof(*sigmaNorm));
    
    eigenVals = malloc(order*sizeof(*eigenVals));
    
    for(j=0; j<order; j++){  
        U[j]         = malloc(order*sizeof(**U));
        eigenVecs[j] = malloc(order*sizeof(**eigenVecs));
        sigmaNorm[j] = malloc(order*sizeof(**sigmaNorm));
    }
    
    // Sigma matrix for pre-whitening.
    for(j=0; j<order-1; j++){
        for(k=0; k<order-1; k++)  sigmaNorm[j][k] = 0.0;
        
        sigmaNorm[j][j] = sqrt(Covariance[j][j]);
    }
    
    // Pre-whitening
    for(i=1; i<order; i++){
        for(j=1; j<order; j++){
            U[i][j] = (float) (1./sigmaNorm[i-1][i-1])*Covariance[i-1][j-1]*(1./sigmaNorm[j-1][j-1]); 
        }
    }    

    // printfCovDiag(chiSq_kmaxIndex, hiMultipoleOrder);
    
    fprintfCov(chiSq_kmaxIndex, hiMultipoleOrder);
    
    // printfCov(chiSq_kmaxIndex, hiMultipoleOrder);
    
    // printfSigma(chiSq_kmaxIndex, hiMultipoleOrder);
    
    // printfNormCov(chiSq_kmaxIndex, hiMultipoleOrder);
    
    jacobi(U, order-1, eigenVals, eigenVecs, &nrotations);
    
    printf("\nJacobi rotations: %d", nrotations); 
    
    // eigsrt(eigenVals, eigenVecs, order-1);
    
    for(j=1; j<order; j++){  
        if((lowestKeptEigenvalue      >= eigenVals[j]) && (eigenVals[j] > 0.01)){
            lowestKeptEigenvalue      = eigenVals[j];
        }
    }
    
    printf("\nSmallest retained eigenvalue: %e", lowestKeptEigenvalue);

    printEigens();
    
    // printfE(kBinNumb, hiMultipoleOrder);
    
    // printfET(kBinNumb, hiMultipoleOrder);
    
    // Numerical recipes ordering. 
    xdata    = malloc(order*sizeof(double));
    xtheory  = malloc(order*sizeof(double));
     
    ydata    = malloc(order*sizeof(double));
    ytheory  = malloc(order*sizeof(double));
    
    xdata[0] = 0.0;
    
    /*
    printf("\n\nx and y data: \n");
    for(j=1; j<order; j++)  printf("\n%e \t %e", xdata[j], ydata[j]);
    
    printf("\n\nx and y theory: \n");
    for(j=1; j<order; j++)  printf("\n%e \t %e", xtheory[j], ytheory[j]);
    */
    
    // double chiSq;
     
    // printChiSq();
    
    return 0;
}


int printChiSq(){
    sprintf(filepath, "%s/Data/Covariance/%s_ChiSq.dat", root_dir, surveyType);

    output = fopen(filepath, "w"); 
    
    double chiSq;

    for(k=2; k<order; k++){
        chiSq = 0.0;
    
        for(j=1; j<k; j++)  chiSq += pow(ydata[j] - ytheory[j], 2.)/eigenVals[j]; 
    
        fprintf(output, "%d \t %e \t %e \n", k, eigenVals[k-1], chiSq);
        
        if(eigenVals[k+1] < 0.01)  break;
    }
    
    fclose(output);

    return 0;
}

int printEigens(){
    printf("\n\nEigen-values: \n");
    for(j=1; j<order; j++)  printf("\n%e", eigenVals[j]);
    
    printf("\n\nEigen-vectors: \n\n");
        
    for(k=1; k<order; k++)  printf("%e \t", eigenVecs[k][1]);
    
    printf("\n");
    
    for(k=1; k<order; k++)  printf("%e \t", eigenVecs[k][2]);
    
    printf("\n");
    
    for(k=1; k<order; k++)  printf("%e \t", eigenVecs[k][3]);
    
    return 0;
}


int printfCovDiag(int kBinNumb){
    sprintf(filepath, "%s/Data/Covariance/%s_DiagCov.dat", root_dir, surveyType);

    output = fopen(filepath, "w"); 

    for(k=0; k<chiSq_kmaxIndex; k++){
        // Only half the modes are independent. 
      // fprintf(output, "%e \t %e \t %e \t %d \n", kMultipoles[k], Covariance[k][k], (*pt2Pk)(kMultipoles[k]) + (*pt2shot)(1.), ModeNumber[k]/2);                  
    }

    fclose(output);

    return 0;
}



int fprintfCov(int kBinNumb, int hiMultipoleOrder){
    sprintf(filepath, "%s/Data/Covariance/%s_Covariance.dat", root_dir, surveyType);

    output = fopen(filepath, "w"); 

    for(k=1; k<order; k++){
        for(j=1; j<order; j++){        
            fprintf(output, "%e \t", U[j][k]);                  
        }
    
        fprintf(output, "\n");
    }

    fclose(output);

    return 0;
}


int printfNormCov(int kBinNumb, int hiMultipoleOrder){
    printf("\n\nNormalised covariance matrix.\n");

    for(k=1; k<order; k++){
        for(j=1; j<order; j++){        
            printf("%e \t", U[j][k]);                  
        }
    
        printf("\n");
    }

    return 0;
}


int printfCov(int kBinNumb, int hiMultipoleOrder){
    printf("\n\nCovariance matrix.\n");

    for(k=0; k<order-1; k++){
        for(j=0; j<order-1; j++){        
            printf("%e \t", Covariance[j][k]);                  
        }
    
        printf("\n");
    }

    return 0;
}


int printfSigma(int kBinNumb, int hiMultipoleOrder){
    printf("\n\nSigma matrix.\n");

    for(k=0; k<order-1; k++){
        for(j=0; j<order-1; j++){        
            printf("%e \t", sigmaNorm[j][k]);                  
        }
    
        printf("\n");
    }

    return 0;
}


int printfE(int kBinNumb, int hiMultipoleOrder){
    printf("\n\nE matrix.\n");

    for(k=1; k<order; k++){
        for(j=1; j<order; j++){        
            printf("%e \t", eigenVecs[k][j]);                  
        }
    
        printf("\n");
    }

    return 0;
}


int printfET(int kBinNumb, int hiMultipoleOrder){
    printf("\n\nE^T matrix.\n");

    for(k=1; k<order; k++){
        for(j=1; j<order; j++){        
            printf("%e \t", eigenVecs[j][k]);                  
        }
    
        printf("\n");
    }

    return 0;
}
