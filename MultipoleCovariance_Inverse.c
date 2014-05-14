int CovarianceInverse(int kBinNumb){
    // Inverse of C_ij using LAPACK.
    
    /*
    int order, order2;

    order  = 2;
    order2 = order*order;
    
    double** A;
    double** A_inv;
    double** A_co;
    
    double*  A_flat;
    
    A        = malloc(order*sizeof(*A));
    A_inv    = malloc(order*sizeof(*A_inv));
    A_co     = malloc(order*sizeof(*A_co));
    
    A_flat   = malloc(order2*sizeof(*A_flat));
    
    for(j=0; j<order; j++){
        A[j]     = malloc(order*sizeof(**A));
        A_inv[j] = malloc(order*sizeof(**A_inv));
        A_co[j]  = malloc(order*sizeof(**A_co));
    }
    
    A[0][0] = 1;
    A[0][1] = 3;
    
    A[1][0] = 4;
    A[1][1] = 9;
    
    CoFactor(A, order, A_co);
    
    Transpose(A_co, order);
    
    double detA = 0;
    
    detA = Determinant(A, order);

    printf("\ndet A: %e", detA);

    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
            A_inv[j][k] = (1./detA)*A_co[j][k];
        
            // A_flat[count] = A[j][k];
            // count        += 1;
        }    
    }
    
    printf("\n\n A inverse.");
    printf("\n%e \t %e", A_inv[0][0], A_inv[0][1]);
    printf("\n%e \t %e", A_inv[1][0], A_inv[1][1]);
    
    double sum = 0.0;
        
    printf("\n\n");
    
    for(k=0; k<order; k++){  
        sum = 0.0;
      
        for(j=0; j<order; j++)  sum += A[k][j]*A_inv[j][k];
        
        printf("\nA*A_inv: %e", sum);
    }
    
    printf("\n\n");
    
    double** V;
    double*  W;
    
    W = malloc(order*sizeof(*W));
    V = malloc(order*sizeof(*V));
    
    for(j=0; j<order; j++)  V[j] = malloc(order*sizeof(**V));
    
    svdcmp(A, order, order, W, V);
    
    // for(j=0; j<order; j++)  printf("\n %e", W[j]);
    
    double Sum;
    
    printf("\n\nU:\n");
    for(i=0; i<order; i++){
        for(j=0; j<order; j++){
            printf("%e \t", A[i][j]);
        }
        
        printf("\n");
    }
    
    
    printf("\n\nA:\n");
    
    for(i=0; i<order; i++){
        for(j=0; j<order; j++){
            Sum = 0.0;
            
            for(k=0; k<order; k++)  Sum += V[j][k]*W[k]*A[i][k];
            
            printf("%e \t", Sum);
        }
        
        printf("\n");
    }
    
    
    printf("\n\nA inverse:\n");
    
    for(i=0; i<order; i++){
        for(j=0; j<order; j++){
            Sum = 0.0;
            
            for(k=0; k<order; k++)  Sum += V[i][k]*pow(W[k], -1.)*A[j][k];
            
            printf("%e \t", Sum);
        }
        
        printf("\n");
    }
    */
    
    
    int count   = 0;
    int order, order2;
    int errorHandler; 

    order  = 2*(kBinNumb-1);
    order2 = order*order;
    
    // length of N+1
    // double pivotArray[order+1];
    // double lapackWorkspace[order2];
    
    // printf("\n\nCovariance matrix.\n");
    // printMatrix(Covariance, order, order);
    
    /*
    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
            flatCov[count] = Covariance[j][k];
		  
   		    count += 1;
	    }
    }
    
    dgetrf_(&order, &order, flatCov, &order, pivotArray, &errorHandler);
    
    dgetri_(&order, flatCov, &order, pivotArray, lapackWorkspace, &order2, &errorHandler);
    
    count = 0;

    for(j=0; j<order; j++){
        for(k=0; k<order; k++){
            invCov[j][k] = flatCov[count];
            
           	count += 1; 
        }
    }
    
    // printf("\n\nInverse covariance matrix.\n");
    // printMatrix(invCov, order, order);
    
    double sum = 0.0;
        
    for(k=0; k<order; k++){  
        sum = 0.0;
      
        for(j=0; j<order; j++)  sum += Covariance[k][j]*invCov[j][k];
        
        printf("\nCov*invCov: %e", sum);
    }
    */
    
    double** U;
    double** V;
    double*  W;
    
    W = malloc(order*sizeof(*W));
    V = malloc(order*sizeof(*V));
    U = malloc(order*sizeof(*U));
    
    for(j=0; j<order; j++){
        V[j] = malloc(order*sizeof(**V));
        U[j] = malloc(order*sizeof(**U));
    }
    
    for(i=0; i<order; i++){
        for(j=0; j<order; j++){        
            U[i][j] = Covariance[i][j];
        }
    }
    
    svdcmp(U, order, order, W, V);

    printf("\n\nLargest singular values:\n");
    
    for(j=0; j<order; j++){  
        if(W[j] > 1.0)  printf("\n%e \t", W[j]);
    }
    
    for(i=0; i<order; i++){
        for(j=0; j<order; j++){
            invCov[i][j] = 0.0;
            
            for(k=0; k<order; k++){  
                if(W[k] > 1.0)  invCov[i][j] += V[i][k]*pow(W[k], -1.)*U[j][k];
            }
        }
    }
    
    // printf("\n\nInverse covariance from SVD:\n");
    // printMatrix(invCov, order, order);
    
    /*
    for(k=0; k<order; k++){  
        sum = 0.0;
      
        for(j=0; j<order; j++)  sum += Covariance[k][j]*invCov[j][k];
        
        printf("\nA*A_inv: %e", sum);
    }
    */
    
    sprintf(filepath, "%s/Data/Covariance/zCube_clipThreshold_1.0e+03_subVols_InverseCovariance.dat", root_dir);

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
