    // Inverse of M using LAPACK.

    double** M;
    
    M = (double **) malloc(3*sizeof(*M));
    
    for(j=0; j<3; j++)  M[j] = malloc(3*sizeof(**M));
    
    for(k=0; k<3; k++){
        for(j=0; j<3; j++){
            M[j][k] = 0.0;
        }
    }
    
    M[0][0] = 4.0;
    M[0][1] = 12.0;
    M[0][2] = -16.0;

    M[1][0] = 12.0;
    M[1][1] = 37.0;
    M[1][2] = -43.0;
    
    M[2][0] = -16.0;
    M[2][1] = -43.0;
    M[2][2] =  98.0;

    int    order  = 3;
    int    order2 = 9;

    // length of N+1
    double pivotArray[4];
    
    int    errorHandler; 

    double lapackWorkspace[9];

    double flattenedM[9];

    int count = 0;

    // Note: FORTRAN Ordering. M = {{a, b, c},
    //                              {d, e, f},
    //                              {g, h, i}};

    // Assigned via: M[0][0] = a, M[0][1] = b, M[0][2] = c, M[1][0] = d, etc.

    for(j=0; j<3; j++){
      for(k=0; k<3; k++){
	flattenedM[count] = M[j][k];
	count            += 1;
      }
    }
    
    // dgetrf_(&order, &order, flattenedM, &order, pivotArray, &errorHandler);
    
    // dgetri_(&order, flattenedM, &order, pivotArray, lapackWorkspace, &order2, &errorHandler);
    
    printf("\n %e \t %e \t %e", flattenedM[0], flattenedM[1], flattenedM[2]);
    printf("\n %e \t %e \t %e", flattenedM[3], flattenedM[4], flattenedM[5]);
    printf("\n %e \t %e \t %e", flattenedM[6], flattenedM[7], flattenedM[8]);
    
    
    // Cholesky decomposition of a Matrix, with LAPACK. 
    /*
    int INFO, size;
    char UPLO;
    
    UPLO ='L';

    size = 3;

    // Cholesky decomposition. Upper triangular matrix, lower triangular
    // left as original matrix. 

    dpotrf_(&UPLO, &size, flattenedM, &size, &INFO);
    
    count = 0;

    for(j=0; j<3; j++){
      printf("\n");
      for(k=0; k<3; k++){
	printf("%e \t", flattenedM[count]);
	count += 1;
      }
    }
