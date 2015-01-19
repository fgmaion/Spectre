int Covariance_eigenVecs(int mocks){        
    // This function computes the eigenvalues of the real symmetric matrix Covariance.  The diagonal and lower triangular part of Covariance are destroyed during the computation.
    gsl_eigen_symmv(Covariance, eval, evec, w);
    
    gsl_eigen_symmv_free(w);
    
    // Descending sort in magnitude. eigenvalues should be non-negative definite.
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_DESC);
    
    for(i=0; i<order; i++){
        double eval_i          = gsl_vector_get(eval, i);

        gsl_vector_view evec_i = gsl_matrix_column(evec, i);
        
        printf("eigenvalue  = %.3e \n", eval_i);

        // printf("eigenvalue  = %.3e, eigenvector = \n", eval_i);
        
        // gsl_vector_fprintf(stdout, &evec_i.vector, "%.3e");
    }
    
    // The diagonal and lower triangular part of covariance are destroyed during the eigenvalue/vector calculation. reassign.
    initialiseCovariance(mocks);
    
    // prewhiten(mocks);
    /*
    // get first eigenvector from evec matrix. eigenvectors are normalised. 
    gsl_matrix_get_col(col, evec, 2);
    
    // first eigenvalue. 
    double eval_0 = gsl_vector_get(eval, 2);
    
    // Covariance x eigenvector_0 = lambda_0 eigenvector_0.
    for(j=0; j<order; j++){
        Interim = 0.0;
      
        for(k=0; k<order; k++)  Interim += gsl_matrix_get(Covariance, j, k)*gsl_vector_get(col, k);  
        
        // printf("\n%e \t %e", Interim, gsl_vector_get(col, j)*eval_0);
    }*/
    
    /*
    double eval_j;
    
    for(j=0; j<order; j++){      
        // jth component of eval which is a vector of the eigenvalues. 
        eval_j = gsl_vector_get(eval, j);
    
        if((smallestEigenvalue >= eval_j) && (eval_j >= 0.5)){
          smallestEigenvalue  = eval_j;
        }
    
        if(eval_j<0.0){
            printf("\nError: covariance is symmetric and non-negative definite. eigenvalues should be real and positive. This is not the case");
        
            return 1;
        }
    }*/
    
    // Store correlated, zero mean, unit variance variables. 
    for(k=0; k<order; k++){
      for(i=0; i<mocks; i++)  Multipoles[i][k] = dMultipoles[i][k]*gsl_matrix_get(sigma_norm, k, k);
    }
    
    printf("\n");
    
    // Decorrelated multipole measurements. 
    for(j=0; j<order; j++){  
        // get eigenvector j from evec matrix. 
        gsl_matrix_get_col(col, evec, j);
    
        for(i=0; i<mocks; i++){            
            dMultipoles[i][j] = 0.0;
        
            for(k=0; k<order; k++){
                dMultipoles[i][j] += gsl_vector_get(col, k)*Multipoles[i][k];
            }
        }
    }
    
    printf("Decorrelated multipoles calculated.\n");
    
    // Reform covariance of decorrelated, zero mean, unit variance variables. Should be diagonal, with non-negative definite eigenvalues corresponding to eval. 
    
    // For running minimisation, the following commands should be commented.  
    // initialiseCovariance(mocks);
    
    // fprintf_dCov();
    
    return 0;
}


int fprintf_dCov(){
    sprintf(filepath, "%s/Data/likelihood/dClipped_mask500s_Covariance.dat", root_dir);

    output = fopen(filepath, "w"); 

    for(j=0; j<order; j++){
        for(k=0; k<order; k++){        
            fprintf(output, "%e \t", gsl_matrix_get(Covariance, j, k));                  
        }
    
        fprintf(output, "\n");
    }

    fclose(output);

    return 0;
}
