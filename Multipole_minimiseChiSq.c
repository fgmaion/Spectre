int minimiseChiSq(){    
    minChiSq   = pow(10., 12.);
    
    // Numerical recipes/Fortran array indexing. 
    order      = hiMultipoleOrder*chiSq_kmaxIndex + 1;
    
    for(j=1; j<order; j++)  xdata[j]                  =   pow(sigmaNorm[j-1][j-1], -1.)*mvGauss[j-1];
    
    for(j=1; j<order; j++){
        ydata[j]    = 0.0;
        
        for(k=1; k<order; k++)    ydata[j]           += eigenVecs[k][j]*xdata[k];
    }
    
    printf("\n\nDecorrelated data.\n");

    for(j=1; j<order; j++)  printf("\n%d %e \t %e", j, xdata[j], ydata[j]);
    
    printf("\n\nbeta: %.2f, sigma: %.2f, A11 sq.: %.2f, kmax: %.2f", beta, velDispersion, A11Sq, ChiSqEval_kmax);
  
    printf("\n\nBeginning Chi sq. minimisation.");
  
    printf("\n\nPriors:");
    printf("\n%.2f < beta  < %.2f", min_beta, max_beta);  
    printf("\n%.2f < sigma < %.2f", min_velDisperse, max_velDisperse);
    printf("\n%.2f < A11 sq. < %.2f", min_A11Sq, max_A11Sq);
    printf("\n");
    
    for(jj=0; jj<Res; jj++){
        beta = min_beta + (max_beta - min_beta)*(jj/dRes);

        for(kk=0; kk<Res; kk++){
            velDispersion = min_velDisperse + (max_velDisperse - min_velDisperse)*(kk/dRes);

            for(ii=0; ii<Res; ii++){
                A11Sq = min_A11Sq + (max_A11Sq - min_A11Sq)*(ii/dRes);
      
                ChiSqGrid[jj][kk][ii]        = ChiSqEval(order);
                
                lnLikelihoodGrid[jj][kk][ii] = -0.5*ChiSqGrid[jj][kk][ii]; 

                if(ChiSqGrid[jj][kk][ii] < minChiSq){
                    minChiSq       = ChiSqGrid[jj][kk][ii];
                    
                    minChiSq_beta  = beta;
                    minChiSq_A11Sq = A11Sq;
                    minChiSq_sigma = velDispersion;
                    
                    printf("\n%e \t %e \t %e \t %e", beta, velDispersion, A11Sq, minChiSq);
                }
	       }
       }
    }
    
    printf("\n\nMinChiSq: %e", minChiSq);
    
    fprintfBestfit(order);
    
    return 0;
}


double ChiSqEval(int order){
    double ChiSq    = 0.0;
    
    double mono_j, quad_j;
    double mono_k, quad_k;
    
    // Numerical recipes/Fortran array indexing. 
    for(j=1; j<chiSq_kmaxIndex + 1; j++){
        xtheory[j]                   =  pow(sigmaNorm[j-1][j-1], -1.)*A11Sq*(*pt2Pk)(kMultipoles[j-1])*(*pt2RSD_k)(kMultipoles[j-1]*velDispersion, beta, 0); 
        xtheory[chiSq_kmaxIndex + j] =  pow(sigmaNorm[chiSq_kmaxIndex + j -1][chiSq_kmaxIndex + j -1], -1.)*A11Sq*(*pt2Pk)(kMultipoles[j-1])*(*pt2RSD_k)(kMultipoles[j-1]*velDispersion, beta, 1); 
    }
    
    for(j=1; j<order; j++){
        ytheory[j]  = 0.0; 
        
        for(k=1; k<order; k++)  ytheory[j]       += eigenVecs[k][j]*xtheory[k];
    }

    for(j=1; j<order; j++){  
        if(eigenVals[j] >= lowestKeptEigenvalue){
            ChiSq += pow(ydata[j] - ytheory[j], 2.)/eigenVals[j]; 
            
            // printf("\n %d \t %e \t %e \t %e \t %e", j, eigenVals[j], ydata[j], ytheory[j], ChiSq);
        }
    }
    
    // Given covariance, Maximum likelihood solution requires minimisation of chi sq. 

    return ChiSq;
}


int fprintfBestfit(int order){
    for(j=1; j<chiSq_kmaxIndex +1; j++){
        xtheory[j]                   = minChiSq_A11Sq*(*pt2Pk)(kMultipoles[j-1])*(*pt2RSD_k)(kMultipoles[j-1]*minChiSq_sigma, minChiSq_beta, 0);
        
        xtheory[chiSq_kmaxIndex + j] = minChiSq_A11Sq*(*pt2Pk)(kMultipoles[j-1])*(*pt2RSD_k)(kMultipoles[j-1]*minChiSq_sigma, minChiSq_beta, 1);
    }
    
    for(j=1; j<order; j++){
        ytheory[j]  = 0.0; 
        
        for(k=1; k<order; k++)  ytheory[j]       += eigenVecs[k][j]*xtheory[k];
    }
    
    
    sprintf(filepath, "%s/Data/Posteriors/bestfit_%s_beta_%.2f_sigma_%.2f_A11Sq_%.2f_kmax_%.2f_%s.dat", root_dir, surveyType, minChiSq_beta, minChiSq_sigma, minChiSq_A11Sq, ChiSqEval_kmax, theoryRSD_flag);

    output = fopen(filepath, "w");

    for(j=1; j<chiSq_kmaxIndex+1; j++)  fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \t %e\n", kMultipoles[j-1], mvGauss[j-1], mvGauss[chiSq_kmaxIndex + j - 1], xtheory[j], xtheory[chiSq_kmaxIndex + j], ydata[j], ydata[chiSq_kmaxIndex + j], ytheory[j], ytheory[chiSq_kmaxIndex + j], eigenVals[j], eigenVals[chiSq_kmaxIndex + j]);
    
    fclose(output);
    
    return 0;
}
