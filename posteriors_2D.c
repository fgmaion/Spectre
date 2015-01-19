int calc_betaSigmaPosterior(){
    PosteriorNorm = 0.0;

    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){        
            for(k=0; k<Res; k++){
                betaSigmaPosterior[i][j] += exp(lnLikelihoodGrid[i][j][k]);
     
                if(betaSigmaPosterior[i][j] > PosteriorNorm){
                    PosteriorNorm   = betaSigmaPosterior[i][j];
                }
            }
        }
    }
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){    
            betaSigmaPosterior[i][j] /= PosteriorNorm;
        }
    }
        
    double   betaInterval = (max_beta - min_beta)/dRes;
    double  sigmaInterval = (max_velDisperse - min_velDisperse)/dRes; 
        
    sprintf(filepath, "%s/Data/likelihood/clipping_mask500s_betaSigma_posterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            beta          = min_beta        + i*betaInterval;
            
            velDispersion = min_velDisperse + j*sigmaInterval;
        
            fprintf(output, "%e \t %e \t %e \n", beta, velDispersion, betaSigmaPosterior[i][j]);
        }
    }
    
    fclose(output);
    
    // 1D beta posterior
    PosteriorNorm = 0.0;
 
    for(i=0; i<Res; i++){
        betaPosterior[i] = 0.0;

        for(j=0; j<Res; j++)  betaPosterior[i] += betaSigmaPosterior[i][j];
        
        if(betaPosterior[i] > PosteriorNorm){
            PosteriorNorm   = betaPosterior[i];
        }
    }


    for(i=0; i<Res; i++)  betaPosterior[i] /= PosteriorNorm;
    
    
    sprintf(filepath, "%s/Data/likelihood/clipped_mask500s_betaPosterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(k=0; k<Res; k++){
        beta = min_beta + betaInterval*k;
        
        fprintf(output, "%e \t %e \n", beta, betaPosterior[k]);
    
        // fprintf(output, "%e \t %e \n", beta, exp(lnLikelihoodGrid[k][4][5]));
    }
    
    fclose(output);
    
    return 0;
}
