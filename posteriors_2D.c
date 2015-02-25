int calc_betaBiasPosterior(){
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
    double   biasInterval = (max_linearbias - min_linearbias)/dRes; 
        
    sprintf(filepath, "%s/Data/500s/betaBias_posterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            beta          = min_beta        + i*betaInterval;
            
            linearBias    = min_linearbias + j*biasInterval;
        
            fprintf(output, "%e \t %e \t %e \n", beta, linearBias, betaSigmaPosterior[i][j]);
        }
    }
    
    fclose(output);
    
    return 0;
}
