int Calc_betaSigmaPosterior(){
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
        
        
    sprintf(filepath, "%s/Data/Posteriors/Multipoles_zCube_clipThreshold_1.0e+03_subVol_betaSigmaPosterior_lowRes_kmax_%.2f.dat", root_dir, ChiSqEval_kmax);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            beta          = min_beta + (max_beta - min_beta)*(i/dRes);
            
            velDispersion = min_velDisperse + (max_velDisperse - min_velDisperse)*(j/dRes);
        
            fprintf(output, "%e \t", betaSigmaPosterior[i][j]);
        }
        
        fprintf(output, "\n");
    }
    
    fclose(output);
    
    return 0;
}
