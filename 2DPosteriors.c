int Calc_betaSigmaPosterior(){
    double** betaSigmaPosterior;
    
    double   maxPosterior  = 0.0;
    double   PosteriorNorm = 0.0;
    
    betaSigmaPosterior = (double **) malloc(Res*sizeof(*betaSigmaPosterior));
    
    for(j=0; j<Res; j++)  betaSigmaPosterior[j] = malloc(Res*sizeof(**betaSigmaPosterior));

    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            betaSigmaPosterior[i][j] = 0.0;
        
            for(k=0; k<Res; k++){
                betaSigmaPosterior[i][j] += ChiSqGrid[i][j][k];
            }
        }
    }
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            betaSigmaPosterior[i][j] /= Res;
            
            betaSigmaPosterior[i][j]  = exp(-2.0*betaSigmaPosterior[i][j]);
     
            PosteriorNorm            += betaSigmaPosterior[i][j];
            
            if(betaSigmaPosterior[i][j] > maxPosterior){
                maxPosterior = betaSigmaPosterior[i][j];
            }
        }
    }
           
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            betaSigmaPosterior[i][j] /= maxPosterior;
        }
    }
    
    sprintf(filepath, "%s/Data/Posteriors/test_betaSigmaPosterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            beta          = 0.25 + 0.2*(i/dRes);
            
            velDispersion = 1.35 + 0.5*(j/dRes);
        
            fprintf(output, "%e \t", betaSigmaPosterior[i][j]);
        }
        
        fprintf(output, "\n");
    }
    
    fclose(output);
    
    return 0;
}
