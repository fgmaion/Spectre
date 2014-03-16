int Calc_betaPosterior(){
    double* betaPosterior;
    double  PosteriorNorm = 0.0;

    betaPosterior = (double *)  malloc(Res*sizeof(*betaPosterior));    

    for(i=0; i<Res; i++){
        betaPosterior[i] = 0.0;

        for(j=0; j<Res; j++){
            for(k=0; k<Res; k++){
                betaPosterior[i] += ChiSqGrid[i][j][k];
            }
        }
    }
    
    for(i=0; i<Res; i++)    betaPosterior[i] /= Res*Res;
    
    for(i=0; i<Res; i++)    betaPosterior[i] = exp(-2.0*betaPosterior[i]);
       
    PosteriorNorm = SumDoubleArray(betaPosterior, Res);
    
    for(i=0; i<Res; i++)    betaPosterior[i] /= PosteriorNorm;
       
    sprintf(filepath, "%s/Data/Posteriors/test_betaPosterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        beta = 0.25 + 0.2*(i/dRes);
        
        fprintf(output, "%e \t %e \n", beta, betaPosterior[i]);
    }
    
    fclose(output);
    
    return 0;
}


int Calc_sigmaPosterior(){
    // Velocity dispersion.
    double* sigmaPosterior;
    double  PosteriorNorm = 0.0;

    sigmaPosterior = (double *)  malloc(Res*sizeof(*sigmaPosterior));    

    for(i=0; i<Res; i++){
        sigmaPosterior[i] = 0.0;

        for(j=0; j<Res; j++){
            for(k=0; k<Res; k++){
                sigmaPosterior[i] += ChiSqGrid[j][i][k];
            }
        }
    }
    
    for(i=0; i<Res; i++)    sigmaPosterior[i] /= Res*Res;
    
    for(i=0; i<Res; i++)    sigmaPosterior[i]  = exp(-2.0*sigmaPosterior[i]);
       
    PosteriorNorm = SumDoubleArray(sigmaPosterior, Res);
    
    for(i=0; i<Res; i++)    sigmaPosterior[i] /= PosteriorNorm;
       
    sprintf(filepath, "%s/Data/Posteriors/test_sigmaPosterior.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        velDispersion = 1.35 + 0.5*(i/dRes);
    
        fprintf(output, "%e \t %e \n", velDispersion, sigmaPosterior[i]);
    }
    
    fclose(output);

    return 0;
}


int Calc_A11SqPosterior(){

    return 0;
}
