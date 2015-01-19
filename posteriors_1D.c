int calc_betaPosterior(){   
    PosteriorNorm = 0.0;
    
    double  betaInterval = (max_beta - min_beta)/dRes;
 
    for(k=0; k<Res; k++){
        betaPosterior[k] = 0.0;

        for(j=0; j<Res; j++){
            for(i=0; i<Res; i++){  
                betaPosterior[k] += exp(lnLikelihoodGrid[k][j][i]);
            }
        }
        
        if(betaPosterior[k] > PosteriorNorm){
            PosteriorNorm   = betaPosterior[k];
        }
    }

    for(k=0; k<Res; k++)  betaPosterior[k] /= PosteriorNorm;
    
    
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

/*
int Calc_betaError_95conf(int maxLikeIndex){
    double Prob          = 0.0;
    double totalProb     = 0.0;
    double confidence_95 = 0.0;
    
    for(i=0; i<Res; i++)  totalProb += betaPosterior[i];
    
    confidence_95        = 0.95*totalProb;
    
    Prob                += betaPosterior[maxLikeIndex];
    
    for(i=1; i<Res; i++){
        Prob            += betaPosterior[maxLikeIndex + i];
        
        Prob            += betaPosterior[maxLikeIndex - i];         
    
        if(Prob > confidence_95){
            break;
        }
    }
      
    printf("\n\nMax of beta posterior: %e", min_beta + (max_beta - min_beta)*(maxLikeIndex/dRes));
        
    printf("\n\n95 percent confidence on:  %f < beta < %f", min_beta + (max_beta - min_beta)*(maxLikeIndex - i)/dRes, min_beta + (max_beta - min_beta)*(maxLikeIndex + i)/dRes);
    
    return 0;
}


int maxLikelihoodBeta(){
    int    maxLikeIndex = 0;
    
    double maxLikelihoodBeta  = 0.0; 
    
    for(i=0; i<Res; i++){
        if(betaPosterior[i] > maxLikelihoodBeta){
            maxLikelihoodBeta = betaPosterior[i];
            
            maxLikeIndex      = i;
        }
    }
    
    return maxLikeIndex;
}*/
