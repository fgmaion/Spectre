int Calc_betaPosterior(){
    PosteriorNorm = 0.0;

    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            for(k=0; k<Res; k++){
                betaPosterior[i]   += exp(lnLikelihoodGrid[i][j][k]);
                
                if(betaPosterior[i] > PosteriorNorm){
                    PosteriorNorm   = betaPosterior[i];
                }
            }
        }
    }
    
    for(i=0; i<Res; i++)    betaPosterior[i] /= PosteriorNorm;
    
    sprintf(filepath, "%s/Data/Posteriors/Multipoles_%s_betaPosterior_Res_%d_kmax_%.2f_%s.dat", root_dir, surveyType, Res, ChiSqEval_kmax, theoryRSD_flag);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        beta = min_beta + (max_beta - min_beta)*(i/dRes);
        
        fprintf(output, "%e \t %e \n", beta, betaPosterior[i]);
    }
    
    fclose(output);
    
    int maxLikeIndex;
    
    maxLikeIndex = maxLikelihoodBeta();
    
    Calc_betaError_95conf(maxLikeIndex);
    
    return 0;
}


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
}

int Calc_sigmaPosterior(){
    // Velocity dispersion.
    PosteriorNorm = 0.0;

    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            for(k=0; k<Res; k++){
                sigmaPosterior[i] += exp(lnLikelihoodGrid[j][i][k]);
                
                if(sigmaPosterior[i] > PosteriorNorm){
                    PosteriorNorm   = sigmaPosterior[i];
                }
            }
        }
    }

    for(i=0; i<Res; i++)    sigmaPosterior[i] /= PosteriorNorm;
    
       
    sprintf(filepath, "%s/Data/Posteriors/Multipoles_zCube_clipThreshold_4.0e+00_subVol_sigmaPosterior_LowRes_kmax_%.2f.dat", root_dir, ChiSqEval_kmax);
    
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
        velDispersion = min_velDisperse + (max_velDisperse - min_velDisperse)*(i/dRes);
    
        fprintf(output, "%e \t %e \n", velDispersion, sigmaPosterior[i]);
    }
    
    fclose(output);

    return 0;
}


int Calc_A11SqPosterior(){
    return 0;
}
