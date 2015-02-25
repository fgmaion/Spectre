int calc_fsigma8Posterior(){   
    PosteriorNorm = 0.0;
    
    int     maxLikelihood_index;
    
    double     fsigma8Interval = (max_fsigma8     - min_fsigma8)/dRes;
 
    for(k=0; k<Res; k++){
        fsigma8Posterior[k] = 0.0;

        for(j=0; j<Res; j++){
            for(i=0; i<Res; i++){      
                fsigma8Posterior[k] += exp(lnLikelihoodGrid[k][j][i]);
            }
        }
        
        if(fsigma8Posterior[k] > PosteriorNorm){
            PosteriorNorm       = fsigma8Posterior[k];
            
            maxLikelihood_index = k;
        }
    }

    for(k=0; k<Res; k++)  fsigma8Posterior[k] /= PosteriorNorm;
    
    
    sprintf(filepath, "%s/Data/500s/clipped_fsigma8Posterior_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    
    output = fopen(filepath, "w");
    
    for(k=0; k<Res; k++){
        fsigma8 = min_fsigma8 + fsigma8Interval*k;
        
        fprintf(output, "%e \t %e \t %e \n", fsigma8, fsigma8Posterior[k], exp(-pow((fsigma8 - 0.4322)/(0.08*sqrt(2.)), 2.)));
    }
    
    fclose(output);
    
    Calc_fsigma8_95conf(maxLikelihood_index);
    
    return 0;
}


int Calc_fsigma8_95conf(int maxLikeIndex){
    double Prob          = 0.0;
    double totalProb     = 0.0;
    double confidence_95 = 0.0;
    
    double     fsigma8Interval = (max_fsigma8     - min_fsigma8)/dRes;

    for(i=0; i<Res; i++)  totalProb += fsigma8Posterior[i];
    
    // by assumption, total probability in the range is 1. 
    confidence_95        = 0.95*totalProb;
    
    Prob                 = fsigma8Posterior[maxLikeIndex];
    
    // printf("\n\nTotal prob. %.4e, 95 limits %.4e", totalProb, confidence_95);
    
    for(i=1; i<Res/2; i++){
        Prob            += fsigma8Posterior[maxLikeIndex + i];
        
        Prob            += fsigma8Posterior[maxLikeIndex - i];         
    
        if(Prob>confidence_95){
            printf("\n\n\nMax of fsigma8 posterior: %e", min_fsigma8 + fsigma8Interval*maxLikeIndex);
        
            printf("\n95 percent confidence interval:  %f", i*fsigma8Interval);
            
            // printf("\n\n95 limits %.4e, prob: %.4e", confidence_95, Prob);
            
            break;
        }
    }
    
    if(Prob < confidence_95)  printf("\n\n95 percent confidence limits lie outside fsigma8 prior range.");
    
    return 0;
}
