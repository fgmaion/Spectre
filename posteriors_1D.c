double calc_bsigma8Posterior(void){
  PosteriorNorm = 0.0;

  int     maxLikelihood_index = 0;

  double  bsigma8Interval = (max_bsigma8 - min_bsigma8)/dRes;

  double  bsigma8_posterior[Res];

  double  result = NAN;
  
  for(k=0; k<Res; k++){
    bsigma8_posterior[k] = 0.0;

    for(j=0; j<Res; j++){
      for(i=0; i<Res; i++)  bsigma8_posterior[k] += exp(-ChiSqGrid[j][k][i][0][0]/2.);
    }

    if(bsigma8_posterior[k] > PosteriorNorm){
      PosteriorNorm       = bsigma8_posterior[k];

      maxLikelihood_index = k;

      result              = min_bsigma8 + bsigma8Interval*maxLikelihood_index;
    }
  }

  // for(k=0; k<Res; k++)  bsigma8_posterior[k] /= PosteriorNorm;

  return result;
}


double calc_velDispPosterior(void){
  PosteriorNorm               = 0.0;

  int     maxLikelihood_index =   0;

  double     sigmaInterval = (max_velDisperse - min_velDisperse)/dRes;

  double sigma_posterior[Res];

  double  result = NAN;
  
  for(k=0; k<Res; k++){
    sigma_posterior[k] = 0.0;

    // fix last index, marginalise over first two - consistent with chiSq_minimisation.c
    for(j=0; j<Res; j++){
      for(i=0; i<Res; i++)  sigma_posterior[k] += exp(-ChiSqGrid[i][j][k][0][0]/2.);
    }

    if(sigma_posterior[k] > PosteriorNorm){
      PosteriorNorm       = sigma_posterior[k];

      maxLikelihood_index = k;

      result              = min_velDisperse + sigmaInterval*maxLikelihood_index;
    }
  }

  // for(k=0; k<Res; k++)  sigma_posterior[k] /= PosteriorNorm;

  return result;
}


double calc_fsigma8Posterior(void){ 
    PosteriorNorm = 0.0;
    
    int     maxLikelihood_index = 0;
    
    double  fsigma8Interval = (max_fsigma8 - min_fsigma8)/dRes;
    
    double  result = NAN; // Gotta gamble.
    
    for(k=0; k<Res; k++){
      fsigma8Posterior[k] = 0.0;

      for(j=0; j<Res; j++){
        for(i=0; i<Res; i++)  fsigma8Posterior[k] += exp(-ChiSqGrid[k][j][i][0][0]/2.);
      }

      if(fsigma8Posterior[k]  > PosteriorNorm){
        PosteriorNorm       = fsigma8Posterior[k];
        
        maxLikelihood_index = k;
        
        result              = min_fsigma8 + fsigma8Interval*maxLikelihood_index;
      }
    }

    // for(k=0; k<Res; k++)  fsigma8Posterior[k] /= PosteriorNorm;
    
    // for(k=0; k<Res; k++)  printf("\n%d \t %e \t %d", k, fsigma8Posterior[k], maxLikelihood_index);    
    
    // double blind_rand;
    // blind_rand = gsl_rng_uniform(gsl_ran_r) - 0.5;
    
    // sprintf(filepath, "%s/Data/500s/clipped_fsigma8Posterior_kmax_%.2f.dat", root_dir, ChiSq_kmax);
    /*
    sprintf(filepath, "%s/W1_Spectro_V5_0/fsigma8_posterior_kmax_%.2f.dat", root_dir, ChiSq_kmax);

    output = fopen(filepath, "w");
    
    for(k=0; k<Res; k++){
        fsigma8 = blind_rand + min_fsigma8 + fsigma8Interval*k;
        
        fprintf(output, "%e \t %e \t %e \n", fsigma8, fsigma8Posterior[k], exp(-pow((fsigma8 - 0.4322)/(0.08*sqrt(2.)), 2.)));
    }
    
    fclose(output);
    */
    
    // printf("\nmax likelihood fsig8: %le", result);
    
    // Calc_fsigma8_68conf(maxLikelihood_index, error);
    
    return result;
}

/*
int Calc_fsigma8_68conf(int maxLikeIndex, double* error){
    double Prob          = 0.0;
    double totalProb     = 0.0;
    double confidence_68 = 0.0;
    
    double     fsigma8Interval = (max_fsigma8 - min_fsigma8)/dRes;

    for(i=0; i<Res; i++)  totalProb += fsigma8Posterior[i];
    
    // by assumption, total probability in the range is 1. 
    confidence_68        = 0.68*totalProb;
    
    Prob                 = fsigma8Posterior[maxLikeIndex];
    
    // printf("\n\nTotal prob. %.4e, 68 limits %.4e", totalProb, confidence_68);
    
    for(i=1; i<Res/2; i++){
        Prob            += fsigma8Posterior[maxLikeIndex + i];
        
        Prob            += fsigma8Posterior[maxLikeIndex - i];         
    
        if(Prob>confidence_68){
            printf("\n\n\nMax of fsigma8 posterior: %e", min_fsigma8 + fsigma8Interval*maxLikeIndex);
        
            printf("\n68 percent confidence interval:  %f", i*fsigma8Interval);
            
            // printf("\n\n95 limits %.4e, prob: %.4e", confidence_95, Prob);
            
	    *error = i*fsigma8Interval;

            return 1;
        }
    }
    
    if(Prob < confidence_68)  printf("\n\n68 percent confidence limits lie outside fsigma8 prior range.");
    
    return 0;
}
*/
