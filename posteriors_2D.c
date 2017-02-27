int calc_fsigma8_bsigma8_posterior(){
    PosteriorNorm = 0.0;
    
    int ll, mm;
    
    double** fsigma8_bsigma8_posterior;

    fsigma8_bsigma8_posterior = (double **) malloc(Res*sizeof(double*));

    for(j=0; j<Res; j++)    fsigma8_bsigma8_posterior[j]    = (double *)  malloc(Res*sizeof(double));
    
    for(i=0; i<Res; i++){
      for(j=0; j<Res; j++)  fsigma8_bsigma8_posterior[i][j] = 0.0;
    }

    
    for(jj=0; jj<Res; jj++){
      // fsigma8 = min_fsigma8 + fsigma8Interval*jj;
        for(kk=0; kk<Res; kk++){  
	  // bsigma8 = min_bsigma8 + bsigma8Interval*kk;
	  for(ii=0; ii<Res; ii++){
	    // velDispersion = min_velDisperse + sigmaInterval*ii;
	    for(ll=0;ll<Res_ap; ll++){
	      for(mm=0; mm<Res_ap; mm++){
		fsigma8_bsigma8_posterior[jj][kk]  +=  exp(lnLikelihoodGrid[jj][kk][ii][ll][mm]);
	      
                if(fsigma8_bsigma8_posterior[jj][kk] > PosteriorNorm){
                    PosteriorNorm   = fsigma8_bsigma8_posterior[jj][kk];
                }
	      }
	    }
          }
       }
    }
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){    
            fsigma8_bsigma8_posterior[i][j] /= PosteriorNorm;
        }
    }
        
    double     fsigma8Interval  = (max_fsigma8     - min_fsigma8)/dRes;
    double     bsigma8Interval  = (max_bsigma8     - min_bsigma8)/dRes;
    double       sigmaInterval  = (max_velDisperse - min_velDisperse)/dRes;
    double  alpha_padInterval   = (max_alpha_pad   - min_alpha_pad)/dRes_ap;
    double  epsilon_padInterval = (max_epsilon_pad - min_epsilon_pad)/dRes_ap;
    double       A11SqInterval  = (max_A11Sq       - min_A11Sq)/dRes;
        
    
    sprintf(filepath, "%s/W1_Spectro_V7_2/Granett/Nagoya_v7data_W%d_%.1lf_%.1lf_fsigma8_bsigma8_posterior.dat", root_dir, fieldFlag, lo_zlim, hi_zlim);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<Res; i++){
      for(j=0; j<Res; j++){
        fsigma8  = min_fsigma8 + i*fsigma8Interval;
            
        bsigma8  = min_bsigma8 + j*bsigma8Interval;
        
        fprintf(output, "%e \t %e \t %e \n", fsigma8, bsigma8, fsigma8_bsigma8_posterior[i][j]);
      }
    }
    
    fclose(output);
    

    printf("\n\nComplete.");
    
    return 0;
}

/*
int Pearson_rankCoefficient(){
  double     fsigma8Interval  = (max_fsigma8     - min_fsigma8)/dRes;                                                                                                                                                                      
  double     bsigma8Interval  = (max_bsigma8     - min_bsigma8)/dRes;   

  double r, n, sum_xi, sum_yi, sum_xiyi, sum_xi2, sum_yi2;

									
  for(jj=0; jj<Res; jj++){                                                                                                                                                                                                                 
    fsigma8 = min_fsigma8 + fsigma8Interval*jj;                                                                                                                                                                                         
    
    for(kk=0; kk<Res; kk++){                                                                                                                                                                                                             
      bsigma8 = min_bsigma8 + bsigma8Interval*kk;                                                                                                                                                                                     

      fsigma8_bsigma8_posterior[jj][kk];  


    }
  }

  return 0;
}
*/
