int set_oldestpriors(){
  min_bsigma8               =      0.05;
  max_bsigma8               =      1.00;

  min_fsigma8               =      0.00;
  max_fsigma8               =      0.80;

  min_velDisperse           =      0.00;
  max_velDisperse           =      7.00;             

  return 0;
}

int set_recordedpriors(){                                  // i.e. that on mull/stacpolly_backups/stacpolly_backup4/ driver_likelihood.c
  min_bsigma8               =      0.05;                  
  max_bsigma8               =      1.00;                  

  min_fsigma8               =      0.05;                  
  max_fsigma8               =      0.80;

  min_velDisperse           =      0.00;              
  max_velDisperse           =      6.00;                                

  return 0;
}

int set_normalpriors(){
  min_bsigma8               =      0.05;                  // FOR GRANETT 2D POSTERIOR:  0.2 < b \sig_8 < 1.6 
  max_bsigma8               =      1.20;                 

  min_fsigma8               =      0.05;                  // Priors on the model params.
  max_fsigma8               =      1.00;

  min_velDisperse           =      0.00;                  // CHANGED FROM 0.00 13/02/2017
  max_velDisperse           =     16.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)

  return 0;
}

int set_clippingpriors(){
  // clipping priors; 0.6 < z < 0.8
  min_bsigma8               =      0.50;                  // FOR GRANETT 2D POSTERIOR.                                                                                                      
  max_bsigma8               =      1.00;                  // Previously 0.2 < b \sig_8 < 1.6                                                                                                  

  min_fsigma8               =      0.05;                  // Priors on the model params.                                                                                                     
  max_fsigma8               =      0.80;                                                                                                                                                      

  min_velDisperse           =      0.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)
  max_velDisperse           =      6.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)

  return 0;
}

int set_widepriors(){                                                                                                                                                                       
  min_bsigma8               =      0.05;
  max_bsigma8               =      3.50;                                                                                                                                                                                                                                                                                                                                                  
  min_fsigma8               =      0.00;
  max_fsigma8               =      1.80;                                                                                                                                                      
  min_velDisperse           =      0.00;                                                                                                                                                     
  max_velDisperse           =     15.00;                                                                                                                                                     

  return 0;
}

int write_priors(){
  printf("\n\nAssumed uniform priors for %d^3 likelihood grid:", Res);

  printf("\n%.3lf < f * sigma_8 < %.3lf (%.3lf resolution)", min_fsigma8,     max_fsigma8,     fsigma8Interval);
  printf("\n%.3lf < b * sigma_8 < %.3lf (%.3lf resolution)", min_bsigma8,     max_bsigma8,     bsigma8Interval);
  printf("\n%.3lf <     sigma_p < %.3lf (%.3lf resolution)", min_velDisperse, max_velDisperse,   sigmaInterval);

  printf("\n");
  
  return 0;
}
