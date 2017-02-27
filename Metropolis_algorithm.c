double likelihood(double alpha, double epsilon, double fsig8, double bsig8, double sigp){
  double exp_alpha = 1.0000, exp_epsilon = 0.0000, exp_fsig8 = 0.5,  exp_bsig8 = 1.2, exp_sigp = 5.0;
  double sig_alpha = 0.0001, sig_epsilon = 0.0001, sig_fsig8 = 0.05, sig_bsig8 = 0.3, sig_sigp = 2.0;

  // assumes multivariate Gaussian, with diagonal covariance. 
  return exp(-0.5*pow((alpha - exp_alpha)/sig_alpha, 2.))*exp(-0.5*pow((epsilon - exp_epsilon)/sig_epsilon, 2.))*exp(-0.5*pow((fsig8 - exp_fsig8)/sig_fsig8, 2.))*exp(-0.5*pow((bsig8 - exp_bsig8)/sig_bsig8, 2.))*exp(-0.5*pow((sigp - exp_sigp)/sig_sigp, 2.));
}

double vipers_chiSq(){
  // return -2.*log(likelihood(alpha_pad, epsilon_pad, fsigma8, bsigma8, velDispersion));
  return ChiSqEval_ap();
}

int new_position(){
  /*
  double alpha_scale   =  0.0001;
  double epsilon_scale =  0.0001;

  double fsig8_scale   =  0.05;  // approx. 10% of expectation                                                                                                                                                                               
  double bsig8_scale   =  0.30;
  double sigp_scale    =  2.00;
  */
  
  double alpha_scale   =  0.01; 
  double epsilon_scale =  0.01;

  double fsig8_scale   =  0.05;  // approx. 10% of expectation
  double bsig8_scale   =  0.05;
  double sigp_scale    =  0.50;
  
  // alpha_pad           += rand_gaussian(gsl_ran_r,   alpha_scale); // (0.5 - gsl_rng_uniform(gsl_ran_r))*alpha_scale;
  epsilon_pad         += rand_gaussian(gsl_ran_r, epsilon_scale); // (0.5 - gsl_rng_uniform(gsl_ran_r))*epsilon_scale;

  fsigma8             += rand_gaussian(gsl_ran_r,   fsig8_scale); // (0.5 - gsl_rng_uniform(gsl_ran_r))*fsig8_scale;
  bsigma8             += rand_gaussian(gsl_ran_r,   bsig8_scale); // (0.5 - gsl_rng_uniform(gsl_ran_r))*bsig8_scale;
 
  velDispersion       += rand_gaussian(gsl_ran_r,    sigp_scale); //(0.5 - gsl_rng_uniform(gsl_ran_r))*sigp_scale;

  return 0;
}


int proposal(double* chiSq_val){
  // Update position in parameter space according to Metropolis algorithm, pg. 14 of Stats by Heavens.
  double theta_alpha      =        alpha_pad;
  double theta_epsilon    =      epsilon_pad;
  double theta_fsig8      =          fsigma8;
  double theta_bsig8      =          bsigma8;
  double theta_sigp       =    velDispersion;

  double theta_chiSq, thetaStar_chiSq;

  double likelihood_ratio =      0.0;

  double draw             =      0.0;

  // evaluate at current position.   
  theta_chiSq             =     *chiSq_val;

  // assign new position. 
  new_position();

  thetaStar_chiSq         = vipers_chiSq();

                           // likelihood at theta*                                  
  likelihood_ratio        = exp(-0.5*(thetaStar_chiSq - theta_chiSq));

  draw                    = gsl_rng_uniform(gsl_ran_r);  

  if(draw < likelihood_ratio){
    // accept new position.  Already assigned. 
    *chiSq_val       = thetaStar_chiSq;
    
    return 0;
  }
    
  else{
    // reject; restore theta.
    alpha_pad        =      theta_alpha;
    epsilon_pad      =    theta_epsilon; 

    fsigma8          =      theta_fsig8;
    bsigma8          =      theta_bsig8; 
    velDispersion    =       theta_sigp; 
      
    *chiSq_val       =      theta_chiSq; 

    return 0;
    
  }
}


int metropolis_mcmc(int choose_mock){
  gsl_rng_env_setup();  // Random variable generation.

  gsl_ran_T = gsl_rng_default;

  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

  printf("\n\nBeginning MCMC.");
  
  // MCMC with symmetric proposal function. 

  // To be recovered. 
  // double exp_alpha = 0.02, exp_epsilon = 0.03, exp_fsig8 = 0.5, exp_bsig8 = 1.2, exp_sigp = 5.0;
  // double sig_alpha = 0.02, sig_epsilon = 0.01, sig_fsig8 = 0.2, sig_bsig8 = 0.3, sig_sigp = 2.0;
  
  // starting values.
  alpha_pad      =    1.0;
  epsilon_pad    =    0.0;

  fsigma8        =   0.50;
  bsigma8        =   1.00;
  velDispersion  =   6.00;

  double chiSq_val =  0.0;

  int jjj, fsig8_index, eps_index;                                                                                                                                                                     

  int N0         = 1000000;

  /*                                                                                                                                                                                                          
  double                      Counts[100][100];                                                                                                                                                                                                                                                                                                                                                                                     
  double fsig8_interval =  (0.70  - 0.20)/100.;                                                                                                                                                                   
  double   eps_interval =  (0.20  + 0.20)/100.;  
  
  for(j=0; j<100; j++){                                                                                                                                                                                          
    for(i=0; i<100; i++)  Counts[i][j] = 0.0;                                                                                                                                                                     
  } 
  */
  // Initialise global variables used by ChiSqEval().
  prepChiSq_minimisation();
  
  // Model signal + hypothetical errors.                                                                                                                                                                          
  // load_aptest_multipoles();

  // Errors 2x better than VIPERS (diagonal covariance). \Chi^2 set to diagonal only.                                                                                                                           
  // scale_Cov();

  FFTlogRes         = 768;

  pt2maskMultipoles = &splint_VIPERS_maskMultipoles;

  // Must have mask multipoles available, uncomment prep_VIPERS_maskMultipoles().                                                                                                                                                            
  FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);

  kvals_matchup();
   
  if(data_mock_flag == 0)  load_mock(choose_mock);
  if(data_mock_flag == 1)  load_data();
  
  // set_meanmultipoles();

  /*
  // all input loaded, destroy lockfile.                                                                                                                                                                                       
  char*       lockfile_path;                                                                                                                                                                                                               
  char delete_lockfile[200];                                                                                                                                                                                                               
                                                                                                                                                                                                                                           
  lockfile_path = malloc(200*sizeof(*lockfile_path));                                                                                                                                                                                      
                                                                                                                                                                                                                                           
  lockfile_path = getenv("LOCKFILEDIR");                                                                                                                                                                                                   
                                                                                                                                                                                                                                           
  sprintf(delete_lockfile, "rm -rf %s", lockfile_path);                                                                                                                                                                                    
                                                                                                                                                                                                                                           
  // remove lockfile. Executes delete_lockfile string as terminal command.                                                                                                                                                                 
  system(delete_lockfile);                                                                                                                                                                                                                 
 
  system("echo $LOCKFILEDIR");
                                                                                                                                                                                                                                          
  system("echo 'lockfile destroyed.'");                                                                                                                                                                                                    
  */
  // new zero mean, unit variance, decorrelated variables.                                                                                                                                                       
  for(j=0; j<order; j++){
    ydata[j] = 0.0;

    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++)  ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
  }

  minChiSq = pow(10., 12.);

  prep_dlnPR_dlnk();

  chiSq_val = vipers_chiSq();
    
  printf("\n\nChi sq.  of start: %le", chiSq_val);

  // Write file to show progress/not hanging at lock file. 
  // if(data_mock_flag == 0)   sprintf(filepath, "%s/W1_Spectro_V7_2/apchain_kmax_%.2lf/vipers_%d_mock_%d.dat",              root_dir, ChiSq_kmax, N0, choose_mock);
  // if(data_mock_flag == 1)   sprintf(filepath, "%s/W1_Spectro_V7_2/apchain_kmax_%.2lf/vipers_%d_data_W%d_%.1lf_%.1lf.dat", root_dir, ChiSq_kmax, N0, fieldFlag, lo_zlim, hi_zlim);
  
  printf("\n\nLikelihood set, beginning to jump.");
  
  if(data_mock_flag == 0)   sprintf(filepath, "%s/W1_Spectro_V7_2/mock_kmax_%.1lf_chains/vipers_%d_mock_%d_W%d_%.1lf_%.1lf.dat", root_dir, ChiSq_kmax, N0, choose_mock, fieldFlag, lo_zlim, hi_zlim);                                   
  if(data_mock_flag == 1)   sprintf(filepath, "%s/W1_Spectro_V7_2/data_kmax_%.1lf_chains/vipers_%d_data_W%d_%.1lf_%.1lf.dat", root_dir, ChiSq_kmax, N0, fieldFlag, lo_zlim, hi_zlim);        

  output = fopen(filepath, "w");

  fclose(output);
  
  //int ij, jk, loop, res;

  //res  = (int) floor(N0/50);
  
  //loop = 0;
  
  // burn in;
  for(jjj=0; jjj<5000; jjj++)  proposal(&chiSq_val);    
  
  // double chain_store[1000][5];

  // printf("\n\n");
  
  double** array;

  array = (double **) malloc(N0*sizeof(double*));
  
  for(j=0; j<N0; j++)  array[j] = (double *) malloc(6*sizeof(double));

  for(jjj=0; jjj<N0; jjj++){
    proposal(&chiSq_val);

    // fsig8_index = (int) floor( (fsigma8     - 0.20)/fsig8_interval);      
    // eps_index   = (int) floor( (epsilon_pad + 0.20)/eps_interval);

    // printf("\n%d \t %d", fsig8_index, eps_index);

    // if((fsig8_index>0) && (fsig8_index<100) && (eps_index>0) && (eps_index<100))  Counts[fsig8_index][eps_index] += 1;
    
    // chain_store[jjj%1000][0] =     alpha_pad;
    // chain_store[jjj%1000][1] =   epsilon_pad;
    // chain_store[jjj%1000][2] =       fsigma8;
    // chain_store[jjj%1000][3] =       bsigma8;
    // chain_store[jjj%1000][4] = velDispersion;
 
    //if(jjj%1000 == 0){
      // output = fopen(filepath, "a");
      
      //for(ij=0; ij<100; ij++){
	//for(jk=0; jk<100; jk++) fprintf(output, "%le \t", Counts[ij][jk]);

	//fprintf(output, "\n");
      //}

      // for(j=0; j<1000; j++)  fprintf(output, "%.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf", chain_store[j][0], chain_store[j][1], chain_store[j][2], chain_store[j][3], chain_store[j][4]);

      // fclose(output);

      //      for(j=0; j<1000; j++)  
      // printf("%.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf \n", chain_store[jjj][0], chain_store[jjj][1], chain_store[jjj][2], chain_store[jjj][3], chain_store[jjj][4]);
      //}
    
    array[jjj][0] =     alpha_pad;
    array[jjj][1] =   epsilon_pad;
    array[jjj][2] =       fsigma8;
    array[jjj][3] =       bsigma8;
    array[jjj][4] = velDispersion;
    array[jjj][5] =     chiSq_val;

    // printf("%.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf \n", array[jjj][0], array[jjj][1], array[jjj][2], array[jjj][3], array[jjj][4]);

    // loop += 1;
  }
  
  output = fopen(filepath, "w");

  for(j=0; j<N0; j++) fprintf(output, "%.3lf \t %.3lf \t %.3lf \t %.3lf \t %.3lf \n", array[j][0], array[j][1], array[j][2], array[j][3], array[j][4]);

  fclose(output);
  
  //output = fopen(filepath, "w");//for(ij=0; ij<100; ij++){
  //  for(jk=0; jk<100; jk++) fprintf(output, "%le \t", Counts[ij][jk]);

  //  fprintf(output, "\n");
  //}

  //fclose(output);
  
  return 0;
}
