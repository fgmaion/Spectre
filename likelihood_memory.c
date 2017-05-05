int assign_chisq_kmaxes(){
  ChiSq_nkmaxes  = 8;

  ChiSq_ikmaxes  = (int *)     malloc(8*sizeof(int));
  ChiSq_kmaxes   = (double *)  malloc(8*sizeof(double));

  for(j=0; j<8; j++)  ChiSq_kmaxes[j] = 0.1*(1 + j);  
  
  // for(j=0; j<ChiSq_nkmaxes; j++)  printf("\nChi. sq. k_max: %.4lf", ChiSq_kmaxes[j]);  
    
  return 0;
}


int assign_LikelihoodMemory(){
  int kk, ll;

  int order = all_order; // 
  
  xdata             = malloc(order*sizeof(double));
  ydata             = malloc(order*sizeof(double));  // Decorrelated data.
  kdata             = malloc(order*sizeof(double));  // (kVals, kVals) 
  dkdata            = malloc(order*sizeof(double));  // 'Decorrelated k'
  
  xtheory           = (double ******) malloc(Res*sizeof(*xtheory));
  ytheory           = (double ******) malloc(Res*sizeof(*ytheory));    
  ChiSqGrid         = (double *****)  malloc(Res*sizeof(*ChiSqGrid));  // Parameter grid of evaluated Chi sq. values.
  
  for(j=0; j<Res; j++){
    xtheory[j] = (double *****) malloc(Res*sizeof(**xtheory));
    ytheory[j] = (double *****) malloc(Res*sizeof(**ytheory));
    ChiSqGrid[j] = (double ****) malloc(Res*sizeof(**ChiSqGrid));
    
    for(k=0; k<Res; k++){
      xtheory[j][k] = (double ****) malloc(Res*sizeof(***xtheory));
      ytheory[j][k] = (double ****) malloc(Res*sizeof(***ytheory));
      
      ChiSqGrid[j][k] = (double ***) malloc(Res*sizeof(***ChiSqGrid));
	  
      for(i=0; i<Res; i++){
        xtheory[j][k][i]  = (double ***) malloc(Res_ap*sizeof(****xtheory));
        ytheory[j][k][i]  = (double ***) malloc(Res_ap*sizeof(****ytheory));
        
        ChiSqGrid[j][k][i]  = (double **) malloc(Res_ap*sizeof(****ChiSqGrid));
	    
        for(ii=0; ii<Res_ap; ii++){
          xtheory[j][k][i][ii] = (double **) malloc(Res_ap*sizeof(*****xtheory));
          ytheory[j][k][i][ii] = (double **) malloc(Res_ap*sizeof(*****ytheory));
          
          ChiSqGrid[j][k][i][ii] = (double *) malloc(Res_ap*sizeof(*****ChiSqGrid));
          
          for(jj=0; jj<Res_ap; jj++){
            xtheory[j][k][i][ii][jj] = (double *) malloc(order*sizeof(******xtheory));
            ytheory[j][k][i][ii][jj] = (double *) malloc(order*sizeof(******ytheory));
            
            ChiSqGrid[j][k][i][ii][jj] = 0.0;
            
            for(kk=0; kk<order; kk++){
              xtheory[j][k][i][ii][jj][kk] = 0.0;
              ytheory[j][k][i][ii][jj][kk] = 0.0;
            }
          }
        }
      }  
    }
  }
    
  fsigma8Posterior = (double *)  malloc(Res*sizeof(*fsigma8Posterior));    
  sigmaPosterior   = (double *)  malloc(Res*sizeof(*sigmaPosterior));
  
  for(i=0; i<Res; i++){
    fsigma8Posterior[i] = 0.0;
    sigmaPosterior[i]   = 0.0;
  }
  
  return 0;
}


int assignCovMat(int mocks){
  // Multipoles is a [CatalogNumber][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin.
   Multipoles                = (double **) malloc(mocks*sizeof(*Multipoles));
  dMultipoles                = (double **) malloc(mocks*sizeof(*dMultipoles));  // Decorrelated multipole moments measurements
        
  for(j=0; j<mocks; j++){ 
     Multipoles[j]           = (double  *) malloc(order*sizeof( *Multipoles));
    dMultipoles[j]           = (double  *) malloc(order*sizeof(*dMultipoles));
  }
    
  kVals                      = (double  *) realloc(kVals, mono_order*sizeof(*kVals));

  MeanMultipoles             = (double  *) malloc(order*sizeof(*MeanMultipoles));    
    
  for(j=0; j<order; j++)  MeanMultipoles[j]   = 0.0;
    
  Covariance = gsl_matrix_alloc(order, order);
  sigma_norm = gsl_matrix_alloc(order, order);
    
  eval       = gsl_vector_alloc(order);
  evec       = gsl_matrix_alloc(order, order);
    
  w          = gsl_eigen_symmv_alloc(order); 
    
  col        = gsl_vector_alloc(order);  // Assign gsl_vector for eigenvector.
    
  return 0;
}


int prep_FFTlog_memory(){
  xi_mu_prefactor  = malloc(FFTlogRes*sizeof(double));  // Pre and Post factors for xi_mu FFT_log transform.
  xi_mu_postfactor = malloc(FFTlogRes*sizeof(double));  // Dropping e.g. transformOrder dependent terms.

  pk_mu_prefactor  = malloc(FFTlogRes*sizeof(double));
  pk_mu_postfactor = malloc(FFTlogRes*sizeof(double));

  FFTlog_Pk        = malloc(FFTlogRes*sizeof(double));  // Input P(k) interpolated to k's used by FFTlog calc.
  FFTlog_W0        = malloc(FFTlogRes*sizeof(double));  // W_0(r) evaluated at FFTlog rvals.
  FFTlog_W2        = malloc(FFTlogRes*sizeof(double));
  FFTlog_W4        = malloc(FFTlogRes*sizeof(double));
  FFTlog_W6        = malloc(FFTlogRes*sizeof(double));

  FFTlog_W0_joint  = malloc(FFTlogRes*sizeof(double)); // Joint W_0(r) evaluated at FFTlog rvals.
  FFTlog_W2_joint  = malloc(FFTlogRes*sizeof(double));
  FFTlog_W4_joint  = malloc(FFTlogRes*sizeof(double));
  FFTlog_W6_joint  = malloc(FFTlogRes*sizeof(double));

  FFTlog_Wk0       = malloc(FFTlogRes*sizeof(double)); // \tilde W_0(k) evaluated at FFTlog rvals.
  FFTlog_Wk2       = malloc(FFTlogRes*sizeof(double)); // \tilde W_2(k) evaluated at FFTlog rvals.
  
  return 0;
}
