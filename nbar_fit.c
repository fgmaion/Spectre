double CSR(double z){
  return 0.5*(1. - erf(17.465*(0.424 - z)));
}

double model_NzGaussian(double z, double lnA, double z0, double sigma){
  double A = exp(lnA);

  return A*CSR(z)*exp(-pow((z - z0)/sigma, 2.));
}

double leastSquares(const gsl_vector *v, void *params){
  double  lnA, z0, sigma, z_minfit; 

  double  chiSq = 0.0;

  double* p = (double *) params;

  lnA     = gsl_vector_get(v, 0);
  z0      = gsl_vector_get(v, 1);
  sigma   = gsl_vector_get(v, 2);

  if(lo_zlim < 0.80)  z_minfit = 0.5; // lo-z slics                                                                                                                                                               
  if(lo_zlim > 0.85)  z_minfit = 0.8;

  // printf("\n%e \t %e \t %e", lnA, z0, sigma);

  for(j=0; j<chibin_no; j++){
    if((Nchi[j] > 0.0) && (zbins[j] > z_minfit)){ 
      Interim = model_NzGaussian(zbins[j], lnA, z0, sigma) - Nchi[j];

      // printf("\n%e \t %e \t %e", zbins[j], Nchi[j], model_NzGaussian(0.7, 8.1, 0.65, 0.05));

      chiSq  += pow(Interim, 2.0);
    }
  }

  return chiSq;
}


int fitted_nbar_calc(){
  // equal intervals in comoving distance, for both W1 and W4 fields.                                                                                                                                            
  prep_nbar();

  double chi;

  int global_fieldFlag;

  double fit_params[3];

  // Store initial fieldFlag.                                                                                                                                                                                    
  global_fieldFlag = fieldFlag;

  // Two field average.                                                                                                                                                                                          
  for(ii=1; ii<5; ii=ii+3){
    fieldFlag = ii;

    // analysis on **mock** catalogues.                                                                                                                                                                          
    if(data_mock_flag == 0){
      if(loopCount<10)        sprintf(filepath, "%s/mocks_v1.7/W%d/mock_00%d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
      else if(loopCount<100)  sprintf(filepath, "%s/mocks_v1.7/W%d/mock_0%d_VAC_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);
      else                    sprintf(filepath, "%s/mocks_v1.7/W%d/mock_%d_VAC_Nagoya_v6_Samhain.dat",    vipersHOD_dir, fieldFlag, loopCount);

      CatalogueInput_500s(filepath);
      
      gal_z = &zobs[0];

      spec_weights();
      
      assignAcceptance_true();
    }

    /*
    // analysis on **data** catalogues.                                                                                                                                                                          
    if(data_mock_flag == 1){
      // W1 catalogue.                                                                                                                                                                                       
      // Also assigns sampling: ESR = TSR x SSR.                                                                                                                                                             
      DataInput();

      // do not impose redshift cuts just yet, obtain full n(z). Some redshifts >2, remove these objects.                                                                                                    
      for(j=0; j<Vipers_Num; j++){
	if((0.1<zobs[j]) && (zobs[j]<1.7)) Acceptanceflag[j] = true;

	else{Acceptanceflag[j] = false;}
      }
    }
    */

    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){
	chi                 = interp_comovingDistance(zobs[j]);

	Index               = (int) floor(chi/chi_interval);

	chibins[Index]     += chi/sampling[j];

	Nchi[Index]        +=  1./sampling[j];
      }
    }
  }

  for(j=0; j<chibin_no; j++){
      chibins[j]  /= Nchi[j];

      if(Nchi[j] == 0) chibins[j] = (j + 0.5)*chi_interval;
  }

  for(j=0; j<chibin_no; j++)  comovVol[j] = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;


  Nz_fit(fit_params);


  for(j=0; j<chibin_no; j++)     nbar[j]  = model_NzGaussian(zbins[j], fit_params[0], fit_params[1], fit_params[2])/comovVol[j];

  if(data_mock_flag == 0)        sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/nbar_fitted/nbar_fitted_Nagoya_v7_Samhain_mock_%d_twofield_avg_%.1lf_%.1lf.dat", root_dir, loopCount, lo_zlim, hi_zlim);
    
  // if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/nbar_150/nbar_smooth_%.1lf_Nagoya_v7_Samhain_twofield_avg.dat", root_dir, kernel_width);

  double z_minfit;

  if(lo_zlim < 0.80)  z_minfit = 0.5; // lo-z slics                                                                                                                                                              
  if(lo_zlim > 0.85)  z_minfit = 0.8;

  output = fopen(filepath, "w");

  for(j=0; j<chibin_no; j++){  
    if(zbins[j] < z_minfit){
      fprintf(output, "%e \t %e \n", chibins[j], 0.0);
    }

    else{fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);}
  }

  fclose(output);

  // Restore global field Flag.                                                                                                                                                                                  
  fieldFlag = global_fieldFlag;

  return 0;
}


int Nz_fit(double output[]){
  // nmsimplex2 -> O(N) rather than O(N^2)                                                                                                                                                                        
  const gsl_multimin_fminimizer_type*  T = gsl_multimin_fminimizer_nmsimplex2;

  // unnecessary.
  double par[3] = {1.0, 2.0, 3.0};

  gsl_multimin_fminimizer*             s = NULL;

  gsl_vector *ss, *x;

  gsl_multimin_function minex_func;

  int      status;  
  double     size;

  size_t iter = 0;
    
  // Starting point in parameter space, x = (A, z0, alpha, beta)                                                                                                                                                  
  x  = gsl_vector_alloc(3);                                                                                                                                                                                         
  // Starting parameters                                                                                                                                                                                         
  gsl_vector_set(x, 0,     8.100);  // ln(A)                                                                                                                                                                      
  gsl_vector_set(x, 1,     0.700);  // z0                                                                                                                                                                        
  gsl_vector_set(x, 2,     0.300);  // sigma                                                                                                                                                                     
  
  ss = gsl_vector_alloc(3);                                                                                                                                                                                     
  
  gsl_vector_set_all(ss, 1.0);   
  
  minex_func.n      =              3;
  minex_func.f      =   leastSquares;
  minex_func.params =            par;

  s = gsl_multimin_fminimizer_alloc(T, 3);
  
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
  
  do{
      iter++;

      status = gsl_multimin_fminimizer_iterate(s);
      
      if(status)
	break;

      size   = gsl_multimin_fminimizer_size(s);

      status = gsl_multimin_test_size(size, 1e-2);

  } while((status == GSL_CONTINUE) && (iter < 100));
  
  double lnA   = gsl_vector_get(s->x, 0);
  double z0    = gsl_vector_get(s->x, 1);
  double sigma = gsl_vector_get(s->x, 2);

  output[0] = lnA;
  output[1] =  z0;
  output[2] = sigma;

  gsl_multimin_fminimizer_free(s);
  
  gsl_vector_free(x);

  gsl_vector_free(ss);
  
  return 0;
}
