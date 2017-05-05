double splintLinearPk(double k){
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);

  else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}


double linearPk_Gaussian(double k){
  return splintLinearPk(k)*exp(-0.5*pow(3.*k, 2.));
}


double splintHODpk(double k){
  // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec).
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);

  else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}


double sigma8_integrand(double k, void* params){
  return pow(2.*pi, -3.)*splintHODpk(k)*4.*pi*k*k*spherical_tophat(k, 8.)*spherical_tophat(k, 8.);
}


double sigma8_calc(){
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double result, error;

  gsl_function F;

  F.function = &sigma8_integrand;

  gsl_integration_qags(&F, 0.0, 200., 0, 1e-7, 1000, w, &result, &error);

  // double precision handles %.15lf for initialisation, e.g. p = PI.
  printf("\n\nsig_8: %.10lf (GSL)", sqrt(result));

  return sqrt(result);
}


double powerlaw(double k){
  return 2.0*pow(10., 4.)*pow(k/0.02, 4.00);
}

int input_powerlaw_n4(){
  if((lo_zlim==0.6) && (hi_zlim==0.9))  app_sigma8 = 0.593;
  if((lo_zlim==0.9) && (hi_zlim==1.2))  app_sigma8 = 0.518;

  pt2Pk      = &powerlaw;

  return 0;
}


int inputLinearPk(){
  double sig8;
  
  sprintf(model_flag, "linear");

  pt2Pk = &splintLinearPk;
  
  if((lo_zlim==0.6) && (hi_zlim<1.0))       sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.593_z_0.75.dat", root_dir);
  else if((lo_zlim==0.8) && (hi_zlim<1.3))  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.518_z_1.05.dat", root_dir);

  else{
    printf("\n\nError during input of real-space P(k) model.");

    exit(EXIT_FAILURE);
  }
  
  inputfile  = fopen(filepath, "r");
   
  line_count(inputfile, &pk_lineNo);
        
  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));   // Interpolated theoretical P(k) on a regular grid in k.
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
  sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.
  
  for(j=0; j<pk_lineNo; j++)  fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);
    
  fclose(inputfile);

  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        
  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for sigma_8 calc.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  sig8 = sigma8_calc();

  for(j=0; j<pk_lineNo; j++)  sdltPk[j] /= pow(sig8, 2.); // normalised to sigma_8 of unity. 
  /*
  for(j=0; j<pk_lineNo; j++){
    if((0.02 < sdltk[j]) && (sdltk[j] < 0.4))  printf("\n%.6le \t %.6le", sdltk[j], sdltPk[j]);
  }
  */
  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);
  
  // if((lo_zlim==0.6) && (hi_zlim<1.0))       app_sigma8 = 0.593139; 
  // else if((lo_zlim==0.8) && (hi_zlim<1.3))  app_sigma8 = 0.518062;  // CHANGED FROM LOZ=0.9 AND app_sigma8 = 0.518 on 22 Mar 2017.  
                                                                                
  // printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 0.75, 0.593, f_Om_545(log(1./(1. + 0.75))), 0.593*f_Om_545(log(1./(1. + 0.75))));
  // printf("\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 1.05, 0.518, f_Om_545(log(1./(1. + 1.05))), 0.518*f_Om_545(log(1./(1. + 1.05))));

  return 0;
}


int inputHODPk(){
  double sig8;

  sprintf(model_flag, "nonlinear");

  pt2Pk = &splintHODpk;
  
  if((lo_zlim == 0.6) && (hi_zlim < 1.0))                           sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/nonlinear_matter_pk_sig8_0.593_z_0.75.dat", root_dir);
  else if((lo_zlim == 0.8) || (lo_zlim == 0.9) && (lo_zlim < 1.3))  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/nonlinear_matter_pk_sig8_0.518_z_1.05.dat", root_dir);
  
  else{
    printf("\n\nError for input of real-space P(k) model.");

    exit(EXIT_FAILURE);
  }

  printf("\n\nLoading %s", filepath);
  
  inputfile  = fopen(filepath, "r");
    
  line_count(inputfile, &pk_lineNo);
        
  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
  sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.

  for(j=0; j<pk_lineNo; j++)  fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);
    
  fclose(inputfile);
    
  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        
  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);  // Add power laws for k<<1 and k>>1 for FFTlog calcs.    
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  sig8 = sigma8_calc();

  sdltPk[j] /= pow(sig8, 2.); // normalised to unit sigma_8.
  
  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);
  
  // aexp       = 1./(1. + z_eff);  
  
  // P(z) = pow(b*sig8/sig8_fid, 2.)*P(0)
  // assuming linear p(k).
    
  // Value added mocks cosmology.   
  // if((lo_zlim==0.6) && (hi_zlim<1.0))       app_sigma8 = 0.638897; // CHANGED FROM app_sigma8 = 0.593 ON 22 Mar 17;  
  // else if((lo_zlim==0.8) && (hi_zlim<1.3))  app_sigma8 = 0.550554; // CHANGED FROM app_sigma8 = 0.518 ON 22 Mar 17;
  
  // ln(a);                                                                                                                         
  // printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 0.75, 0.593, f_Om_545(log(1./(1. + 0.75))), 0.593*f_Om_545(log(1./(1. + 0.75))));
  // printf("\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 1.05, 0.518, f_Om_545(log(1./(1. + 1.05))), 0.518*f_Om_545(log(1./(1. + 1.05))));

  // for(j=0; j<pk_lineNo; j++)  printf("\n %.6le \t %.6le", sdltk[j], splintHODpk(sdltk[j] + 0.01));
  
  return 0;
}


int print_fsigma8(){
  double redshift;
  double fD;

  void* null_point = NULL;
  
  // Based on: 
  // D+(0.75)   = 0.680582
  // sig8(0.75) = 0.5570   (from CAMB) 
  // sig8(z)    = D+ sig8(0.0)
  // 
  // -> sig8(0) = 0.8184 

  double sig8_0 = 0.8184;

  sprintf(filepath, "%s/W1_Spectro_V7_0/f_sig8.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<100; j++){
    redshift = j*0.01;

    // sig8(0) x f(z) x D+(z).
    fD = sig8_0*f_Om_545(log(1./(1. + redshift)), null_point)*linearGrowth_factor(log(1./(1. + redshift)));

    fprintf(output, "%.6lf \t %.6lf \n", redshift, fD);
  }
  
  fclose(output);

  return 0;
}

/*
int print_matterpk(){
    sprintf(filepath, "%s/W1_Spectro_V5_0/camb/powerlawfits_linear_matter_pk_sig8_0.593_z_%.2f.dat", root_dir, z_eff);

    output = fopen(filepath, "w");
    
    for(j=0; j<1000; j++)  fprintf(output, "%e \t %e \n", j*0.1, splintHODpk(j*0.1));

    fclose(output);

    return 0;
}*/


double HODPk_Gaussian(double k){
    return  splintHODpk(k)*exp(-0.5*pow(3.*k, 2.));
}


double splint_maskedRSD_pk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);
    
    else if(k<0.001)   return pk_loA*pow(k, 3. + pk_lon); 

    else{
        double Interim;
    
        splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);
    
        return Interim;
    }
}


int set_maskedRSDpaper_pk(){
    sprintf(filepath, "%s/Data/HODTheoryPk/outdated_cosmology/cambExtendedPk_hod_20.00.dat", root_dir);
    
    inputfile  = fopen(filepath, "r");
    
    ch         = 0;
    pk_lineNo  = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  pk_lineNo += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    // Interpolated theoretical P(k) on a regular grid in k. 
    sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
    sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
    
    // Second derivates of HOD P(k) for cubic spline.           
    sdlt2d = realloc(sdlt2d,         pk_lineNo*sizeof(*sdlt2d));

    for(j=0; j<pk_lineNo; j++)       fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);
    
    fclose(inputfile);
    
    spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        
    powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);
    
    powerlaw_regression(pk_lineNo,    8.,  10., 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

    pt2Pk = &splint_maskedRSD_pk; 

    return 0;
}
