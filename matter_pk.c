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
    if((lo_zlim==0.6) && (hi_zlim==0.9))  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.593_z_0.75.dat", root_dir);
    else                                  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.518_z_1.05.dat", root_dir);

    inputfile  = fopen(filepath, "r");
   
    line_count(inputfile, &pk_lineNo);
        
    sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));   // Interpolated theoretical P(k) on a regular grid in k.
    sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
    sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.

    for(j=0; j<pk_lineNo; j++)  fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);
    
    fclose(inputfile);
    
    spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        
    powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
    powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

    pt2Pk      = &splintLinearPk;
    
    // aexp       = 1./(1. + z_eff);  

    // P(z) = pow(b*sig8/sig8_fid, 2.)*P(0)
    // assuming linear p(k).
    
    if((lo_zlim==0.6) && (hi_zlim==0.9))  app_sigma8 = 0.593; 
    if((lo_zlim==0.9) && (hi_zlim==1.2))  app_sigma8 = 0.518;
    // app_sigma8 = 0.518;  // 0.9 < z < 1.2;  Will need changed.  

    // sigma8_calc();
                                                                                 // ln(a);    
    // printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 0.75, 0.593, f_Om_545(log(1./(1. + 0.75))), 0.593*f_Om_545(log(1./(1. + 0.75))));
    // printf("\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 1.05, 0.518, f_Om_545(log(1./(1. + 1.05))), 0.518*f_Om_545(log(1./(1. + 1.05))));

    return 0;
}


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


double sigma8_integrand(double k){
    return pow(2.*pi, -3.)*splintHODpk(k)*4.*pi*k*k*spherical_tophat(k, 8.)*spherical_tophat(k, 8.);
}


double sigma8_calc(){
    double Interim;
    
    Interim = qromb(&sigma8_integrand, 0.0000001, 200.);
    
    return sqrt(Interim);
}


int inputHODPk(){
    // sprintf(filepath, "%s/Data/500s/EisensteinHu_halofit_pk_%.2f.dat", root_dir, z_eff);    
    if((lo_zlim == 0.6) && (hi_zlim == 0.9))  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/nonlinear_matter_pk_sig8_0.593_z_0.75.dat", root_dir, z_eff);
    if((lo_zlim == 0.8) || (lo_zlim == 0.9))  sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/nonlinear_matter_pk_sig8_0.518_z_1.05.dat", root_dir, z_eff);
    
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

    pt2Pk      = &splintHODpk;
    
    aexp       = 1./(1. + z_eff);  

    // P(z) = pow(b*sig8/sig8_fid, 2.)*P(0)
    // assuming linear p(k).
    
    // Planck Cosmology
    // app_sigma8 = 0.59; // sigma8_calc();
                      
    // Value added mocks cosmology.   

    if((lo_zlim==0.6) && (hi_zlim==0.9))  app_sigma8 = 0.593;
    if((lo_zlim==0.8) || (hi_zlim==0.9))  app_sigma8 = 0.518;

    // ln(a);                                                                                                                         
    // printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 0.75, 0.593, f_Om_545(log(1./(1. + 0.75))), 0.593*f_Om_545(log(1./(1. + 0.75))));
    // printf("\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 1.05, 0.518, f_Om_545(log(1./(1. + 1.05))), 0.518*f_Om_545(log(1./(1. + 1.05))));
    
    return 0;
}


int print_fsigma8(){
  double redshift;
  double fD;

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
    fD = sig8_0*f_Om_545(log(1./(1. + redshift)))*linearGrowth_factor(log(1./(1. + redshift)));

    fprintf(output, "%.6lf \t %.6lf \n", redshift, fD);
  }
  
  fclose(output);

  return 0;
}


int print_matterpk(){
    sprintf(filepath, "%s/W1_Spectro_V5_0/camb/powerlawfits_linear_matter_pk_sig8_0.593_z_%.2f.dat", root_dir, z_eff);

    output = fopen(filepath, "w");
    
    for(j=0; j<1000; j++)  fprintf(output, "%e \t %e \n", j*0.1, splintHODpk(j*0.1));

    fclose(output);

    return 0;
}

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
