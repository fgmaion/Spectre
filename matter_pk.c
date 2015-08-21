double powerlaw(double k){
  return pow(k, 4.00);
}


int inputLinearPk(){
    // sprintf(filepath, "%s/Data/500s/EisensteinHu_halofit_pk_%.2f.dat", root_dir, z_eff);
    sprintf(filepath, "%s/W1_Spectro_V5_0/camb/linear_matter_pk_sig8_0.593_z_%.2f.dat", root_dir, z_eff);
    
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

    
    pt2Pk      = &splintLinearPk;
    
    aexp       = 1./(1. + z_eff);  

    // P(z) = pow(b*sig8/sig8_fid, 2.)*P(0)
    // assuming linear p(k).
    app_sigma8 = 0.59; // sigma8_calc();
                                                                                 // ln(a);    
    printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", z_eff, app_sigma8, f_Om_545(log(1./(1.+z_eff))), app_sigma8*f_Om_545(log(1./(1.+z_eff))));
    
    // print_matterpk();
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
    //** Must match Cosmology declared in cosmology_planck2015.h, cosmology_valueaddedmocks.h **//
    //** This is NOT automatically ensured. **//
    
    // sprintf(filepath, "%s/Data/500s/EisensteinHu_halofit_pk_%.2f.dat", root_dir, z_eff);
    
    sprintf(filepath, "%s/W1_Spectro_V5_0/camb/nonlinear_matter_pk_sig8_0.593_z_%.2f.dat", root_dir, z_eff);

    //** HOD cube.  
    // sprintf(filepath, "%s/Data/HODCube/nonlinear_matter_pk_0.76.dat", root_dir, z_eff);
    
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
    
    // for(j=0; j<20; j++)  printf("\n%e \t %e", sdltk[j], sdltPk[j]);
    
    fclose(inputfile);
    
    
    spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        

    powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);
    
    powerlaw_regression(pk_lineNo,    8.,  10., 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

    
    pt2Pk      = &splintHODpk;
    
    aexp       = 1./(1. + z_eff);  

    // P(z) = pow(b*sig8/sig8_fid, 2.)*P(0)
    // assuming linear p(k).
    
    // Planck Cosmology
    // app_sigma8 = 0.59; // sigma8_calc();
                      
    // Value added mocks cosmology.   
    app_sigma8 = 0.557;
                                                                                 // ln(a);    
    printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", z_eff, app_sigma8, f_Om_545(log(1./(1.+z_eff))), app_sigma8*f_Om_545(log(1./(1.+z_eff))));

    printf("\n\nlinear growth rate: D_+(z_eff): %.6lf", linearGrowth_factor(log(1./(1. + z_eff))));
   
    // print_matterpk();
    
    print_fsigma8();

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
    return  splintHODpk(k)*exp(-0.5*pow(4.*k, 2.));
}
