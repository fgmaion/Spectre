double splintHODpk(double k){
  // Interpolated matter power spectrum evaluated at |k_vec|.
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);
  else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}

double HODPk_Gaussian(double k){
  return  splintHODpk(k)*exp(-0.5*pow(3.*k, 2.));
}

int nonlinear_pk(){
  char buffer[200];

  get_linearsig8();

  pt2Pk = &splintHODpk;
  
  sprintf(model_flag, "nonlinear");
  
  sprintf(filepath,   "/home/mjw/CAMB/models/camb_matter_pk_linearity_1_Om_cdm_%.3lf_Om_v_%.3lf_H0_%.2lf_z_%.3lf.dat", Om_m - Om_b, Om_v, 100.*h, z_eff);
  // sprintf(filepath,   "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/pkmodels/nonlinear_matter_pk_sig8_0.593_z_0.75.dat");
  
  printf("\n\nLoading %s", filepath);

  while((inputfile = fopen(filepath, "r")) == NULL){  // Generate file if it can't be found. 
    camb_call(1, z_eff);      // nonlinear/linear flag, redshift;
  }

  linecount_header(inputfile, 1, &pk_lineNo);
  // line_count(inputfile, &pk_lineNo);
  
  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
  sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.

  fgets(buffer, 200, inputfile);                        // Skip header after rewind. 
  
  for(j=0; j<pk_lineNo; j++){
    fscanf(inputfile, "    %lE    %lE \n", &sdltk[j], &sdltPk[j]);
    // printf("\n%le \t %le", sdltk[j], sdltPk[j]);
  }
  
  fclose(inputfile);
    
  // spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        
  // powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);  // Add power laws for k<<1 and k>>1 for FFTlog calcs.    
  // powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  // camb_sig8 = 0.593;
  // camb_sig8 = sigma8_calc();    // Linear sigma_8 is required.

  for(j=0; j<pk_lineNo; j++)  sdltPk[j] /= pow(camb_sig8, 2.); // normalised to unit sigma_8.
  
  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  // for(j=0; j<pk_lineNo; j++)  printf("\n %.6le \t %.6le", sdltk[j], sdltPk[j]);

  print_fsig8();  // needs camb sig_8(z).
  
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

