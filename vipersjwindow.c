double splint_VIPERS_jmaskMono(double r){
  if(r<0.15)   return 1.;
  if(r>900.)   return 0.;

  //if(r>4000.)  return 0.;

  if(r<jhiRes_hihiRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hihi, VIPERS_jmaskMono_hihi, VIPERS_jmaskMono2D_hihi, VIPERS_jmask_lineNo_hihi, r, &Interim);

    // amplitude determined by powerlaw_regression at r=1. -> A parameter.                                                                                                                                     
    return Interim;
  }

  else if(r<jloRes_highRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hi, VIPERS_jmaskMono_hi, VIPERS_jmaskMono2D_hi, VIPERS_jmask_lineNo_hi, r, &Interim);

    // amplitude determined by powerlaw_regression at r=1. -> A parameter.                                                                                                                                     
    return Interim;
  }

  else{
    double Interim;

    splint(VIPERS_jmaskr_lo, VIPERS_jmaskMono_lo, VIPERS_jmaskMono2D_lo, VIPERS_jmask_lineNo_lo, r, &Interim);

    // amplitude determined by powerlaw_regression at r=1. -> A parameter.                                                                                                                                     
    return Interim;
  }
}

double splint_VIPERS_jmaskQuad(double r){
  if(r<0.2)    return 0.;
  if(r>900.0)  return 0.;

  //if(r>4000.0) return 0.;

  if(r<jhiRes_hihiRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hihi, VIPERS_jmaskQuad_hihi, VIPERS_jmaskQuad2D_hihi, VIPERS_jmask_lineNo_hihi, r, &Interim);

    return Interim;
  }

  else if(r<jloRes_highRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hi, VIPERS_jmaskQuad_hi, VIPERS_jmaskQuad2D_hi, VIPERS_jmask_lineNo_hi, r, &Interim);

    return Interim;
  }

  else{
    double Interim;

    splint(VIPERS_jmaskr_lo, VIPERS_jmaskQuad_lo, VIPERS_jmaskQuad2D_lo, VIPERS_jmask_lineNo_lo, r, &Interim);

    return Interim;
  }
}

double splint_VIPERS_jmaskHex(double r){
  if(r<0.2)    return 0.;
  if(r>4000.0) return 0.;

  if(r<jhiRes_hihiRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hihi, VIPERS_jmaskHex_hihi, VIPERS_jmaskHex2D_hihi, VIPERS_jmask_lineNo_hihi, r, &Interim);

    return Interim;
  }

  else if(r<jloRes_highRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hi, VIPERS_jmaskHex_hi, VIPERS_jmaskHex2D_hi, VIPERS_jmask_lineNo_hi, r, &Interim);

    return Interim;
  }

  else{
    double Interim;

    splint(VIPERS_jmaskr_lo, VIPERS_jmaskHex_lo, VIPERS_jmaskHex2D_lo, VIPERS_jmask_lineNo_lo, r, &Interim);

    return Interim;
  }
}

double splint_VIPERS_jmaskOct(double r){
  if(r<0.2)    return 0.;
  if(r>4000.0) return 0.;

  if(r<jhiRes_hihiRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hihi, VIPERS_jmaskOct_hihi, VIPERS_jmaskOct2D_hihi, VIPERS_jmask_lineNo_hihi, r, &Interim);

    return Interim;
  }

  else if(r<jloRes_highRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hi, VIPERS_jmaskOct_hi, VIPERS_jmaskOct2D_hi, VIPERS_jmask_lineNo_hi, r, &Interim);

    return Interim;
  }

  else{
    double Interim;

    splint(VIPERS_jmaskr_lo, VIPERS_jmaskOct_lo, VIPERS_jmaskOct2D_lo, VIPERS_jmask_lineNo_lo, r, &Interim);

    return Interim;
  }
}

double splint_VIPERS_jmaskDec(double r){
  if(r<0.2)    return 0.;
  if(r>4000.0) return 0.;

  if(r<jhiRes_hihiRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hihi, VIPERS_jmaskDec_hihi, VIPERS_jmaskDec2D_hihi, VIPERS_jmask_lineNo_hihi, r, &Interim);

    return Interim;
  }

  else if(r<jloRes_highRes_join){
    double Interim;

    splint(VIPERS_jmaskr_hi, VIPERS_jmaskDec_hi, VIPERS_jmaskDec2D_hi, VIPERS_jmask_lineNo_hi, r, &Interim);

    return Interim;
  }

  else{
    double Interim;

    splint(VIPERS_jmaskr_lo, VIPERS_jmaskDec_lo, VIPERS_jmaskDec2D_lo, VIPERS_jmask_lineNo_lo, r, &Interim);

    return Interim;
  }
}

double splint_VIPERS_jmaskMultipoles(double r, int transformOrder){
  switch(transformOrder){
    case 0:
      return  splint_VIPERS_jmaskMono(r);
    case 2:
      return  splint_VIPERS_jmaskQuad(r);
    case 4:
      return  splint_VIPERS_jmaskHex(r);
  }
}

int print_jwindowCorrfn(){
  double r;
  char   windowxi_filepath[200];

  double maxlog   = log10(4000.);
  double zerolog  = log10(0.001);
  double logbinsz = log10(1.050);

  int    nbins    = (int) ceil((maxlog - zerolog)/logbinsz);
  
  
  sprintf(windowxi_filepath, "%s_mixedRes_hex.dat", filepath);

  output = fopen(windowxi_filepath, "w");

  for(j=0; j<nbins; j++){
    r = pow(10., zerolog + j*logbinsz);
    
    fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", r,  splint_VIPERS_jmaskMono(r), splint_VIPERS_jmaskQuad(r), splint_VIPERS_jmaskHex(r), splint_VIPERS_jmaskOct(r), splint_VIPERS_jmaskDec(r));
  }

  fclose(output);

  return 0;
}

int prep_VIPERS_jmaskMultipoles(){
  char           header[200];
  
  char   loRes_filepath[200];
  char   hiRes_filepath[200];
  char hihiRes_filepath[200];

  if(mull==0){
    sprintf(filepath, "%s/Qmultipoles/Ql_W1W4_Nag_v7_specweight_nbar_Pfkp_8000_%.1lf_%.1lf_thread_1", maskmultipoles_path, lo_zlim, hi_zlim);
  }

  else if(mull == 1){
    sprintf(filepath, "%s/Qmultipoles/maskmultipoles_W1W4_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf", maskmultipoles_path, lo_zlim, hi_zlim);
  }

  else{
    printf("Error onjoint Q-multipoles load.");
  }
  
  printf("\nand joint ... %s", filepath);
  
  sprintf(hihiRes_filepath, "%s_hihiRes_hex.dat", filepath);
  sprintf(  hiRes_filepath, "%s_hiRes_hex.dat",   filepath);
  sprintf(  loRes_filepath, "%s_loRes_hex.dat",   filepath);    
  
  // Super high resolution: very large number density of pairs counting for r ~ 1. 
  inputfile = fopen(hihiRes_filepath, "r");

  line_count(inputfile, &VIPERS_jmask_lineNo_hihi);
  // linecount_header(inputfile, 1, &VIPERS_jmask_lineNo_hihi);
  
  // Memory allocation
  VIPERS_jmaskr_hihi        = realloc( VIPERS_jmaskr_hihi,      VIPERS_jmask_lineNo_hihi*sizeof(double));
    
  VIPERS_jmaskMono_hihi     = realloc( VIPERS_jmaskMono_hihi,   VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskQuad_hihi     = realloc( VIPERS_jmaskQuad_hihi,   VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskHex_hihi      = realloc( VIPERS_jmaskHex_hihi,    VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskOct_hihi      = realloc( VIPERS_jmaskOct_hihi,    VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskDec_hihi      = realloc( VIPERS_jmaskDec_hihi,    VIPERS_jmask_lineNo_hihi*sizeof(double));
    
  VIPERS_jmaskMono2D_hihi   = realloc( VIPERS_jmaskMono2D_hihi, VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskQuad2D_hihi   = realloc( VIPERS_jmaskQuad2D_hihi, VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskHex2D_hihi    = realloc( VIPERS_jmaskHex2D_hihi,  VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskOct2D_hihi    = realloc( VIPERS_jmaskOct2D_hihi,  VIPERS_jmask_lineNo_hihi*sizeof(double));
  VIPERS_jmaskDec2D_hihi    = realloc( VIPERS_jmaskDec2D_hihi,  VIPERS_jmask_lineNo_hihi*sizeof(double));

  fgets(header, 200, inputfile); // skip header
    
  for(j=0; j<VIPERS_jmask_lineNo_hihi; j++){  
    fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le \n", &VIPERS_jmaskr_hihi[j], &VIPERS_jmaskMono_hihi[j], &VIPERS_jmaskQuad_hihi[j],  &VIPERS_jmaskHex_hihi[j], 
                                                                                                                         &VIPERS_jmaskOct_hihi[j],  &VIPERS_jmaskDec_hihi[j]);
  }

  fclose(inputfile);
    
  // Divide out volume factor, r^2 for vol of shell at r, r for dlogr width of shell -> r^3
  for(j=0; j<VIPERS_jmask_lineNo_hihi; j++){
    // printf("%le \t %le \t %le \t %le \n", VIPERS_maskr_hihi[j], VIPERS_maskMono_hihi[j], VIPERS_maskQuad_hihi[j], VIPERS_maskHex_hihi[j]);
    VIPERS_jmaskMono_hihi[j] /= pow(VIPERS_jmaskr_hihi[j], 3.); 
    VIPERS_jmaskQuad_hihi[j] /= pow(VIPERS_jmaskr_hihi[j], 3.); 
    VIPERS_jmaskHex_hihi[j]  /= pow(VIPERS_jmaskr_hihi[j], 3.);
    VIPERS_jmaskOct_hihi[j]  /= pow(VIPERS_jmaskr_hihi[j], 3.);
    VIPERS_jmaskDec_hihi[j]  /= pow(VIPERS_jmaskr_hihi[j], 3.);
  }
  
  // Determine amplitude of pair counts for r << 1 where the correlation function is flat. 
  flatSlope_amp(VIPERS_jmask_lineNo_hihi, 0.1, 0.2, 1.0, VIPERS_jmaskr_hihi, VIPERS_jmaskMono_hihi, &jmask_monopolenorm_hihi);
    
  // Normalise pair counts such that the amplitude of the monopole tends to 1 for r<<1.
  for(j=0; j<VIPERS_jmask_lineNo_hihi; j++){
    VIPERS_jmaskMono_hihi[j] /= jmask_monopolenorm_hihi;
    VIPERS_jmaskQuad_hihi[j] /= jmask_monopolenorm_hihi;
    VIPERS_jmaskHex_hihi[j]  /= jmask_monopolenorm_hihi;
    VIPERS_jmaskOct_hihi[j]  /= jmask_monopolenorm_hihi;
    VIPERS_jmaskDec_hihi[j]  /= jmask_monopolenorm_hihi;
  }
    
  // Spline for interpolation.
  spline(VIPERS_jmaskr_hihi, VIPERS_jmaskMono_hihi, VIPERS_jmask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_jmaskMono2D_hihi);
  spline(VIPERS_jmaskr_hihi, VIPERS_jmaskQuad_hihi, VIPERS_jmask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_jmaskQuad2D_hihi);
  spline(VIPERS_jmaskr_hihi, VIPERS_jmaskHex_hihi,  VIPERS_jmask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_jmaskHex2D_hihi);
  spline(VIPERS_jmaskr_hihi, VIPERS_jmaskOct_hihi,  VIPERS_jmask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_jmaskOct2D_hihi);
  spline(VIPERS_jmaskr_hihi, VIPERS_jmaskDec_hihi,  VIPERS_jmask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_jmaskDec2D_hihi);
  
  // Outdated.  Set the r value at which the super high resolution points will join the high resolution counts.  
  // jhiRes_hihiRes_join = VIPERS_jmaskr_hihi[VIPERS_jmask_lineNo_hihi - 1];
  jhiRes_hihiRes_join = 1.8;
  
  // High resolution
  inputfile = fopen(hiRes_filepath, "r");

  // line_count(inputfile, &VIPERS_jmask_lineNo_hi);
  linecount_header(inputfile, 1, &VIPERS_jmask_lineNo_hi);
    
  // Memory allocation
  VIPERS_jmaskr_hi        = realloc( VIPERS_jmaskr_hi, VIPERS_jmask_lineNo_hi*sizeof(double));
    
  VIPERS_jmaskMono_hi     = realloc( VIPERS_jmaskMono_hi,   VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskQuad_hi     = realloc( VIPERS_jmaskQuad_hi,   VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskHex_hi      = realloc( VIPERS_jmaskHex_hi,    VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskOct_hi      = realloc( VIPERS_jmaskOct_hi,    VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskDec_hi      = realloc( VIPERS_jmaskDec_hi,    VIPERS_jmask_lineNo_hi*sizeof(double));
    
  VIPERS_jmaskMono2D_hi   = realloc( VIPERS_jmaskMono2D_hi, VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskQuad2D_hi   = realloc( VIPERS_jmaskQuad2D_hi, VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskHex2D_hi    = realloc( VIPERS_jmaskHex2D_hi,  VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskOct2D_hi    = realloc( VIPERS_jmaskOct2D_hi,  VIPERS_jmask_lineNo_hi*sizeof(double));
  VIPERS_jmaskDec2D_hi    = realloc( VIPERS_jmaskDec2D_hi,  VIPERS_jmask_lineNo_hi*sizeof(double));
  
  fgets(header, 200, inputfile); // skip header
    
  for(j=0; j<VIPERS_jmask_lineNo_hi; j++)  fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le \n", &VIPERS_jmaskr_hi[j],   &VIPERS_jmaskMono_hi[j], &VIPERS_jmaskQuad_hi[j], 
  						                                                        &VIPERS_jmaskHex_hi[j], &VIPERS_jmaskOct_hi[j],    &VIPERS_jmaskDec_hi[j]);
 
  fclose(inputfile);
  
  // Divide out volume factor, r^2 for vol of shell at r, r for dlogr width of shell -> r^3
  for(j=0; j<VIPERS_jmask_lineNo_hi; j++){
    VIPERS_jmaskMono_hi[j] /= pow(VIPERS_jmaskr_hi[j], 3.);
    VIPERS_jmaskQuad_hi[j] /= pow(VIPERS_jmaskr_hi[j], 3.);
    VIPERS_jmaskHex_hi[j]  /= pow(VIPERS_jmaskr_hi[j], 3.);
    VIPERS_jmaskOct_hi[j]  /= pow(VIPERS_jmaskr_hi[j], 3.);
    VIPERS_jmaskDec_hi[j]  /= pow(VIPERS_jmaskr_hi[j], 3.);
  }
  
  lowRes_amplitudeCalc(VIPERS_jmask_lineNo_hihi, VIPERS_jmask_lineNo_hi, 1.0, 1.5, VIPERS_jmaskr_hihi, VIPERS_jmaskr_hi, VIPERS_jmaskMono_hihi, VIPERS_jmaskMono_hi, &jmask_monopolenorm_hi);
        
  // printf("\n\n%e", jmask_monopolenorm_hi);

  for(j=0; j<VIPERS_jmask_lineNo_hi; j++){
    VIPERS_jmaskMono_hi[j] *= jmask_monopolenorm_hi;
    VIPERS_jmaskQuad_hi[j] *= jmask_monopolenorm_hi;
    VIPERS_jmaskHex_hi[j]  *= jmask_monopolenorm_hi;
    VIPERS_jmaskOct_hi[j]  *= jmask_monopolenorm_hi;
    VIPERS_jmaskDec_hi[j]  *= jmask_monopolenorm_hi;
  }
    
  // Spline for interpolation.
  spline(VIPERS_jmaskr_hi, VIPERS_jmaskMono_hi, VIPERS_jmask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_jmaskMono2D_hi);
  spline(VIPERS_jmaskr_hi, VIPERS_jmaskQuad_hi, VIPERS_jmask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_jmaskQuad2D_hi);
  spline(VIPERS_jmaskr_hi, VIPERS_jmaskHex_hi,  VIPERS_jmask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_jmaskHex2D_hi);
  spline(VIPERS_jmaskr_hi, VIPERS_jmaskOct_hi,  VIPERS_jmask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_jmaskOct2D_hi);
  spline(VIPERS_jmaskr_hi, VIPERS_jmaskDec_hi,  VIPERS_jmask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_jmaskDec2D_hi);
    
  // Outdated. Set the r value at which the  high resolution counts will join the low resolution counts.  
  // jloRes_highRes_join = VIPERS_jmaskr_hi[VIPERS_jmask_lineNo_hi - 1];
  jloRes_highRes_join = 8.5;
  
                                                //--//
                                                
  // Low resolution on larger scales.
  inputfile = fopen(loRes_filepath, "r");

  // line_count(inputfile, &VIPERS_jmask_lineNo_lo);
  linecount_header(inputfile, 1, &VIPERS_jmask_lineNo_lo);
    
  // Memory allocation
  VIPERS_jmaskr_lo        = realloc( VIPERS_jmaskr_lo,      VIPERS_jmask_lineNo_lo*sizeof(double));
    
  VIPERS_jmaskMono_lo     = realloc( VIPERS_jmaskMono_lo,   VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskQuad_lo     = realloc( VIPERS_jmaskQuad_lo,   VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskHex_lo      = realloc( VIPERS_jmaskHex_lo,    VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskOct_lo      = realloc( VIPERS_jmaskOct_lo,    VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskDec_lo      = realloc( VIPERS_jmaskDec_lo,    VIPERS_jmask_lineNo_lo*sizeof(double));
  
  VIPERS_jmaskMono2D_lo   = realloc( VIPERS_jmaskMono2D_lo, VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskQuad2D_lo   = realloc( VIPERS_jmaskQuad2D_lo, VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskHex2D_lo    = realloc( VIPERS_jmaskHex2D_lo,  VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskOct2D_lo    = realloc( VIPERS_jmaskOct2D_lo,  VIPERS_jmask_lineNo_lo*sizeof(double));
  VIPERS_jmaskDec2D_lo    = realloc( VIPERS_jmaskDec2D_lo,  VIPERS_jmask_lineNo_lo*sizeof(double));
  
  fgets(header, 200, inputfile); // skip header
    
  for(j=0; j<VIPERS_jmask_lineNo_lo; j++){  fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le \n", &VIPERS_jmaskr_lo[j],   &VIPERS_jmaskMono_lo[j], &VIPERS_jmaskQuad_lo[j], 
                                                                                                         &VIPERS_jmaskHex_lo[j], &VIPERS_jmaskOct_lo[j],  &VIPERS_jmaskDec_lo[j]);
  }
  
  fclose(inputfile);
  
  // Divide out volume factor, r^2 for vol of shell at r, r for dlogr width of shell -> r^3
  for(j=0; j<VIPERS_jmask_lineNo_lo; j++){
    // printf("\n%e \t %e \t %e \t %e", VIPERS_maskr_lo[j], VIPERS_maskMono_lo[j], VIPERS_maskQuad_lo[j], VIPERS_maskHex_lo[j]);    
    VIPERS_jmaskMono_lo[j] /= pow(VIPERS_jmaskr_lo[j], 3.);
    VIPERS_jmaskQuad_lo[j] /= pow(VIPERS_jmaskr_lo[j], 3.);
    VIPERS_jmaskHex_lo[j]  /= pow(VIPERS_jmaskr_lo[j], 3.);
    VIPERS_jmaskOct_lo[j]  /= pow(VIPERS_jmaskr_lo[j], 3.);
    VIPERS_jmaskDec_lo[j]  /= pow(VIPERS_jmaskr_lo[j], 3.);
  }
    
  // Calculate amplitude factor to join lo res counts to hihi res, which are properly normalised at this point. 
  // Previously commented 19 NOV: 
  lowRes_amplitudeCalc(VIPERS_jmask_lineNo_hi, VIPERS_jmask_lineNo_lo, 8.,   9., VIPERS_jmaskr_hi, VIPERS_jmaskr_lo, VIPERS_jmaskMono_hi, VIPERS_jmaskMono_lo, &jmask_monopolenorm_lo);
        
  for(j=0; j<VIPERS_jmask_lineNo_lo; j++){
    VIPERS_jmaskMono_lo[j] *= jmask_monopolenorm_lo;
    VIPERS_jmaskQuad_lo[j] *= jmask_monopolenorm_lo;
    VIPERS_jmaskHex_lo[j]  *= jmask_monopolenorm_lo;
    VIPERS_jmaskOct_lo[j]  *= jmask_monopolenorm_lo;
    VIPERS_jmaskDec_lo[j]  *= jmask_monopolenorm_lo;
  }
  
  // Spline for interpolation.
  spline(VIPERS_jmaskr_lo, VIPERS_jmaskMono_lo, VIPERS_jmask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_jmaskMono2D_lo);
  spline(VIPERS_jmaskr_lo, VIPERS_jmaskQuad_lo, VIPERS_jmask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_jmaskQuad2D_lo);
  spline(VIPERS_jmaskr_lo, VIPERS_jmaskHex_lo,  VIPERS_jmask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_jmaskHex2D_lo);
  spline(VIPERS_jmaskr_lo, VIPERS_jmaskOct_lo,  VIPERS_jmask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_jmaskOct2D_lo);
  spline(VIPERS_jmaskr_lo, VIPERS_jmaskDec_lo,  VIPERS_jmask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_jmaskDec2D_lo);
  
  print_jwindowCorrfn();
  
  return 0;
}
