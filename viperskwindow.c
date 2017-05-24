double splint_VIPERS_kSpaceMono(double k){
  if(k<0.0001)   return 1.;
  if(k>1.0000)   return 0.;

  else{
    double Interim;

    splint(VIPERS_k, VIPERS_kMono, VIPERS_kMono2D,  VIPERS_kSpace_multipoles_lineNo, k, &Interim);

    return Interim;
  }
}

double splint_VIPERS_kSpaceQuad(double k){
  if(k<0.0001)   return 0.0;
  if(k>1.0000)   return 0.0;

  else{
    double Interim;

    splint(VIPERS_k, VIPERS_kQuad, VIPERS_kQuad2D,  VIPERS_kSpace_multipoles_lineNo, k, &Interim);

    return Interim;
  }
}

double get_kMask_norm(){
  //  Hankel transform pair counts of the window.
  FFTLog_initialise_mask(mono_config);

  // correlation_fns already assigned to the window in FFTlog_memory.
  pk_mu(mono_config);
  
  double norm = 0.0;  // norm = 4.700981*1.823239*pow(10., 6.);

  for(j=0; j<mono_config->N; j++){
    if(mono_config->krvals[j][0] > 0.0001){
      norm = mono_config->pk[j][0];

      break;
    }
  }

  return norm;
}

int printf_kMask_multipoles(){  
  double norm;  // norm = 4.700981*1.823239*pow(10., 6.);

  norm = get_kMask_norm();

  FFTLog_initialise_mask(mono_config);
  FFTLog_initialise_mask(quad_config);

  // correlation_fns already assigned to the window in FFTlog_memory.
  pk_mu(mono_config);
  pk_mu(quad_config);

  sprintf(filepath, "%s/Qmultipoles/Qlk_W%d_Nag_v7_specweight_nbar_Pfkp_4000_%.1lf_%.1lf_thread_0", maskmultipoles_path, fieldFlag, lo_zlim, hi_zlim);  
//sprintf(filepath,"%s/Qmultipoles/mask_kmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf.dat", maskmultipoles_path, fieldFlag, lo_zlim, hi_zlim);

  output = fopen(filepath, "w");

  for(j=0; j<mono_config->N; j++){
    if((mono_config->krvals[j][0] > 0.0001) && (mono_config->krvals[j][0] < 100.)){
      fprintf(output, "%le \t %le \t %le \n", mono_config->krvals[j][0], mono_config->pk[j][0]/norm, quad_config->pk[j][0]/norm);
    }
  }

  fclose(output);
  
  return 0;
}

int prepVIPERS_kSpaceMultipole(){  
  printf_kMask_multipoles();

  sprintf(filepath, "%s/Qmultipoles/Qlk_W%d_Nag_v7_specweight_nbar_Pfkp_4000_%.1lf_%.1lf_thread_0", maskmultipoles_path, fieldFlag, lo_zlim, hi_zlim);
//sprintf(filepath, "%s/Qmultipoles/mask_kmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf.dat", maskmultipoles_path, fieldFlag, lo_zlim,hi_zlim);

  inputfile = fopen(filepath, "r");

  line_count(inputfile, &VIPERS_kSpace_multipoles_lineNo);

  VIPERS_k               = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
  
  VIPERS_kMono           = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
  VIPERS_kMono2D         = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));

  VIPERS_kQuad           = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
  VIPERS_kQuad2D         = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));

  // printf("\n\nVIPERS k-space window:");
  
  for(j=0; j<VIPERS_kSpace_multipoles_lineNo; j++){
    fscanf(inputfile, "%le \t %le \t %le \n", &VIPERS_k[j], &VIPERS_kMono[j], &VIPERS_kQuad[j]);
    
    // set limits due to noise. 
    if(VIPERS_k[j] < 0.0001){
      VIPERS_kMono[j] = 1.0;
      VIPERS_kQuad[j] = 0.0;
    }

    if(VIPERS_k[j] > 1.0000){
      VIPERS_kMono[j] = 0.0;
      VIPERS_kQuad[j] = 0.0;
    }
  }
  
  fclose(inputfile);

  spline(VIPERS_k, VIPERS_kMono, VIPERS_kSpace_multipoles_lineNo, 1.0e31, 1.0e31, VIPERS_kMono2D);
  spline(VIPERS_k, VIPERS_kQuad, VIPERS_kSpace_multipoles_lineNo, 1.0e31, 1.0e31, VIPERS_kQuad2D);
  
  return 0;
}
