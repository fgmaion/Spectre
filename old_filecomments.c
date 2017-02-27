//** Added 08/02/2017 **//
int spec_weights(){
  // analysis on **mock** catalogues.

  // local TSR weights for mock spec cats I made, e.g. non-Poisson sampling on a per cell basis of a grid overlaid across the survey.
  // if(loopCount<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_spec_weights.dat", root_dir, loopCount);
  // else if(loopCount<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_spec_weights.dat",  root_dir, loopCount);
  // else                    sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_spec_weights.dat",   root_dir, loopCount);

  // if(loopCount<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_mockTSR_weights.dat", root_dir, loopCount);
  // else if(loopCount<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_mockTSR_weights.dat",  root_dir, loopCount);
  // else                    sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_mockTSR_weights.dat",   root_dir, loopCount);

  // if(loopCount<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_mockSpec_nonpoisson_weights.dat", root_dir, loopCount);
  // else if(loopCount<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_mockSpec_nonpoisson_weights.dat",  root_dir, loopCount);
  // else                    sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_mockSpec_nonpoisson_weights.dat",   root_dir, loopCount);

  // if(loopCount<10)        sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_00%d_W%d_spec_weights.dat", root_dir, fieldFlag, loopCount, fieldFlag);
  // else if(loopCount<100)  sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_0%d_W%d_spec_weights.dat",  root_dir, fieldFlag, loopCount, fieldFlag);
  // else                    sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_%d_W%d_spec_weights.dat",   root_dir, fieldFlag, loopCount, fieldFlag);

  // New TSR weights by Delauney Tesselation. 7 Dec 2015.
  // if(loopCount<10)        sprintf(filepath, "%s/W%d_mocks/TSR_W%d_mock_00%d_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, fieldFlag, loopCount);
  // else if(loopCount<100)  sprintf(filepath, "%s/W%d_mocks/TSR_W%d_mock_0%d_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, fieldFlag, loopCount);
  // else                    sprintf(filepath, "%s/W%d_mocks/TSR_W%d_mock_%d_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, fieldFlag, loopCount);

  // DT TSR, for new mocks in W1_Spectro_V7_2.
  if(loopCount<10)        sprintf(filepath, "%s/mocks_v1.7/TSR/TSR_W%d_mock_00%d_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, loopCount);
  else if(loopCount<100)  sprintf(filepath, "%s/mocks_v1.7/TSR/TSR_W%d_mock_0%d_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
  else                    sprintf(filepath, "%s/mocks_v1.7/TSR/TSR_W%d_mock_%d_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);

  inputfile     = fopen(filepath, "r");

  for(j=0; j<Vipers_Num; j++){
    fscanf(inputfile, "%*d \t %*lf \t %*lf \t %lf \n", &sampling[j]);
  }

  fclose(inputfile);

  return 0;
}


int prep_fftw(){
  H_k                = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(fftw_complex));
  H2_k               = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(fftw_complex));

  p                  = fftw_plan_dft_3d(n0, n1, n2, overdensity, H_k, FFTW_FORWARD,  FFTW_ESTIMATE);
                                                                                                                                                        
  PkArray            = (double **)     malloc(n0*n1*n2*sizeof(double*));           // rows                                                                
  muIntervalPk       = (double **)     malloc(n0*n1*n2*sizeof(*muIntervalPk));                                                                              
                                                                                                                                                           
  for(j=0; j<n0*n1*n2; j++){                                                                                                                                
    PkArray[j]     = (double *)      malloc(2*sizeof(double));                     // columns                                                             
                                                                                                                                                            
    PkArray[j][0]  = 0.0;                                                                                                                                   
    PkArray[j][1]  = 0.0;                                                                                                                                
                                                                                                                                                            
    // muIntervalPk[j]    = (double *)     malloc(2*sizeof(double));                                                                                      
                                                                                                                                                             
    // muIntervalPk[j][0] = 0.0;                                                                                                                          
    // muIntervalPk[j][1] = 0.0;                                                                                                                          
  }                                                                                                                                                         
   
  // iplan              = fftw_plan_dft_3d(n0, n1, n2, H_k, overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);

  // W2_veck            = malloc(n0*n1*n2*sizeof(double));
  // W2_vecr            = malloc(n0*n1*n2*sizeof(double));

  // FFTW2_vecr_re      = malloc(n0*n1*n2*sizeof(double));
  // FFTW2_vecr_im      = malloc(n0*n1*n2*sizeof(double));

  return 0;
}

// added 10/02/2017
int  assignAcceptance(){
  int accepted = 0;

  for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = false;

  for(j=0; j<Vipers_Num; j++){
    if((lo_zlim<=gal_z[j]) && (gal_z[j]<=hi_zlim)){
      if(data_mock_flag ==0){
	// dec problem in mocks.
	if(dec[j] >= -5.97)  Acceptanceflag[j]  = true;
      }

      else{
	Acceptanceflag[j]  = true;
      }
    }
  }
  /*                                                                                                                                                                                             daccepted_gals     = 0.0;                                                                                                                                                                

    for(j=0; j<Vipers_Num; j++){                                                                                                                                                                   if(Acceptanceflag[j] == true){                                                                                                                                                                 accepted         += 1;                                                                                                                                                               
    
        daccepted_gals   += 1./sampling[j];		                                                                                                                                     
	
	// daccepted_gals   += clip_galweight[j]/sampling[j];	                                                                                                                                   }																							     }                                                                                                                                                                                         																							     
   printf("\n\nTotal number of galaxies on input: %d, accepted: %d, accepted (weighted) %.2lf", Vipers_Num, accepted, daccepted_gals); \
  */
  
  return 0;
}


int load_homogeneous_rands_window(double sampling){
  // -- randoms follow solely the geometry -- //
  // printf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_parent_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);

  // same size as Stefano's.
  // sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_stefano.cat", root_dir, fieldFlag, lo_zlim, hi_zlim);

  // redshifts reassigned below.
  sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_stefano.cat", root_dir, fieldFlag, 0.6, 0.9);

  // Checking rand_number and angular limit dependence of P(k) measured in high redshift slice of W1. -> more rigorous angular limits and more randoms.
  // sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W1_Nagoya_v6_Samhain_xyz_0.9_1.2_600000.cat", root_dir);

  // Stefano, (ra, dec) only.
  // sprintf(filepath, "%s/W1_Spectro_V7_2/stefano_dat/random/random_W1.txt", root_dir);

  printf("\n\nMask file: %s", filepath);

  inputfile   = fopen(filepath, "r");

  line_count(inputfile, &rand_number);

  lowerSampling_randomisedCatalogue(sampling);  // assumes catalogue order is random.

  assign_randmemory();

  for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le \n", &rand_ra[j], &rand_dec[j], &rand_chi[j], &rand_x[j], &rand_y[j], &rand_z[j]);

  fclose(inputfile);

  rand_chiReassignment();  // Given <n(z)>, assign chi for randoms such that they satisfy this.

  // randoms_nbar_calc();

  // randoms_maskGen_GriddingError();  // add perturbation to position of random, length scale: grid size, to account for FFT gridding.

  assignAcceptance_rand();  // redshift limit cuts

  for(j=0; j<rand_number; j++)  rand_weight[j] = 1./(1. + interp_nz(rand_chi[j])*fkpPk);  // rand fkp weights.

  StefanoReflection(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);

  StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);

  printf("\n\nStefano basis, randoms co-ordinates.");

  printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));

  Jenkins_foldRand();

  return 0;
}
