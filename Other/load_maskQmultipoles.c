int set_rand_rng(void){
  for(j=0; j<rand_number; j++)  rand_rng[j] = gsl_rng_uniform(gsl_ran_r);

  return 0;
}

int rotate_ra(void){
  // Stefano basis has ra aligned with y and redshift along z; R1 rotation can be done simply by changing ra coordinate,
  // this is NOT true of dec.
  for(j=0; j<rand_number; j++)  rand_ra[j] -= CentreRA;

  return 0;
}

int rand_chiReassignment(){
  // With nbar specified according to interp_nz(chi), assign chi for randoms such that they satisfy this bar.
  // Achieved with the transformation method, see pg. 287 of NR and smooth_nbar.c

  double  cos_dec;
  double   x2, z2;
  double    c_dec;

  c_dec   = CentreDec*(pi/180.);

  rotate_ra();
  
  printf("\n\nRandoms chi reassignment; centre RA: %.3lf \t centre DEC: %.3lf", CentreRA, CentreDec);
  
  // #pragma omp parallel for private(j, cos_dec, x2, z2) if(thread == 1)
  for(j=0; j<rand_number; j++){
    rand_chi[j]    = inverse_cumulative_nbar(rand_rng[j]);

    rand_ra[j]    *= (pi/180.0); // Converted to radians.  No need to convert back.
    rand_dec[j]   *= (pi/180.0);

    cos_dec        = cos(rand_dec[j]);
      
    rand_x[j]      =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
    rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
    rand_z[j]      = -rand_chi[j]*sin(rand_dec[j]);

    rand_weight[j] = 1./(1. + (*pt2nz)(rand_chi[j])*fkpPk);   // rand fkp weights.
    
    x2             = -sin(c_dec)*rand_x[j]  - cos(c_dec)*rand_z[j];
    z2             =  cos(c_dec)*rand_x[j]  - sin(c_dec)*rand_z[j];

    rand_x[j]      = x2 + stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
    rand_y[j]     +=      stefano_trans_y;
    rand_z[j]      = z2 + stefano_trans_z;   
  }

  printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));     
  
  return 0;
}

int make_fastread_randomCats(int rand_number){
  //  Output ra and dec only, in binary. assumes loaded already.
  output = fopen(filepath, "wb");

  fwrite(rand_ra,   sizeof(double),  rand_number,   output);
  fwrite(rand_dec,  sizeof(double),  rand_number,   output);

  fclose(output);

  // load_fastread_randomCats(rand_number);

  return 0;
}

int load_fastread_randomCats(int rand_number){
  inputfile = fopen(filepath, "rb");

  fread(rand_ra,  sizeof(double), rand_number, inputfile);
  fread(rand_dec, sizeof(double), rand_number, inputfile);

  fclose(inputfile);

  return 0;
}

int load_ascii_randomCats(int rand_number){
  inputfile = fopen(filepath, "r");

  for(j=0; j<rand_number; j++)  fscanf(inputfile, "%le \t %le \t %*le \t %*le \t %*le \t %*le \n", &rand_ra[j], &rand_dec[j]);

  fclose(inputfile);

  return 0;
}

int lowerSampling_randomisedCatalogue(double sampling){
    rand_number = (int) ceil(rand_number*sampling);

    return 0;
}

int load_homogeneous_rands_window(double sampling, int count_res){  
  if(count_res < 3){
    sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIG.dat", root_dir, fieldFlag, 0.6, 0.9); // chis to be reset. 
  
    if(strcmp(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/randoms/randoms_W1_xyz_0.6_0.9_Nagoya_v6_Samhain_BIG.dat") == 0)       rand_number = 27670652;      // Data: 50829
    else if(strcmp(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/randoms/randoms_W4_xyz_0.6_0.9_Nagoya_v6_Samhain_BIG.dat") == 0)  rand_number = 25593777;      //     : 26377   
    else{
      printf("\n\nFile not found: %s", filepath);

      exit(EXIT_FAILURE);
    }
  }
  
  else{
    sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W1W4_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIGSHUF.dat", root_dir, 0.6, 0.9);       // joint field

    if(strcmp(filepath, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/randoms/randoms_W1W4_xyz_0.6_0.9_Nagoya_v6_Samhain_BIGSHUF.dat") == 0) rand_number = 53264429;     // Data: 50829
    else{
      printf("\n\nFile not found: %s", filepath);

      exit(EXIT_FAILURE);
    }
  }

  lowerSampling_randomisedCatalogue(sampling);  // assumes catalogue order is random.
 
  printf("\n\nSampling %.2e, %d randoms number", sampling, rand_number);
  
  walltime("Wall time at start of mask loading.");

  assign_rand_radecmemory();

  assign_randmemory();
  
  walltime("Wall time at end of mask loading.");

  // load_ascii_randomCats(rand_number);  
  load_fastread_randomCats(rand_number);
  
  // make_fastread_randomCats(rand_number); 
  
  set_rand_rng();
  
  rand_chiReassignment();
     
  return 0;
}

int load_maskfits(double sampling, int count_res){
  // e.g. nagoya_v6_samhain_W1W4_all_deccut.fits
  if(count_res<3) sprintf(filepath, "/home/mjw/venice-4.0.2/nagoya_v6_samhain_W%d_all_deccut.fits", fieldFlag);

  else{
    sprintf(filepath, "/home/mjw/venice-4.0.2/nagoya_v6_samhain_W1W4_all_deccut.fits");          
  }

  walltime("Wall time at start of mask loading.");
  
  read_maskfits(filepath, sampling);  // rand_number set by nrows.
  
  printf("\n\nMask: %s with sampling %.4lf has %d randoms", filepath, sampling, rand_number);
  
  assign_randmemory();

  walltime("Wall time at end of mask loading.");
  
  set_rand_rng();

  rand_chiReassignment();
  
  return 0;
}


int load_maskedRSD_rands(double sampling){
  sprintf(filepath, "/home/mjw/maskedRSD/dat/randoms_W1_Nagoya_xyz_0.7_0.8_gridded.cat");

  inputfile = fopen(filepath, "r");

  line_count(inputfile, &rand_number);  

  lowerSampling_randomisedCatalogue(sampling);
  
  assign_randmemory();
  
  for(j=0; j<rand_number; j++){
    fscanf(inputfile, "%le \t %le \t %le", &rand_x[j], &rand_y[j], &rand_z[j]);
    
    rand_weight[j]  = 1.0;    
  }

  return 0;
}
