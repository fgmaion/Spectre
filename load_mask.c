int rotate_ra(void){
  // Stefano basis has ra aligned with y and redshift along z; R1 rotation can be done simply by changing ra coordinate,
  // this is NOT true of dec. 
  for(j=0; j<rand_number; j++)  rand_ra[j] -= CentreRA;

  return 0;
}
  
int load_rands_radec(double sampling){
  rand_number = accepted_rand = (int) ceil(1382582*sampling);  // Hard coded catalogue max row number.

  assign_rand_radecmemory();
  
  assign_randmemory();
    
  // load_ascii_randomCats(sampling);    
  load_fastread_randomCats(rand_number);

  // make_fastread_randomCats();

  rotate_ra();
  
  // To do: save randoms in radians.
  for(j=0; j<rand_number; j++){
    rand_ra[j]  *= (pi/180.0); // Converted to radians.  No need to convert back.
    rand_dec[j] *= (pi/180.0);
  }

  set_rand_rng();
    
  walltime("Wall time after randoms load");
  
  return 0;
}

int set_rand_rng(void){
  for(j=0; j<rand_number; j++)  rand_rng[j] = gsl_rng_uniform(gsl_ran_r);

  return 0;
}

int set_cnst_nbar(void){
  // Reassign randoms to have constant \bar n(z); this will (much) better sample high z tail.
  pt2nz = &unity; // constant <n(z)>

  prep_inverseCumulative_nbar();
  
  return 0;
}

int rand_newchi_newbasis(void){
  // With nbar specified according to interp_nz(chi), assign chi for randoms such that they satisfy this bar.
  // Achieved with the transformation method, see pg. 287 of NR and smooth_nbar.c
  
  // printf("\n\nAngular limits of randoms: %.4lf < ra < %.4lf, %.4lf < dec < %.4lf", arrayMin(rand_ra, rand_number),  arrayMax(rand_ra,  rand_number),
  //                                                                                  arrayMin(rand_dec, rand_number), arrayMax(rand_dec, rand_number));  

  // printf("\n\nNew basis for randoms.");
  
  double     cos_dec;
  double      x2, z2;  
  double       c_dec;
  
  c_dec   = CentreDec*(pi/180.);
  
  #pragma omp parallel for private(j, cos_dec, x2, z2) if(thread == 1)
  for(j=0; j<rand_number; j++){
    rand_chi[j]    = inverse_cumulative_nbar(rand_rng[j]);
      
    cos_dec        = cos(rand_dec[j]);
    
    rand_x[j]      =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
    rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
    rand_z[j]      = -rand_chi[j]*sin(rand_dec[j]);           // Stefano reflection included. 
      
    rand_weight[j] = 1./(1. + (*pt2nz)(rand_chi[j])*fkpPk);   // rand fkp weights.
      
    // basis formed by: normal spherical co-ordinates subject to inversion through xy plane, then R1 and finally R2.
    // R1: rotation about z such that in the new basis, (x',y',z'), x' hat lies in x-y plane at an angle centreRA to x.
    // x1  =     cos(c_ra)*rand_x[j] + sin(c_ra)*rand_y[j]; 
    // y1  =    -sin(c_ra)*rand_x[j] + cos(c_ra)*rand_y[j]; //** RA rotation can be done explicitly. 
    // z1  =               rand_z[j];

    // R2: rotation about y such that in the new basis, (x'', y'', z''), z'' hat lies in (x', z') plane at an angle -CentreDec to x' hat.
    x2  = -sin(c_dec)*rand_x[j]  - cos(c_dec)*rand_z[j];
    z2  =  cos(c_dec)*rand_x[j]  - sin(c_dec)*rand_z[j];

    rand_x[j]  = x2 + stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
    rand_y[j] +=      stefano_trans_y;
    rand_z[j]  = z2 + stefano_trans_z;
  }
 
  // printf("\n\nRandoms: Stefano basis.");                                                                                                      
  // printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  // printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  // printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));                                                   

  // walltime("Wall time after randoms chi reassignment");
  
  return 0;
}

int assign_randbox(){
  for(j=0; j<rand_number; j++){
    xlabel     = (int)  floor((rand_x[j] - min_x)/dx);
    ylabel     = (int)  floor((rand_y[j] - min_y)/dy);
    zlabel     = (int)  floor((rand_z[j] - min_z)/dz);

    rand_box[j] = (int)  xlabel + n2*ylabel + n2*n1*zlabel;
  }

  return 0;
}

int lowerSampling_randomisedCatalogue(double sampling){
  rand_number = (int) ceil(rand_number*sampling);

  return 0;
}

int make_fastread_randomCats(void){
  //  Output ra and dec only, in binary. assumes loaded already. 
  sprintf(filepath, "%s/W1_Spectro_V7_4/randoms/randoms_W%d_Nagoya_v6_Samhain_stefano.cat", root_dir, fieldFlag);
  
  output = fopen(filepath, "wb");

  fwrite(rand_ra,   sizeof(double),  rand_number,   output);
  fwrite(rand_dec,  sizeof(double),  rand_number,   output);
  
  fclose(output);

  // load_fastread_randomCats(rand_number);
  
  return 0;
}

int load_fastread_randomCats(int rand_number){  
  sprintf(filepath, "%s/W1_Spectro_V7_4/randoms/randoms_W%d_Nagoya_v6_Samhain_stefano.cat", root_dir, fieldFlag);

  inputfile = fopen(filepath, "rb");

  fread(rand_ra,  sizeof(double), rand_number, inputfile);
  fread(rand_dec, sizeof(double), rand_number, inputfile);

  fclose(inputfile);
  
  return 0;
}

double cmpfunc(const void* a, const void* b){
  return ( *(double*) a - *(double*) b);
}

int cut_rand_bydec(){
  // Sorting only needs to be done once, ever.
  gsl_sort2(rand_dec, 1, rand_ra, 1, rand_number);

  double store_ra, store_dec;
  
  for(j=0; j<rand_number/2; j++){ // reverse order. 
    store_ra                      = rand_ra[j]; 
    rand_ra[j]                    = rand_ra[rand_number - (j+1)];
    rand_ra[rand_number - (j+1)]  = store_ra;

    store_dec                     = rand_dec[j];
    rand_dec[j]                   = rand_dec[rand_number - (j+1)];
    rand_dec[rand_number - (j+1)] = store_dec;
  }

  assignAcceptance_rand();
  
  // Sort by dec to ensure easy application of dec flag. 
  make_fastread_randomCats();  // Create binary for faster reading.                                                                                                                           
  return 0;
}

int load_ascii_randomCats(double sampling){
  sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_stefano.cat", root_dir, fieldFlag, 0.6, 0.9);  // redshifts reassigned below.                 
                                                                                                                                                                                            
  printf("\n\nMask file: %s", filepath);                                                                                                                                                    
                                                                                                                                                                                             
  inputfile   = fopen(filepath, "r");                                                                                                                                                       
                                                                                                                                                                                            
  line_count(inputfile, &rand_number);          // Provide line count to save time.                                                                                                                                                                                                                                                                                                
  lowerSampling_randomisedCatalogue(sampling);  // assumes catalogue order is random.                                                                                                       
                                                                                                                                                                                             
  for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %*e \t %*e \t %*e \t %*e \n", &rand_ra[j], &rand_dec[j]);                                                                                                                                                                                                                                                       
  fclose(inputfile);                                                                                                                                                                     

  // cut_rand_bydec();
  
  return 0;
}

int del_lockfile(){
  // Delete lockfile if necessary, e.g. for AP MCMC.                                                                                                        
  char*       lockfile_path;                                                                                                                               
  char delete_lockfile[200];                                                                                                                              

  lockfile_path = malloc(200*sizeof(*lockfile_path));                                                                                                        
                                                                                                                                                            
  lockfile_path = getenv("LOCKFILEDIR");                                                                                                                    
  
  sprintf(delete_lockfile, "rm -rf %s", lockfile_path);                                                                                                      
  
  system(delete_lockfile);                                                                                                                                  

  system("echo 'lockfile destroyed.'");                                                                                                                      

  return 0;
}
