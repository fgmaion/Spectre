int load_rands_radec(double sampling){
  rand_number = accepted_rand = (int) ceil(1382582*sampling);  // Hard coded catalogue row number.
    
  assign_randmemory();
    
  // load_ascii_randomCats(sampling);    
  load_fastread_randomCats(rand_number);

  // make_fastread_randomCats();
  
  // To do: save randoms in radians.
  for(j=0; j< rand_number; j++){
    rand_ra[j]    *= (pi/180.0);                  // Converted to radians.  No need to convert back.
    rand_dec[j]   *= (pi/180.0);
  }

  set_rand_rng();
    
  walltime("Wall time after randoms load");
  
  return 0;
}


int set_rand_rng(){
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
  
  double F, cos_dec;

  double x1,  y1, z1;
  double x2,  y2, z2;
  double c_ra, c_dec;

  int new;
  
  c_ra    =  CentreRA*(pi/180.);
  c_dec   = CentreDec*(pi/180.);
  
  #pragma omp parallel for private(j, cos_dec, x1, y1, z1, x2, y2, z2) if(thread == 1)
  for(j=0; j<rand_number; j++){
    rand_chi[j]    = inverse_cumulative_nbar(rand_rng[j]);
      
    cos_dec        = cos(rand_dec[j]);
    
    rand_x[j]      =  rand_chi[j]*cos(rand_ra[j])*cos_dec;
    rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos_dec;
    rand_z[j]      = -rand_chi[j]*sin(rand_dec[j]);            // Stefano reflection included. 
      
    rand_weight[j] = 1./(1. + interp_nz(rand_chi[j])*fkpPk);  // rand fkp weights.
      
    // basis formed by: normal spherical co-ordinates subject to inversion through xy plane, then R1 and finally R2.
    // R1: rotation about z such that in the new basis, (x',y',z'), x' hat lies in x-y plane at an angle centreRA to x.
    x1  =     cos(c_ra)*rand_x[j] + sin(c_ra)*rand_y[j];
    y1  =    -sin(c_ra)*rand_x[j] + cos(c_ra)*rand_y[j];
    z1  =               rand_z[j];

    // R2: rotation about y such that in the new basis, (x'', y'', z''), z'' hat lies in (x', z') plane at an angle -CentreDec to x' hat.
    x2  = -sin(c_dec)*x1  - cos(c_dec)*z1;
    y2  =  y1;
    z2  =  cos(c_dec)*x1  - sin(c_dec)*z1;

    rand_x[j] = x2 + stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
    rand_y[j] = y2 + stefano_trans_y;
    rand_z[j] = z2 + stefano_trans_z;
  }
 
  printf("\n\nRandoms: Stefano basis.");                                                                                                      

  printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
  printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
  printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));                                                     
  
  // walltime("Wall time after randoms chi reassignment");
  
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
                                                                                                                                                                                             
  for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %*le \t %*le \t %*le \t %*le \n", &rand_ra[j], &rand_dec[j]);                                                                                                                                                                                                                                                       
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


int load_maskedRSDpaper_mask(double sampling){
    sprintf(filepath, "%s/Data/maskedRSD_draftwork/randoms_W1_Nagoya_xyz_0.7_0.8_gridded.cat", root_dir);

    inputfile   = fopen(filepath, "r");

    ch          = 0;
    rand_number = 0;

    do{
        ch = fgetc(inputfile);
        
        if(ch == '\n')
            rand_number += 1;
    } while(ch != EOF);

    printf("\n\n%d randoms number", rand_number);

    rewind(inputfile);

    lowerSampling_randomisedCatalogue(sampling);

    assign_randmemory();

    for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %le \n", &rand_x[j], &rand_y[j], &rand_z[j]);
    
    fclose(inputfile);
    
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    /*
    for(j=0; j<n0*n1*n2; j++)  surveyMask[j] = 0.0;
    
    for(j=0; j<rand_number; j++){
        boxlabel = boxCoordinates(rand_x, rand_y, rand_z, j);
    
        surveyMask[boxlabel] = 1.;
    }
    */
    for(j=0; j<rand_number; j++)  rand_weight[j] = 1.;


    sprintf(filepath, "");

    return 0;
}


int zeff_calc(){
    double Norm = 0.0;
    double zz;
    
    Interim     = 0.0;
    
    for(j=0; j<rand_number; j++){
        zz      = interp_inverseComovingDistance(rand_chi[j]);
        
        Interim += rand_weight[j]*zz;
    
        Norm    += rand_weight[j];
    }

    z_eff = Interim/Norm;

    printf("\n\nEffective redshift: %.2lf", z_eff);

    return 0;
}


int mean_CellChi(){
    // assign memory for grid representation of mask, Cell_SurveyLimitsMask. 
    double* meanchi;
    double*surveyMask;
    
    double  x, y, z, chi;
    
    meanchi = malloc(n0*n1*n2*sizeof(double));    
    
    for(j=0; j<n0*n1*n2; j++)  meanchi[j] = 0.0;
    
    for(j=0; j<rand_number; j++){
        boxlabel             = boxCoordinates(rand_x, rand_y, rand_z, j);
        
        meanchi[boxlabel]   += rand_chi[j];
    }
    
    for(j=0; j<n0*n1*n2; j++)  meanchi[j] /= surveyMask[j];


    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               Index = i + n2*j + n2*n1*k;
            
               if(surveyMask[Index] > 0.0){
                  x = (i + 0.5)*xCellSize;
                  y = (j + 0.5)*yCellSize;
                  z = (k + 0.5)*zCellSize;
               
                  chi   = invert_StefanoBasis(CentreRA, CentreDec, &x, &y, &z);
            
                  printf("\n%e \t %e", meanchi[Index], chi);         
              }
           }
        }
    }
    
    free(meanchi);
    
    return 0;
}


int load_randStefanoCoordinates(int load, double sampling){
    sprintf(filepath, "%s/W1_Spectro_V5_0/randoms_W1_Nagoya_v5_Samhain_realdata_incnbar_StefanoCoordinates_xyz_%.1lf_%.1lf.cat", root_dir, lo_zlim, hi_zlim);

    inputfile   = fopen(filepath, "r");

    ch          = 0;
    rand_number = 0;

    do{
        ch = fgetc(inputfile);
        
        if(ch == '\n')
            rand_number += 1;
    } while(ch != EOF);
    
    // assumes catalogue is in a random order.
    lowerSampling_randomisedCatalogue(sampling);

    printf("\n\nsampling %.2e, %d randoms number", sampling, rand_number);

    if(load == 1){
        rewind(inputfile);

        assign_randmemory();

        for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le \n", &rand_ra[j], &rand_dec[j], &rand_chi[j], &rand_x[j], &rand_y[j], &rand_z[j]);
    }

    fclose(inputfile);

    printf("\n\nStefano basis, randoms co-ordinates.");
    
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));

    assignAcceptance_rand();

    return 0;
}

