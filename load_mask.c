int load_homogeneous_rands_window(int load, double sampling, int flag_data_mock){
    // -- randoms follow solely the geometry -- //
    // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_parent_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);
    // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);
    // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_gridded_3.12_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim); 
    
    
    // -- randoms follow the number density, specific nbar implies distinction of mask between mocks and data -- //
    if(flag_data_mock == 0){
        // analysis on **mock** catalogues. 
    
        // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_incnbar_xyz_0.6_0.9.cat");
        // sprintf(filepath, "%s/W1_Spectro_V5_0/randoms_W1_Nagoya_v5_Samhain_realdata_incnbar_xyz_%.1lf_%.1lf.cat", root_dir, lo_zlim, hi_zlim);
        // sprintf(filepath, "%s/W1_Spectro_V5_0/randoms_W1_Nagoya_v5_Samhain_incmocknbar_xyz_0.6_0.9.cat", root_dir, lo_zlim, hi_zlim);
    
        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/randoms/randoms_W1_Nagoya_v6_samhain_incmock_specweight_nbar_xyz_%.1lf_%.1lf.cat", root_dir, lo_zlim, hi_zlim);
    
        // Angular coverage of mocks and data should be identical.
        sprintf(filepath, "%s/W1_Spectro_V7_0/randoms_W%d_Nagoya_V7_Samhain_incmock_specweight_nbar_xyz_%.1lf_%.1lf.cat", root_dir, fieldFlag, lo_zlim, hi_zlim);
    }
    
    if(flag_data_mock == 1){
        // analysis on **VIPERS DATA** catalogues. 
        
        // Nagoya v7., Samhain, local TSR weighted nbar.         
        sprintf(filepath, "%s/W1_Spectro_V7_0/randoms_W%d_Nagoya_V7_Samhain_incmock_specweight_nbar_xyz_%.1lf_%.1lf.cat", root_dir, fieldFlag, lo_zlim, hi_zlim);
    
        printf("\n\nMask file: %s", filepath);
    }
    
    
    if(flag_data_mock == 2){
        // analysis on **VIPERS DATA** catalogues, combined W1 and W4.
        
        // Nagoya v7., Samhain, local TSR weighted nbar.         
        sprintf(filepath, "%s/W1_Spectro_V7_0/randoms_W1W4_Nagoya_V7_Samhain_incmock_specweight_nbar_xyz_%.1lf_%.1lf.cat", root_dir, lo_zlim, hi_zlim);
    
        printf("\n\nMask file: %s", filepath);
    }
    
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
    
    // With nbar specified according to interp_nz(chi), assign chi for randoms such that they satisfy this n bar.
    rand_chiReassignment();
    
    // randoms_nbar_calc();
    
    // add a random component to the position vector of the random, on the scale of the grid size. account for FFT gridding in window calc. 
    // randoms_maskGen_GriddingError();
    
    // rand fkp weights. 
    // for(j=0; j<rand_number; j++)  rand_weight[j] = 1.;
    for(j=0; j<rand_number; j++)  rand_weight[j] = 1./(1. + interp_nz(rand_chi[j])*fkpPk);
    
    StefanoReflection(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
    
    StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
    
    printf("\n\nStefano basis, randoms co-ordinates.");
    
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    // mean_CellChi();
    
    // redshift limit cuts. 
    assignAcceptance_rand();
    
    return 0;
}


int rand_chiReassignment(){
    // With nbar specified according to interp_nz(chi), assign chi for randoms such that they satisfy this bar.
    // Achieved with the transformation method, see pg. 287 of NR and smooth_nbar.c

    double F; 

    for(j=0; j<rand_number; j++){ 
        rand_chi[j]   = 0.0;
           
        while((rand_chi[j]<loChi) || (hiChi<rand_chi[j])){
          F           = gsl_rng_uniform(gsl_ran_r);
          
          rand_chi[j] = inverse_cumulative_nbar(F);
        }
            
        rand_ra[j]   *= (pi/180.0);                                 // Converted to radians.
        rand_dec[j]  *= (pi/180.0);                                 // Converted to radians.
    
        rand_x[j]     = rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);        
        rand_y[j]     = rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
        rand_z[j]     = rand_chi[j]*sin(rand_dec[j]);
      
        rand_ra[j]   /= (pi/180.0);                                 // Converted to radians.
        rand_dec[j]  /= (pi/180.0);                                 // Converted to radians.
    }

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


int lowerSampling_randomisedCatalogue(double sampling){
    rand_number = (int) ceil(rand_number*sampling);

    return 0;
}
