int load_homogeneous_rands_window(double sampling){
    // sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIG.cat", root_dir, fieldFlag, lo_zlim, hi_zlim);
    sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIG.cat", root_dir, fieldFlag, 0.6,         0.9);
  
    // Joint-field catalogue. 
    // sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W1W4_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIGSHUF.cat", root_dir, fieldFlag, lo_zlim, hi_zlim);
    // sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W1W4_xyz_%.1lf_%.1lf_Nagoya_v6_Samhain_BIGSHUF.cat", root_dir, fieldFlag, 0.6, 0.9);  

    printf("\n\nMask file: %s", filepath);
    
    inputfile = fopen(filepath, "r");

    line_count(inputfile, &rand_number);
    
    lowerSampling_randomisedCatalogue(sampling);  // assumes catalogue order is random.

    printf("\n\nsampling %.2e, %d randoms number", sampling, rand_number);

    assign_randmemory();

    for(j=0; j<rand_number; j++)  fscanf(inputfile, "%le \t %le \t %*le \t %*le \t %*le \t %*le \n", &rand_ra[j], &rand_dec[j]);

    fclose(inputfile);
    
    rand_chiReassignment();
    
    // randoms_nbar_calc();
    
    // randoms_maskGen_GriddingError();  // add a random component to the position vector of the random, on the scale of the grid size. account for FFT gridding in window calc.
        
    //** assignAcceptance_rand();  // redshift limit cuts.
    
    return 0;
}


int rand_chiReassignment(){
    // With nbar specified according to interp_nz(chi), assign chi for randoms such that they satisfy this bar.
    // Achieved with the transformation method, see pg. 287 of NR and smooth_nbar.c

    printf("\n\nRandoms chi reassignment.");

    double F; 

    for(j=0; j<rand_number; j++){      
      F              = gsl_rng_uniform(gsl_ran_r);
          
      rand_chi[j]    = inverse_cumulative_nbar(F);

      rand_weight[j] = 1./(1. + interp_nz(rand_chi[j])*fkpPk); 
      
      rand_ra[j]    *= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   *= (pi/180.0);                                 // Converted to radians.

      rand_x[j]      = rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);        
      rand_y[j]      = rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
      rand_z[j]      = rand_chi[j]*sin(rand_dec[j]);
    }

    // for(j=0; j<20; j++)  printf("\n%.6lf \t %.6lf \t %.6lf \t %.6lf", rand_x[j], rand_y[j], rand_z[j], rand_chi[j]);
    
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
