int rand_Gen(){  
    // Given distribution of redshifts and angular distributin create randoms catalogue.
    /*
    int     redshiftsNumber = 0;
    double* redshifts;

    sprintf(filepath, "%s/redshiftDistribution_NagoyaRandoms.dat", root_dir);

    inputfile = fopen(filepath, "r");

    ch         = 0;
       
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  redshiftsNumber += 1;
    } while (ch != EOF);

    printf("\n\n%d redshifts. up to date", redshiftsNumber);
    
    rewind(inputfile);
    
    redshifts =  malloc(redshiftsNumber*sizeof(double));

    for(j=0; j<redshiftsNumber; j++)  fscanf(inputfile, "%le \n", &redshifts[j]);
    
    fclose(inputfile);
    */


    // load randoms to be assigned redshifts. 
    // sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/mock_W1_rand_Nagoya_join.cat");
    // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Venice/venice_3.8.3/mock_W1_rand_Nagoya_newnew.cat");
    /*
    inputfile = fopen(filepath, "r");
    
    rand_number = 0;
    ch          = 0;

    do{
      ch = fgetc(inputfile);
        if(ch == '\n')  rand_number += 1;
    } while (ch != EOF);

    rewind(inputfile);

    double* rand_ra;
    double* rand_dec;
    double* rand_chi;
    
    printf("\n\n%d rands \n", rand_number);
    */
    
    rand_number = 10000000;
    
    rand_ra     =  malloc(rand_number*sizeof(double));
    rand_dec    =  malloc(rand_number*sizeof(double));
    rand_chi    =  malloc(rand_number*sizeof(double));
    
    rand_x      =  malloc(rand_number*sizeof(double));
    rand_y      =  malloc(rand_number*sizeof(double));
    rand_z      =  malloc(rand_number*sizeof(double));
    
    // for(j=0; j<rand_number; j++)  fscanf(inputfile, "%le \t %le \n", &rand_ra[j], &rand_dec[j]);

    // fclose(inputfile);
    
    // 500' mocks, 30.1 < right ascension < 38.8;
    //            -5.95 <     dec         < -4.15; 
    
    for(j=0; j<rand_number; j++){
         rand_ra[j] =  LowerRAlimit  + (UpperRAlimit  -  LowerRAlimit)*gsl_rng_uniform(gsl_ran_r);
        rand_dec[j] =  LowerDecLimit + (UpperDecLimit - LowerDecLimit)*gsl_rng_uniform(gsl_ran_r);
    }
    
    // redshift set 0.7 to 0.8, change the cosmology! 
    for(j=0; j<rand_number; j++){  
      // Index     = (int) gsl_rng_uniform_int(gsl_ran_r, redshiftsNumber);
      // rand_z[j] = redshifts[Index]; 

      Interim    = pow(interp_comovingDistance(redshiftLowLimit), 3.) + (pow(interp_comovingDistance(redshiftHiLimit), 3.) - pow(interp_comovingDistance(redshiftLowLimit), 3.))*gsl_rng_uniform(gsl_ran_r);

      rand_chi[j] = pow(Interim, 1./3.);
    }
    
    // if(loopCount<10)  sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms_W1_Nagoya_0.65_0.75_v1.2/randoms_W1_Nagoya_ra_dec_z_00%d.cat", loopCount);
    // else sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms_W1_Nagoya_0.65_0.75_v1.2/randoms_W1_Nagoya_ra_dec_z_0%d.cat", loopCount);

    printf("\n\nWriting randoms co-ordinates.");
    
    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/likelihood/randoms_W1_500s_xyz_0.7_0.8.cat");

    output = fopen(filepath, "w");
    
    for(j=0; j<rand_number; j++){  
      rand_ra[j]    *= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   *= (pi/180.0);                                 // Converted to radians.
    
      rand_x[j]      =     rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);        
      rand_y[j]      =     rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
      rand_z[j]      =    -rand_chi[j]*sin(rand_dec[j]);
      
      rand_ra[j]    /= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   /= (pi/180.0);                                 // Converted to radians.
    }
    
    // centre
    StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
    
    for(j=0; j<rand_number; j++)  fprintf(output, "%e \t %e \t %e \t %e \t %e \n", rand_ra[j], rand_dec[j], rand_x[j], rand_y[j], rand_z[j]);
    
    // fprintf(output, "%e \t %e \t %e \n", rand_x[j]+CellSize*(gsl_rng_uniform(gsl_ran_r)-0.5), rand_y[j]+CellSize*(gsl_rng_uniform(gsl_ran_r)-0.5), rand_z[j]+CellSize*(gsl_rng_uniform(gsl_ran_r)-0.5));
    
    fclose(output);
    
    printf("\n%e \t %e", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\n%e \t %e", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\n%e \t %e", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    return 0;
}
