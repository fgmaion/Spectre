int assign_randmemory(){
    rand_ra       = (double *) malloc(rand_number*sizeof(*rand_ra));
    rand_dec      = (double *) malloc(rand_number*sizeof(*rand_dec));
    rand_chi      = (double *) malloc(rand_number*sizeof(*rand_chi));

    rand_x        = (double *) malloc(rand_number*sizeof(*rand_x));
    rand_y        = (double *) malloc(rand_number*sizeof(*rand_y));
    rand_z        = (double *) malloc(rand_number*sizeof(*rand_z));

    // rand_redshift = (double *) malloc(rand_number*sizeof(*rand_redshift));
    
    return 0;
}


int randoms_maskGen(){   
    rand_number = 20000000;
                  
    assign_randmemory();
    
    for(j=0; j<rand_number; j++){
         rand_ra[j] =  LowerRAlimit  + (UpperRAlimit  -  LowerRAlimit)*gsl_rng_uniform(gsl_ran_r);
        rand_dec[j] =  LowerDecLimit + (UpperDecLimit - LowerDecLimit)*gsl_rng_uniform(gsl_ran_r);
    }
    
    for(j=0; j<rand_number; j++){  
      // Index     = (int) gsl_rng_uniform_int(gsl_ran_r, redshiftsNumber);
      // rand_z[j] = redshifts[Index]; 

      Interim    = pow(loChi, 3.) + (pow(hiChi, 3.) - pow(loChi, 3.))*gsl_rng_uniform(gsl_ran_r);

      rand_chi[j] = pow(Interim, 1./3.);
    }
    
    for(j=0; j<rand_number; j++){  
      rand_ra[j]    *= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   *= (pi/180.0);                                 // Converted to radians.
    
      rand_x[j]      = rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);        
      rand_y[j]      = rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
      rand_z[j]      = rand_chi[j]*sin(rand_dec[j]);
      
      rand_ra[j]    /= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   /= (pi/180.0);                                 // Converted to radians.
    }
    
    
    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_parent_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);

    output = fopen(filepath, "w");
    
    for(j=0; j<rand_number; j++)  fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \n", rand_ra[j], rand_dec[j], rand_chi[j], rand_x[j], rand_y[j], rand_z[j]);
    
    fclose(output);
    
    printf("\n\nrandoms co-ordinates.");
    
    printf("\n%e \t %e", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\n%e \t %e", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\n%e \t %e", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    return 0;
}


int randoms_zdist(){
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

    return 0;
}
