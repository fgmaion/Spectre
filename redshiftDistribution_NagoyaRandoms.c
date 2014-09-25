int randoms_redshiftDistribution(){  
    int redshiftsNumber = 0;
    
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
  
    // for(j=0; j<10; j++)  printf("\n%e", redshifts[j]);

    // load randoms to be assigned redshifts. 

    sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/mock_W1_rand_Nagoya_join.cat");

    inputfile = fopen(filepath, "r");
    
    int howmanyRandoms;

    howmanyRandoms = 0;

    ch         = 0;

    do{
      ch = fgetc(inputfile);
        if(ch == '\n')  howmanyRandoms += 1;
    } while (ch != EOF);

    rewind(inputfile);

    double* rand_ra;
    double* rand_dec;
    double* rand_z;
    
    printf("\n\n%d rands", howmanyRandoms);
    
    rand_ra  =  malloc(howmanyRandoms*sizeof(double));
    rand_dec =  malloc(howmanyRandoms*sizeof(double));
    rand_z   =  malloc(howmanyRandoms*sizeof(double));
    
    for(j=0; j<howmanyRandoms; j++)  fscanf(inputfile, "%le \t %le \n", &rand_ra[j], &rand_dec[j]);

    fclose(inputfile);
   
    // for(j=0; j<10; j++)  printf("\n%e \t %e ", rand_ra[j], rand_dec[j]);
    
    // printf("\n\n%d redshifts", redshiftsNumber);

    // for(j=0; j<10; j++) printf("\n%d", (int) gsl_rng_uniform_int(gsl_ran_r, 66548));

    for(j=0; j<howmanyRandoms; j++){  
      Index     = (int) gsl_rng_uniform_int(gsl_ran_r, redshiftsNumber);

      rand_z[j] = redshifts[Index];
    
      // printf("\n%d \t %e", Index, rand_z[j]);
    }
   
    // if(loopCount<10)  sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms_W1_Nagoya_0.65_0.75_v1.2/randoms_W1_Nagoya_ra_dec_z_00%d.cat", loopCount);
    // else sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms_W1_Nagoya_0.65_0.75_v1.2/randoms_W1_Nagoya_ra_dec_z_0%d.cat", loopCount);

    sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms_W1_Nagoya_0.65_0.75_v1.2/randoms20_W1_Nagoya_ra_dec_z_.cat");

    output = fopen(filepath, "w");

    for(j=0; j<howmanyRandoms; j++){  
      fprintf(output, "%e \t %e \t %e \n", rand_ra[j], rand_dec[j], rand_z[j]);
    }

    fclose(output);
    
    return 0;
}
