int prep_CatalogueInput_500s(){
  // Maximum number of galaxies present in any mock of the collection (i.e. those for covariance estimate).
  if(strcmp(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W1") == 0)  max_gals = 61765;  // no z cuts.
  
  else max_gals = max_gal();

  // Initialise to max. memory required to process all mocks.
  ra             =  (double *)  malloc(max_gals*sizeof(*ra));
  dec            =  (double *)  malloc(max_gals*sizeof(*dec));
  zobs           =  (double *)  malloc(max_gals*sizeof(*zobs));
  Acceptanceflag =  (bool  *)   malloc(max_gals*sizeof(*Acceptanceflag));
  rDist          =  (double *)  malloc(max_gals*sizeof(*rDist));
  xCoor          =  (double *)  malloc(max_gals*sizeof(*xCoor));
  yCoor          =  (double *)  malloc(max_gals*sizeof(*yCoor));
  zCoor          =  (double *)  malloc(max_gals*sizeof(*zCoor));
  sampling       =  (double *)  malloc(max_gals*sizeof(*sampling));
  fkp_galweight  =  (double *)  malloc(max_gals*sizeof(*fkp_galweight));
  clip_galweight =  (double *)  malloc(max_gals*sizeof(*clip_galweight));

  gal_z = &zobs[0];  // Choice of redshift from zcos, zpec, zphot, zobs.       

  return 0;
}


int spec_weights(){
  // DT TSR, for new mocks in W1_Spectro_V7_2.
  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/TSR/TSR_W%d_mock_%03d_Nagoya_v6_Samhain.dat", root_dir, fieldFlag, loopCount);

  inputfile = fopen(filepath, "r");

  for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%*d \t %*lf \t %*lf \t %lf \n", &sampling[j]);

  fclose(inputfile);

  return 0;
}


int CatalogueInput_500s(){  
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile   = fopen(filepath, "r");

    line_count(inputfile, &Vipers_Num);
    
    for(j=0; j<Vipers_Num; j++){  
      clip_galweight[j] = 1.0;
        
      fscanf(inputfile, "%*d \t %le \t %le \t %le \t %*le \t %*le \t %*le \t %*d \t %*d \t %*d \n", &ra[j], &dec[j], &zobs[j]);
      
       ra[j]               *= (pi/180.0);                                 // Converted to radians.
      dec[j]               *= (pi/180.0);                                 // Converted to radians. 

      rDist[j]              = interp_comovingDistance(gal_z[j]);          // Comoving distances in h^-1 Mpc

      xCoor[j]              =  rDist[j]*cos(dec[j])*cos(ra[j]);
      yCoor[j]              =  rDist[j]*cos(dec[j])*sin(ra[j]);
      zCoor[j]              = -rDist[j]*sin(dec[j]);                      // reflection of spherical coordinates. 

      ra[j]                /= (pi/180.0);                                 // Converted to degrees.
      dec[j]               /= (pi/180.0);                                 // Converted to degrees.  
    }
    
    fclose(inputfile);
    
    // Load ESR weights.  
    spec_weights();    // load sampling according to local TSR.

    return 0;
}


int DataInput(){
    // Also loads sampling: ESR = TSR x SSR.  
    double TSR, SSR;

    sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/vipers_lss_v7_W%d.txt", root_dir, fieldFlag);

    printf("\n\nOpening catalogue: %s\n", filepath);
    
    inputfile     = fopen(filepath, "r"); 

    line_count(inputfile, &Vipers_Num);
    
    // (id, ra, dec, zspec, TSR, SSR).                                                                                                            
    printf("\n\nNumber of galaxies in catalogue:  %d", Vipers_Num);
    
    id             =  (int    *)  realloc(id,             Vipers_Num*sizeof(*id));
    ra             =  (double *)  realloc(ra,             Vipers_Num*sizeof(*ra));
    dec            =  (double *)  realloc(dec,            Vipers_Num*sizeof(*dec));
    zobs           =  (double *)  realloc(zobs,           Vipers_Num*sizeof(*zobs)); 
    
    // derived parameters. 
    Acceptanceflag =  (bool  *)   realloc(Acceptanceflag, Vipers_Num*sizeof(*Acceptanceflag));
    rDist          =  (double *)  realloc(rDist,          Vipers_Num*sizeof(*rDist));
    xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));

    sampling       =  (double *)  malloc(Vipers_Num*sizeof(*sampling));
    fkp_galweight  =  (double *)  malloc(Vipers_Num*sizeof(*fkp_galweight));
    clip_galweight =  (double *)  malloc(Vipers_Num*sizeof(*clip_galweight));
    
    // Acceptanceflag set to zflag currently.
    for(j=0; j<Vipers_Num; j++){   
      fscanf(inputfile, "%d %lf %lf %lf %lf %lf\n", &id[j], &ra[j], &dec[j], &zobs[j], &TSR, &SSR);
      
      // ESR = TSR x SSR
      sampling[j] = TSR*SSR;

      clip_galweight[j] = 1.0;
    }

    for(j=0; j<10; j++)  printf("\n%d \t %le \t %le \t %le \t %le", id[j], ra[j], dec[j], zobs[j], sampling[j]);
    
    printf("\n\nVipers W%d catalogue input successful.", fieldFlag);
    
    fclose(inputfile);
    
    return 0;
}


int ParentInput_500s(char filepath[]){
  printf("\n\nOpening catalogue: %s", filepath);

  inputfile     = fopen(filepath, "r");

  if(inputfile == NULL){
    printf("\nError opening %s\n", filepath);
    return 1;
  }

  ch         = 0;
  Vipers_Num = 0;

  do{
    ch = fgetc(inputfile);
    if(ch == '\n')
      Vipers_Num += 1;
  } while (ch != EOF);

  rewind(inputfile);

  id             =  (int   *)   realloc(id, Vipers_Num*sizeof(*id));
  ra             =  (double *)  realloc(ra, Vipers_Num*sizeof(*ra));
  dec            =  (double *)  realloc(dec, Vipers_Num*sizeof(*dec));
  zobs           =  (double *)  realloc(zobs, Vipers_Num*sizeof(*zobs));
  M_B            =  (double *)  realloc(M_B, Vipers_Num*sizeof(*M_B));

  // derived parameters.                                                                                                                                                                                                                   
  Acceptanceflag =  (bool  *)   realloc(Acceptanceflag, Vipers_Num*sizeof(*Acceptanceflag));
  rDist          =  (double *)  realloc(rDist,          Vipers_Num*sizeof(*rDist));
  xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
  yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
  zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));
  sampling       =  (double *)  realloc(sampling,       Vipers_Num*sizeof(*sampling));
  fkp_galweight  =  (double *)  realloc(fkp_galweight,  Vipers_Num*sizeof(*fkp_galweight));
  clip_galweight =  (double *)  realloc(clip_galweight, Vipers_Num*sizeof(*clip_galweight));

  // redshift range 0.7<z<0.8, as traced by randoms. magnitude cut for known linear bias. volume limited to z=0.85                                                                                                                         
  for(j=0; j<Vipers_Num; j++){
    id[j] = j;

    fscanf(inputfile, "%*d \t %le \t %le \t %le \t %*le \t %*le \t %*le \t %*d \t %*d \t %*d \n", &ra[j], &dec[j], &zobs[j]);

    // add redshift errors. 
    zobs[j] += rand_gaussian(gsl_ran_r, 4.7*pow(10., -4.)*(1. + zobs[j]));

    sampling[j] = 1.;

    clip_galweight[j] = 1.;
  }

  fclose(inputfile);

  printf("\nHOD 500s catalogue input successful.");

  return 0;
}

/*
int CatalogueInput(char filepath[]){
    // Value added mocks input. //
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r"); 
     
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

    // Column  0: id                                                                                                         
    // Column  1: ra                                                                                           
    // Column  2: dec                                                                                               
    // Column  3: zcos from parent HOD catalogue.    
    // Column  4: zpec, cos + peculiar velocity.                                                                                                    
    // Column  5: zobs,  zpec + z_spec err. err_spec  = 0.0005*(1+z)                                                                          
    // Column  6: zphot, zpec + z_photo err, err_phot = 0.03*(1+z)
    // Column  7: MB, abs mag from the parent HOD. 
    // Column  8: type, 1:red, 2:blue. 
    // Column  9: csr, colour sampling rate. 
    // Column 10: sampling, TSRxSSR per quadrant.
    // Column 11: sampling35, TSRxSSR per quadrant, after diluting SSR to obtain an average TSRxSSR of 35%. 
    // Column 12: pointing, string.  Name of the pointing, e.g. W1P056
    // Column 13: quadrant, Q+number of quadrant. e.g. Q2.
    // Column 14: flag_inside,  1: inside Nagoya,   2: outside Nagoya. 
    // Column 15: flag_SSPOC,   1: inside SSR mock, 2: outside SSR mock. 
    // Column 16: flag_SSPOC35, 1: inside sample with global sampling rate of 35%, 2: all other galaxies. 
    // Column 17: rand_sel, random number from 0 to 1.  Assigned to all galaxies. 
 
    // Calculate number of Galaxies in the ZADE catalogue (line number);
    
    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    // ZADE Catalogue parameters.
    
    id             =  (int   *)   realloc(id,             Vipers_Num*sizeof(*id));
    ra             =  (double *)  realloc(ra,             Vipers_Num*sizeof(*ra));
    dec            =  (double *)  realloc(dec,            Vipers_Num*sizeof(*dec));
    zcos           =  (double *)  realloc(zcos,           Vipers_Num*sizeof(*zcos));
    zpec           =  (double *)  realloc(zpec,           Vipers_Num*sizeof(*zpec));
    zobs           =  (double *)  realloc(zobs,           Vipers_Num*sizeof(*zobs)); 
    zphot          =  (double *)  realloc(zphot,          Vipers_Num*sizeof(*zphot));
    M_B            =  (double *)  realloc(M_B,            Vipers_Num*sizeof(*M_B));
    type           =  (int   *)   realloc(type,           Vipers_Num*sizeof(*type));
    csr            =  (double *)  realloc(csr,            Vipers_Num*sizeof(*csr));
    sampling       =  (double *)  realloc(sampling,       Vipers_Num*sizeof(*sampling));
    sampling35     =  (double *)  realloc(sampling35,     Vipers_Num*sizeof(*sampling35));

    pointing       =  (char **)   realloc(pointing,       Vipers_Num*sizeof(char*));
    quadrant       =  (char **)   realloc(quadrant,       Vipers_Num*sizeof(char*));

    flag_Nagoya    =  (int   *)   realloc(flag_Nagoya,    Vipers_Num*sizeof(*flag_Nagoya));
    flag_SSPOC     =  (int   *)   realloc(flag_SSPOC,     Vipers_Num*sizeof(*flag_SSPOC));
    flag_SSPOC35   =  (int   *)   realloc(flag_SSPOC35,   Vipers_Num*sizeof(*flag_SSPOC35));
    rand_sel       =  (double *)  realloc(rand_sel,       Vipers_Num*sizeof(*rand_sel));
     
    // derived parameters. 
    Acceptanceflag =  (bool  *)   realloc(Acceptanceflag, Vipers_Num*sizeof(*Acceptanceflag));
    rDist          =  (double *)  realloc(rDist,          Vipers_Num*sizeof(*rDist));
    xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));
    
    
    for(j=0; j<1000; j++){
       pointing[j] =  (char *)    realloc(pointing[j], 20*sizeof(char));
       quadrant[j] =  (char *)    realloc(quadrant[j], 20*sizeof(char));
    }
    
    
    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf \t %lf \t %*s \t %*s \t %d \t %d \t %d \t %lf \n", &id[j], &ra[j], &dec[j], &zcos[j], &zpec[j], &zobs[j], &zphot[j], &M_B[j], &type[j], &csr[j], &sampling[j], &sampling35[j], &flag_Nagoya[j], &flag_SSPOC[j], &flag_SSPOC35[j], &rand_sel[j]);
    
    fclose(inputfile);
    
    printf("\nHOD catalogue input successful.");
    printf("\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    return 0;
}
*/
