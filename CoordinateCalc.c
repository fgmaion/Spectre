int prep_CatalogueInput_500s(){
  // Maximum number of galaxies present in any mock of the collection (i.e. those for covariance estimate).
  if(strcmp(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W1") == 0)       max_gals = 61765;      // Data: 50829
  else if(strcmp(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W4") == 0)  max_gals = 30261;      //     : 26377  
  else    max_gals = max_gal();

  // Initialise to max. memory required to process all mocks.
  ra               =  (double *)  calloc(max_gals, sizeof(*ra));
  dec              =  (double *)  calloc(max_gals, sizeof(*dec));
  zobs             =  (double *)  calloc(max_gals, sizeof(*zobs));
  M_B              =  (double *)  calloc(max_gals, sizeof(*M_B));
  Acceptanceflag   =  (bool   *)  calloc(max_gals, sizeof(*Acceptanceflag));
  rDist            =  (double *)  calloc(max_gals, sizeof(*rDist));
  xCoor            =  (double *)  calloc(max_gals, sizeof(*xCoor));
  yCoor            =  (double *)  calloc(max_gals, sizeof(*yCoor));
  zCoor            =  (double *)  calloc(max_gals, sizeof(*zCoor));
  sampling         =  (double *)  calloc(max_gals, sizeof(*sampling));
  fkp_galweight    =  (double *)  calloc(max_gals, sizeof(*fkp_galweight));
  clip_galweight   =  (double *)  calloc(max_gals, sizeof(*clip_galweight));

  gal_z = &zobs[0];  // Choice of redshift from zcos, zpec, zphot, zobs.       

  return 0;
}


int spec_weights(){
  // DT TSR, for new mocks in W1_Spectro_V7_2.
  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/TSR/TSR_W%d_mock_%03d_Nagoya_v6_Samhain.dat", root_dir, fieldFlag, loopCount);

  inputfile = fopen(filepath, "r");

  for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%*d \t %*f \t %*f \t %lf \n", &sampling[j]);

  fclose(inputfile);

  return 0;
}


int CatalogueInput_500s(){  
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile   = fopen(filepath, "r");

    line_count(inputfile, &Vipers_Num);
    
    for(j=0; j<Vipers_Num; j++){  
      fscanf(inputfile, "%*d \t %le \t %le \t %le \t %le \t %*e \t %*e \t %*d \t %*d \t %*d \n", &ra[j], &dec[j], &zobs[j], &M_B[j]);
      
       ra[j]           -=   CentreRA;                                 // R1 rotation of stefano basis, remains to rotate by R2. 
      
       ra[j]           *= (M_PI/180.0);                               // Convert to radians.
      dec[j]           *= (M_PI/180.0);                               // Convert to radians. 

      rDist[j]          = interp_comovingDistance(gal_z[j]);          // Comoving distances in h^-1 Mpc
      
      xCoor[j]          =  rDist[j]*cos(dec[j])*cos(ra[j]);
      yCoor[j]          =  rDist[j]*cos(dec[j])*sin(ra[j]);
      zCoor[j]          = -rDist[j]*sin(dec[j]);                       // reflection of spherical coordinates. 

       ra[j]           /= (M_PI/180.0);                               // Converted to degrees.
      dec[j]           /= (M_PI/180.0);                               // Converted to degrees.  
    }
    
    fclose(inputfile);
    
    // Load ESR weights.  
    spec_weights();    // load sampling according to local TSR.
    
    return 0;
}


int DataInput(){
  double TSR, SSR;

  sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/vipers_lss_v7_W%d.txt", root_dir, fieldFlag);

  printf("\n\nOpening catalogue: %s\n", filepath);
  
  inputfile = fopen(filepath, "r"); 

  line_count(inputfile, &Vipers_Num);
    
  for(j=0; j<Vipers_Num; j++){   
    fscanf(inputfile, "%*d %lf %lf %lf %lf %lf\n", &ra[j], &dec[j], &zobs[j], &TSR, &SSR);

    sampling[j] = TSR*SSR;  // ESR = TSR x SSR  
    
    ra[j]      -= CentreRA; // R1 rotation of stefano basis, remains to rotate by R2.
      
    ra[j]      *= (M_PI/180.0);                               // Convert to radians.
    dec[j]     *= (M_PI/180.0);                               // Convert to radians.

    rDist[j]    = interp_comovingDistance(gal_z[j]);          // Comoving distances in h^-1 Mpc

    xCoor[j]    =  rDist[j]*cos(dec[j])*cos(ra[j]);
    yCoor[j]    =  rDist[j]*cos(dec[j])*sin(ra[j]);
    zCoor[j]    = -rDist[j]*sin(dec[j]);                      // reflection of spherical coordinates.

    ra[j]      /= (M_PI/180.0);                               // Converted to degrees.
    dec[j]     /= (M_PI/180.0);                               // Converted to degrees.  
  }
    
  fclose(inputfile);
    
  return 0;
}
