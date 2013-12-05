 int CoordinateCalc(char filepath[]){
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
    
    id           =  (int   *)  realloc(id,           Vipers_Num*sizeof(*id));
    ra           =  (float *)  realloc(ra,           Vipers_Num*sizeof(*ra));
    dec          =  (float *)  realloc(dec,          Vipers_Num*sizeof(*dec));
    zcos         =  (float *)  realloc(zcos,         Vipers_Num*sizeof(*zcos));
    zpec         =  (float *)  realloc(zpec,         Vipers_Num*sizeof(*zpec));
    zobs         =  (float *)  realloc(zobs,         Vipers_Num*sizeof(*zobs)); 
    zphot        =  (float *)  realloc(zphot,        Vipers_Num*sizeof(*zphot));
    M_B          =  (float *)  realloc(M_B,          Vipers_Num*sizeof(*M_B));
    type         =  (int   *)  realloc(type,         Vipers_Num*sizeof(*type));
    csr          =  (float *)  realloc(csr,          Vipers_Num*sizeof(*csr));
    sampling     =  (float *)  realloc(sampling,     Vipers_Num*sizeof(*sampling));
    sampling35   =  (float *)  realloc(sampling35,   Vipers_Num*sizeof(*sampling35));

    pointing     =  (float **) realloc(pointing,     Vipers_Num*sizeof(char*));
    quadrant     =  (float **) realloc(quadrant,     Vipers_Num*sizeof(char*));

    flag_Nagoya  =  (int   *)  realloc(flag_Nagoya,  Vipers_Num*sizeof(*flag_Nagoya));
    flag_SSPOC   =  (int   *)  realloc(flag_SSPOC,   Vipers_Num*sizeof(*flag_SSPOC));
    flag_SSPOC35 =  (int   *)  realloc(flag_SSPOC35, Vipers_Num*sizeof(*flag_SSPOC35));
    rand_sel     =  (float *)  realloc(rand_sel,     Vipers_Num*sizeof(*rand_sel));
     
    // derived parameters. 
    polarAngle   =  (float *)  realloc(polarAngle,   Vipers_Num*sizeof(*polarAngle));
    rDist        =  (float *)  realloc(rDist,        Vipers_Num*sizeof(*rDist));
    xCoor        =  (float *)  realloc(xCoor,        Vipers_Num*sizeof(*xCoor));
    yCoor        =  (float *)  realloc(yCoor,        Vipers_Num*sizeof(*yCoor));
    zCoor        =  (float *)  realloc(zCoor,        Vipers_Num*sizeof(*zCoor));
    
    for(j=0; j<Vipers_Num; j++){
        pointing[j] =  (float *) realloc(pointing[j], 20*sizeof(char));
        quadrant[j] =  (float *) realloc(quadrant[j], 20*sizeof(char));
    }

    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %d \t %f \t %f \t %f \t %s \t %s \t %d \t %d \t %d \t %f \n", &id[j], &ra[j], &dec[j], &zcos[j], &zpec[j], &zobs[j], &zphot[j], &M_B[j], &type[j], &csr[j], &sampling[j], &sampling35[j], pointing[j], quadrant[j], &flag_Nagoya[j], &flag_SSPOC[j], &flag_SSPOC35[j], &rand_sel[j]);

    // Note: &pointing[j] must be passed to any printf statement. 
    fclose(inputfile);

    printf("\nHOD catalogue input successful.");

    for(j=0; j<Vipers_Num; j++){
            //  Derived parameters 
            polarAngle[j]         =  pi/2.0 - (pi/180.0)*dec[j];                // Converted to radians. 
            ra[j]                *= (pi/180.0);                                 // Converted to radians.

            //  Cosmology dependent, HOD mock parameters assumed - see header.h
            rDist[j]              = interp_comovingDistance(zobs[j]);           // Comoving distances in h^-1 Mpc
            xCoor[j]              = rDist[j]*sin(polarAngle[j])*cos(ra[j]);
            yCoor[j]              = rDist[j]*sin(polarAngle[j])*sin(ra[j]);
            zCoor[j]              = rDist[j]*cos(polarAngle[j]);
            ra[j]                /= (pi/180.0);                                 // Converted to degrees  
    }

    printf("\n\nOn input...");
    printf("\nNumber of galaxies in catalogue:  %d", Vipers_Num);
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));

    printf("\n\nRedshift max:  %f \t Redshift min:  %f", arrayMax(zobs, Vipers_Num), arrayMin(zobs, Vipers_Num));
    return 0;
}


int VIPERSbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len){
  // Rotate the z Cartesian axis to line along the LOS.                                                                 

  float theta  = pi/2.0 - (pi/180.0)*centerDec;
  float phi    = centerRA*pi/180.0;

  float gamma1, gamma2, gamma3;

  for(j=0; j<len; j++){
    gamma1 =  xCoors[j]*sin(theta)*cos(phi)   + yCoors[j]*sin(theta)*sin(phi)   + zCoors[j]*cos(theta);
    gamma2 = -xCoors[j]*sin(phi)              + yCoors[j]*cos(phi)              + 0.;
    gamma3 = -xCoors[j]*cos(theta)*cos(phi)   - yCoors[j]*cos(theta)*sin(phi)   + zCoors[j]*sin(theta);

    xCoors[j] = gamma1;
    yCoors[j] = gamma2;
    zCoors[j] = gamma3;
  }

  return 0;
}


int Celestialbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len){
  // Rotate the z Cartesian axis to line along the LOS.                                                                                              
  float theta  = pi/2.0 - (pi/180.0)*centerDec;
  float phi    = centerRA*pi/180.0;

  float x, y, z;

  for(j=0; j<len; j++){
    x = xCoors[j]*sin(theta)*cos(phi) - yCoors[j]*sin(phi) - zCoors[j]*cos(theta)*cos(phi);
    y = xCoors[j]*sin(theta)*sin(phi) + yCoors[j]*cos(phi) - zCoors[j]*cos(theta)*sin(phi);
    z = xCoors[j]*cos(theta)          + 0.                 + zCoors[j]*sin(theta);

    xCoors[j] = x;
    yCoors[j] = y;
    zCoors[j] = z;
  }

  return 0;
}


// Returns comoving distance at redshift z in h^-1 Mpc. 
float interp_comovingDistance(float z){
    float InterimInterp_yVal;
    splint(z_Array, ComovingDistance_z, z_ComovingDistance_2derivatives, nPoints, z, &InterimInterp_yVal);
    return InterimInterp_yVal;
}


float interp_inverseComovingDistance(float r){
    float InterimInterp_yVal;
    splint(ComovingDistance_z, z_Array, ComovingDistance_z_2derivatives, nPoints, r, &InterimInterp_yVal);
    return InterimInterp_yVal;
}   // Returns z at comoving distance, [h^-1 Mpc]. 


float arrayMax(float a[], int n){
  float max = a[0];

  for(j=0; j< n; j++){
    if(a[j] > max){
      max = a[j]; 
    }
  }
  return max;
}


float arrayMin(float a[], int n){
  float min = a[0];

  for(j=0; j<n; j++){
    if(a[j] < min){
      min = a[j];
    }
  }
  return min;
}


double DoubleArrayMax(double a[], int n){
  double max = a[0];

  for(j=0; j< n; j++){
    if(a[j] > max){
      max = a[j]; 
    }
  }
  return max;
}


double DoubleArrayMin(double a[], int n){
  double min = a[0];

  for(j=0; j< n; j++){
    if(a[j] < min){
      min = a[j];
    }
  }
  return min;
}