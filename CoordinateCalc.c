int CatalogueInput(char filepath[]){
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
    polarAngle     =  (double *)  realloc(polarAngle,     Vipers_Num*sizeof(*polarAngle));
    rDist          =  (double *)  realloc(rDist,          Vipers_Num*sizeof(*rDist));
    xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));
    
    /*
    for(j=0; j<1000; j++){
       pointing[j] =  (char *)    realloc(pointing[j], 20*sizeof(char));
       quadrant[j] =  (char *)    realloc(quadrant[j], 20*sizeof(char));
    }
    */
    
    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf \t %lf \t %*s \t %*s \t %d \t %d \t %d \t %lf \n", &id[j], &ra[j], &dec[j], &zcos[j], &zpec[j], &zobs[j], &zphot[j], &M_B[j], &type[j], &csr[j], &sampling[j], &sampling35[j], &flag_Nagoya[j], &flag_SSPOC[j], &flag_SSPOC35[j], &rand_sel[j]);
    /*
    for(j=0; j<Vipers_Num; j++){  
      fscanf(inputfile, "%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t \
%lf \t %lf \t %d \t %d \t %d \t %lf \n", &id[j], &ra[j], &dec[j], &zcos[j], &zpec[j], &zobs[j], &zphot[j], &M_B[j], &type[j], &csr[j], &sampling[j], &sampling35[j], &flag_Nagoya[j], &flag_SSPOC[j], &flag_SSPOC35[j], &rand_sel[j]);
    
      if(flag_Nagoya[j] == 1){
        printf("\n %e \t %e \t %e \t %d", ra[j], dec[j], zcos[j], flag_Nagoya[j]);
      }
    }
    */
    // Note: &pointing[j] must be passed to any printf statement. 
    
    fclose(inputfile);
    
    printf("\nHOD catalogue input successful.");
    printf("\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    return 0;
}


int CatalogueInput_500s(char filepath[]){
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    
    if(inputfile == NULL){  
        printf("\nError opening %s\n", filepath); 
        return 1;
    }

    // Column  0: id                                                                                                         
    // Column  1: ra                                                                                           
    // Column  2: dec                                                                                                                                                                
    // Column  3: zobs presumably.                                                                           
    // Column  4: MB, abs mag from the parent HOD. 
    // Column  5: apparent mags. 
    
    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    id             =  (int   *)   malloc(Vipers_Num*sizeof(*id));
    ra             =  (double *)  malloc(Vipers_Num*sizeof(*ra));
    dec            =  (double *)  malloc(Vipers_Num*sizeof(*dec));
    zobs           =  (double *)  malloc(Vipers_Num*sizeof(*zobs)); 
    zcos           =  (double *)  malloc(Vipers_Num*sizeof(*zcos)); 
    M_B            =  (double *)  malloc(Vipers_Num*sizeof(*M_B));

    // derived parameters. 
    Acceptanceflag =  (bool  *)   malloc(Vipers_Num*sizeof(*Acceptanceflag));
    polarAngle     =  (double *)  malloc(Vipers_Num*sizeof(*polarAngle));
    rDist          =  (double *)  malloc(Vipers_Num*sizeof(*rDist));
    xCoor          =  (double *)  malloc(Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  malloc(Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  malloc(Vipers_Num*sizeof(*zCoor));
    
    // redshift range 0.7<z<0.8, as traced by randoms. magnitude cut for known linear bias. volume limited to z=0.85
    for(j=0; j<Vipers_Num; j++){  
        id[j] = j;
    
        fscanf(inputfile, "%le \t %le \t %le \t %le \t %*le \t %*le \t %*le \t %*le \t %*le \t %*le \t %le \t %*le %*le \t %*d \t %*d \t %*d \n", &ra[j], &dec[j], &zobs[j], &zcos[j], &M_B[j]);
    
        // if((-20.2<M_B[j]) && (M_B[j] < -19.8) && (0.7<zobs[j]) && (zobs[j]<0.8))  Acceptanceflag[j] = true;
    }
    
    // for(j=0; j<10; j++)  printf("\n%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf", id[j], ra[j], dec[j], zobs[j], zcos[j], M_B[j]);
    
    fclose(inputfile);
    
    printf("\nHOD 500s catalogue input successful.");
    
    return 0;
}


double invert_StefanoBasis(double centreRA, double centreDec, double* xval, double* yval, double* zval){
    double x,  y,   z;
    double x1, y1, z1;
    double x2, y2, z2;
    
    double c_ra, c_dec;
    
    double rmod;
    
    x       = *xval;
    y       = *yval;
    z       = *zval;
    
    c_ra    =  centreRA*(pi/180.);
    c_dec   = centreDec*(pi/180.);

    // reverse translation. 
    x  -= stefano_trans_x;
    y  -= stefano_trans_y;
    z  -= stefano_trans_z;   
        
    // invert R2, x'' to x'
    x2         =  -sin(c_dec)*x + cos(c_dec)*z; 
    y2         =              y;
    z2         =  -cos(c_dec)*x - sin(c_dec)*z;
    
    // invert R1, x' to x.
    x1         =  cos(c_ra)*x2 - sin(c_ra)*y2;
    y1         =  sin(c_ra)*x2 + cos(c_ra)*y2;
    z1         =  z2;  
    
    // invert z inversion through xy plane. 
    *xval      =  x1;
    *yval      =  y1;
    *zval      = -z1; 
    
    rmod       = sqrt(x1*x1 + y1*y1 + z1*z1);

    return rmod;
}


int StefanoBasis(int Num, double ra[], double dec[], double rDist[], double xCoor[], double yCoor[], double zCoor[]){
    // Convert (ra, dec, z) into (x, y, z) for each catalogue in Stefano's basis. 

    for(j=0; j<Num; j++){
         ra[j]               *= (pi/180.0);                                 // Converted to radians.
        dec[j]               *= (pi/180.0);                                 // Converted to radians.
        
        rDist[j]              = interp_comovingDistance(gal_z[j]);          // Comoving distances in h^-1 Mpc
    
        xCoor[j]              = rDist[j]*cos(dec[j])*cos(ra[j]);        
        yCoor[j]              = rDist[j]*cos(dec[j])*sin(ra[j]);
        zCoor[j]              = rDist[j]*sin(dec[j]);                   // usual spherical co-ordinates.
        
        ra[j]                /= (pi/180.0);                                 // Converted to degrees.
        dec[j]               /= (pi/180.0);                                 // Converted to degrees.
    }
    
    
    printf("\n\nOn input...");
    printf("\nx max:  %.3f \t x min:  %.3f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %.3f \t y min:  %.3f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %.3f \t z min:  %.3f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    printf("\nr max:  %.3f \t r min:  %.3f", arrayMax(rDist, Vipers_Num), arrayMin(rDist, Vipers_Num));
    
    printf("\n\nRedshift max:  %f \t Redshift min:  %f", arrayMax(gal_z, Vipers_Num), arrayMin(gal_z, Vipers_Num));

    printf("\n\nAbs. mag. max:  %f \t Abs. mag. min:  %f", arrayMax(M_B, Vipers_Num), arrayMin(M_B, Vipers_Num));
    
    
    StefanoReflection(Vipers_Num, CentreRA, CentreDec, xCoor, yCoor, zCoor);
    
    // Rotate the input co-ordinates such that the ra direction is aligned more or less with the y axis, dec direction with x, and redshift along z. 
    StefanoRotated(Vipers_Num, CentreRA, CentreDec, xCoor, yCoor, zCoor);
    
    printf("\n\ninverted, rotated & translated");                                                                                                                                   
    printf("\nx max:  %.3f \t x min:  %.3f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %.3f \t y min:  %.3f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %.3f \t z min:  %.3f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    printf("\nr max:  %.3f \t r min:  %.3f", arrayMax(rDist, Vipers_Num), arrayMin(rDist, Vipers_Num));
                                                                                                                                                                                 
    printf("\n\naccepted. inverted, rotated & translated");                                                                  
    printf("\nx max:  %.3f \t x min:  %.3f", AcceptedMax(xCoor, Acceptanceflag, Vipers_Num), AcceptedMin(xCoor, Acceptanceflag, Vipers_Num));                       
    printf("\ny max:  %.3f \t y min:  %.3f", AcceptedMax(yCoor, Acceptanceflag, Vipers_Num), AcceptedMin(yCoor, Acceptanceflag, Vipers_Num));         
    printf("\nz max:  %.3f \t z min:  %.3f", AcceptedMax(zCoor, Acceptanceflag, Vipers_Num), AcceptedMin(zCoor, Acceptanceflag, Vipers_Num));                                     
    printf("\nr max:  %.3f \t r min   %.3f", AcceptedMax(rDist, Acceptanceflag, Vipers_Num), AcceptedMin(rDist, Acceptanceflag, Vipers_Num));                         
    
    return 0;
}


int StefanoReflection(int Number, double centreRA, double centreDec, double xCoors[], double yCoors[], double zCoors[]){
    // inversion through the xy plane.
    for(j=0; j<Number; j++){
        xCoors[j] *=  1.;
        yCoors[j] *=  1.;
        zCoors[j] *= -1.;   
    }

    return 0;
}


int StefanoRotated(int Number, double centreRA, double centreDec, double xCoors[], double yCoors[], double zCoors[]){
    double x1, y1, z1;
    double x2, y2, z2;
    double c_ra, c_dec;
    
    c_ra    =  centreRA*(pi/180.);
    c_dec   = centreDec*(pi/180.);

    // basis formed by: normal spherical co-ordinates subject to inversion through xy plane, then R1 and finally R2. 
    for(j=0; j<Number; j++){
        // R1: rotation about z such that in the new basis, (x',y',z'), x' hat lies in x-y plane at an angle centreRA to x.
        x1  =     cos(c_ra)*xCoors[j] + sin(c_ra)*yCoors[j];
        y1  =    -sin(c_ra)*xCoors[j] + cos(c_ra)*yCoors[j];
        z1  =               zCoors[j];
        
        // R2: rotation about y such that in the new basis, (x'', y'', z''), z'' hat lies in (x', z') plane at an angle -CentreDec to x' hat.
        x2  = -sin(c_dec)*x1  - cos(c_dec)*z1;
        y2  =  y1;
        z2  =  cos(c_dec)*x1  - sin(c_dec)*z1;
        
        // obsolete
        // x2  = cos(c_dec + pi/2.)*x1  - sin(c_dec + pi/2.)*z1;
        // y2  = y1;
        // z2  = sin(c_dec + pi/2.)*x1  + cos(c_dec + pi/2.)*z1;
        
        // finally translation in the box. P(k) unaffected. 
        xCoors[j] = x2 + stefano_trans_x;
        yCoors[j] = y2 + stefano_trans_y;
        zCoors[j] = z2 + stefano_trans_z;   
    }

    return 0;
}
