 int ObservedDataInput(){
    printf("\n\nOpening catalogue: %s", filepath);
    
    sprintf(filepath, "/disk1/mjw/VIPERS/ZADE_VIPERSv3_w1_allphot.cat");
    
    inputfile     = fopen(filepath, "r");  
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

    // Column  0: id                                                                                                         
    // Column  1: ra                                                                                           
    // Column  2: dec                                                                                               
    // Column  3: zobs                                                                     
    // Column  4: zflag.
    // Column  5: MB, abs mag from the parent HOD. 
    // Column  6: sampling, TSRxSSR per quadrant.
    // Column  7: ZADE weight, TSRxSSR per quadrant.
    // Column  8: ZADE flag.

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
    zobs           =  (double *)  realloc(zobs,           Vipers_Num*sizeof(*zobs)); 
    zflag          =  (double *)  realloc(zflag,          Vipers_Num*sizeof(*zflag));
 
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
    
    for(j=0; j<Vipers_Num; j++){
        pointing[j] =  (char *) realloc(pointing[j], 20*sizeof(char));
        quadrant[j] =  (char *) realloc(quadrant[j], 20*sizeof(char));
    }

    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf \t %lf \t %s \t %s \t %d \t %d \t %d \t %lf \n", &id[j], &ra[j], &dec[j], &zcos[j], &zpec[j], &zobs[j], &zphot[j], &M_B[j], &type[j], &csr[j], &sampling[j], &sampling35[j], pointing[j], quadrant[j], &flag_Nagoya[j], &flag_SSPOC[j], &flag_SSPOC35[j], &rand_sel[j]);

    // Note: &pointing[j] must be passed to any printf statement. 
    
    fclose(inputfile);

    printf("\nHOD catalogue input successful.");
    printf("\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    return 0;
}


int CoordinateCalc(){
    for(j=0; j<Vipers_Num; j++){
            //  Derived parameters 
            polarAngle[j]         =  pi/2.0 - (pi/180.0)*dec[j];                // Converted to radians. 
            ra[j]                *= (pi/180.0);                                 // Converted to radians.

            //  Cosmology dependent, HOD mock parameters assumed - see header.h
            rDist[j]              = interp_comovingDistance(zUtilized[j]);      // Comoving distances in h^-1 Mpc
            
            xCoor[j]              = rDist[j]*sin(polarAngle[j])*cos(ra[j]);
            yCoor[j]              = rDist[j]*sin(polarAngle[j])*sin(ra[j]);
            zCoor[j]              = rDist[j]*cos(polarAngle[j]);
            ra[j]                /= (pi/180.0);                                 // Converted to degrees  
    }

    /*
    printf("\n\nOn input...");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    printf("\n\nRedshift max:  %f \t Redshift min:  %f", arrayMax(zUtilized, Vipers_Num), arrayMin(zUtilized, Vipers_Num));*/
    //printf("\n\nAbs. mag. max:  %f \t Abs. mag. min:  %f", arrayMax(M_B, Vipers_Num), arrayMin(M_B, Vipers_Num));

    return 0;
}
