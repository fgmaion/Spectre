int CoordinateCalcCube(char filepath[]){
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    
    // Column 0: x co-ordinate                                                                                                     
    // Column 1: y co-ordinate                                                                                           
    // Column 2: z co-ordinate                                                                                             

    // Calculate number of in the catalogue (line number);

    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);     
           
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while(ch != EOF);
    
    // subsample the density to increase shot noise.
    Vipers_Num = (int) ceil(Vipers_Num/depletion_factor); 

    printf("\n\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    rewind(inputfile);
    
    // ZADE Catalogue parameters.
    xCoor          =  (double *)  realloc(xCoor,          Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  realloc(yCoor,          Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  realloc(zCoor,          Vipers_Num*sizeof(*zCoor));
    
    sampling       =  (double *)  realloc(sampling,       Vipers_Num*sizeof(*sampling));
    
    Acceptanceflag =  (bool  *)   realloc(Acceptanceflag, Vipers_Num*sizeof(*Acceptanceflag));
    fkp_galweight  =  (double *)  realloc(fkp_galweight,  Vipers_Num*sizeof(*fkp_galweight));
    clip_galweight =  (double *)  realloc(clip_galweight, Vipers_Num*sizeof(*clip_galweight));
    
    for(j=0; j<Vipers_Num; j++)   fscanf(inputfile, "%lf \t %lf \t %lf \n", &xCoor[j], &yCoor[j], &zCoor[j]);
    
    fclose(inputfile);
   
    printf("\n\nCatalogue input successful.");
    
    // for(j=0; j<Vipers_Num; j++)  xCoor[j] = 1000.*gsl_rng_uniform(gsl_ran_r);
    // for(j=0; j<Vipers_Num; j++)  yCoor[j] = 1000.*gsl_rng_uniform(gsl_ran_r);
    // for(j=0; j<Vipers_Num; j++)  zCoor[j] = 1000.*gsl_rng_uniform(gsl_ran_r);
    
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));

    // assign acceptance
    for(j=0; j<Vipers_Num; j++)  Acceptanceflag[j] = true;

    accepted_gals = 0;
    
    for(j=0; j<Vipers_Num; j++){
             sampling[j] = 1.;
        fkp_galweight[j] = 1.;
        
        if(Acceptanceflag[j] == true){
            accepted_gals += 1;
        }
    }
    
    // initialise.     
    for(j=0; j<Vipers_Num; j++) sampling[j] = 1.;
    
    return 0;
}
