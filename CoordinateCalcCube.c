int CoordinateCalcCube(){
    sprintf(filepath, "%s/Data/HODCube/cube_gal_-20.0.dat", root_dir);

    printf("\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

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
    } while (ch != EOF);

    printf("\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    rewind(inputfile);
    
    // ZADE Catalogue parameters.
    xCoor      =  (float *) realloc(xCoor,       Vipers_Num*sizeof(*xCoor));
    yCoor      =  (float *) realloc(yCoor,       Vipers_Num*sizeof(*yCoor));
    zCoor      =  (float *) realloc(zCoor,       Vipers_Num*sizeof(*zCoor));
    
    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%f \t %f \t %f \n", &xCoor[j], &yCoor[j], &zCoor[j]);
    
    fclose(inputfile);
   
    printf("\nCatalogue input successful.");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));

    return 0;
}