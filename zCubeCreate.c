int zCubeCreate(){
    sprintf(filepath, "%s/Data/HODCube/cube_gal_-20.0_vel.dat", root_dir);
    
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

    // Column  0: x coordinate.                                                                                                         
    // Column  1: y coordinate.                                                                                         
    // Column  2: z coordinate.                                                                                                
    // Column  3: x peculiar vel.    
    // Column  4: y peculiar vel.                                                                               
    // Column  5: z peculiar vel.

    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while (ch != EOF);

    printf("\n\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    rewind(inputfile);
    
    xCoor = (float *) realloc(xCoor, Vipers_Num*sizeof(*xCoor));
    yCoor = (float *) realloc(yCoor, Vipers_Num*sizeof(*yCoor));
    zCoor = (float *) realloc(zCoor, Vipers_Num*sizeof(*zCoor));

    xVel  = (float *) realloc(xVel, Vipers_Num*sizeof(*xVel));
    yVel  = (float *) realloc(yVel, Vipers_Num*sizeof(*yVel));
    zVel  = (float *) realloc(zVel, Vipers_Num*sizeof(*zVel));


    for(j=0; j<Vipers_Num; j++){
        fscanf(inputfile, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", &xCoor[j], &yCoor[j], &zCoor[j], &xVel[j], &yVel[j], &zVel[j]);
    }

    printf("\nInput of catalogue for z cube creation successful.");

    // R_0 dr = (c/H(z)) dz = dv/H(z), given that 1+z = 1+(v/c).  
    //        = dv [km s^-1]/ (100 h(z) [h]) h^-1 Mpc
    
    // Here we presume the cube is at redshift 0.8

    hz = pow(Om_v + Om_m*pow(1. + 0.8, 3.), 0.5); // Units of h, [h].
    printf("\nh(z) at redshift 0.8:  %f", hz);
    
    // Add reshift space distortion in the plane parallel approximation, along the x axis direction. 
    
    for(j=0; j<Vipers_Num; j++){
        xCoor[j]                      += (1. + 0.8)*xVel[j]/(100.*hz);
        
        if(xCoor[j] > 1000.) xCoor[j] -= 1000.;
        if(xCoor[j] <    0.) xCoor[j] += 1000.;
    }
    
    sprintf(filepath, "%s/Data/HODCube/zcube_gal_-20.0_2.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<Vipers_Num; j++){
        fprintf(output, "%f \t %f \t %f\n", xCoor[j], yCoor[j], zCoor[j]);
    }
    
    fclose(output);

    return 0;
}
