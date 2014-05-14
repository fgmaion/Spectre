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
    
    xCoor = (double *) realloc(xCoor, Vipers_Num*sizeof(*xCoor));
    yCoor = (double *) realloc(yCoor, Vipers_Num*sizeof(*yCoor));
    zCoor = (double *) realloc(zCoor, Vipers_Num*sizeof(*zCoor));

    xVel  = (double *) realloc(xVel,  Vipers_Num*sizeof(*xVel));
    yVel  = (double *) realloc(yVel,  Vipers_Num*sizeof(*yVel));
    zVel  = (double *) realloc(zVel,  Vipers_Num*sizeof(*zVel));

    for(j=0; j<Vipers_Num; j++){
        fscanf(inputfile, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n", &xCoor[j], &yCoor[j], &zCoor[j], &xVel[j], &yVel[j], &zVel[j]);
    }
    
    // for(j=0; j<10; j++)  printf("\n %le \t %le \t %le \t %le \t %le \t %le", xCoor[j], yCoor[j], zCoor[j], xVel[j], yVel[j], zVel[j]);

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
    
    sprintf(filepath, "%s/Data/HODCube/zcube_xvel_gal_-20.0.dat", root_dir);
    
    
    double xprime, yprime, zprime;
    
    /*
    // along the y axis direction. then relabel to x direction.  
    for(j=0; j<Vipers_Num; j++){
        yCoor[j]                      += (1. + 0.8)*yVel[j]/(100.*hz);
        
        if(yCoor[j] > 1000.) yCoor[j] -= 1000.;
        if(yCoor[j] <    0.) yCoor[j] += 1000.;
    
        xprime                         =   yCoor[j];
        yprime                         =  -xCoor[j] + 1000.;
        zprime                         =   zCoor[j];
    
        xCoor[j]                       =   xprime;
        yCoor[j]                       =   yprime;
        zCoor[j]                       =   zprime;
    }
    
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    
    sprintf(filepath, "%s/Data/HODCube/zcube_yvel_gal_-20.0.dat", root_dir);
    */
    /*
    // along the z axis direction. then relabel to x direction.  
    for(j=0; j<Vipers_Num; j++){
        zCoor[j]                      += (1. + 0.8)*zVel[j]/(100.*hz);
        
        if(zCoor[j] > 1000.) zCoor[j] -= 1000.;
        if(zCoor[j] <    0.) zCoor[j] += 1000.;
    
        xprime                         =   zCoor[j];
        yprime                         =   yCoor[j];
        zprime                         =  -xCoor[j] + 1000.;
    
        xCoor[j]                       =   xprime;
        yCoor[j]                       =   yprime;
        zCoor[j]                       =   zprime;
    }
    
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    
    sprintf(filepath, "%s/Data/HODCube/zcube_zvel_gal_-20.0.dat", root_dir);
    */
    
    output = fopen(filepath, "w");
    
    for(j=0; j<Vipers_Num; j++){
        fprintf(output, "%lf \t %lf \t %lf\n", xCoor[j], yCoor[j], zCoor[j]);
    }
    
    fclose(output);
    
    return 0;
}
