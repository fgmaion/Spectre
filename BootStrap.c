int prepBootStrap(int objectNumber, double* xVal, double* yVal, double* zVal, double CubeSize){
    BootStrap_flag = realloc(BootStrap_flag, objectNumber*sizeof(int));
    BootStrap_Wght = realloc(BootStrap_Wght, objectNumber*sizeof(double));

    // for(j=0; j<objectNumber; j++)  BootStrap_flag[j] = SubVolAssign(xVal[j], yVal[j], zVal[j], CubeSize);
    
    for(j=0; j<objectNumber; j++)  BootStrap_Wght[j] = 1.0;

    return 0;
}


int BootStrapGen(int objectNumber, double* xVal, double* yVal, double* zVal, double CubeSize){
    // flag denoting subsample in which the galaxy resides. 

    for(j=0; j<objectNumber; j++)  BootStrap_Wght[j] = 0.0;

    // xBoot     = realloc(xBoot, 2*objectNumber*sizeof(double));
    // yBoot     = realloc(yBoot, 2*objectNumber*sizeof(double));
    // zBoot     = realloc(zBoot, 2*objectNumber*sizeof(double));
    
    // deltaBoot = realloc(deltaBoot, 2*objectNumber*sizeof(double));
    
    // for(j=0; j<n0*n1*n2; j++) deltaBoot[j] = densityArray[j]; 

    // Number of sub volumes.
    int C =  8;
    
    int sample;
    
    // int BootStrap_Count = 0;
    
    // Reweight original catalogue in the BootStrap method, with Nr = Nsub.
    printf("\n\nBoot strapping catalogue.\n\n");
    
    // int printGal;
    
    // oversampling
    for(k=0; k<3*C; k++){
        sample = rand() % C;
    
        printf("%d \n", sample);
    
        // printGal = 0;
        
        for(j=0; j<objectNumber; j++){
            if(BootStrap_flag[j] == sample){
                // 
	            // xBoot[BootStrap_Count] = fmod(xVal[j], CubeSize/2.) + 0.5*CubeSize*(k%2); 
	            // yBoot[BootStrap_Count] = fmod(yVal[j], CubeSize/2.) + 0.5*CubeSize*((k-k%2)%4)/2; 
	            // zBoot[BootStrap_Count] = fmod(zVal[j], CubeSize/2.) + 0.5*CubeSize*(k-(k%2)-((k-k%2)%4))/4; 
                
                // ii = (int) floor(fmod(xVal[j], CubeSize/2.)/CellSize);
                // jj = (int) floor(fmod(yVal[j], CubeSize/2.)/CellSize);
                // kk = (int) floor(fmod(zVal[j], CubeSize/2.)/CellSize);
                
                // boxlabel = ((k%2) + ii) + (((k-k%2)%4)/2  +  jj)*n2 + ((k-(k%2)-((k-k%2)%4))/4 + kk)*n1*n2;
                
                // densityArray[boxlabel] = densityArray[j];
                
                BootStrap_Wght[j] += 1.0;
                
                /*
                if(printGal<10){
                    printf("\n %d \t %d \t %e \t %e \t %e \t %e \t %e \t %e", BootStrap_flag[j], k, xVal[j],  xBoot[BootStrap_Count], yVal[j], yBoot[BootStrap_Count], zVal[j], zBoot[BootStrap_Count]);
    
                    printGal += 1;
                }
                */
                
                // BootStrap_Count                 += 1;
            }
        }
        
        // Note due to sample/Cosmic variance number of galaxies in boot strap catalogue need not match number
        // in input catalogue.
    }

    // Re-assignment of number of galaxies in original catalogue to number in boot strap catalogue.
    /*
    objectNumber = BootStrap_Count;
    
    printf("\n\nNumber of galaxies in boot strap catalogue: %d", BootStrap_Count);
    
    xVal           =  (double *)  realloc(xVal,          objectNumber*sizeof(*xVal));
    yVal           =  (double *)  realloc(yVal,          objectNumber*sizeof(*yVal));
    zVal           =  (double *)  realloc(zVal,          objectNumber*sizeof(*zVal));
    
    Acceptanceflag =  (bool  *)   realloc(Acceptanceflag, objectNumber*sizeof(*Acceptanceflag));
    
    for(j=0; j<objectNumber; j++){
        xVal[j] = xBoot[j]; 
        yVal[j] = yBoot[j]; 
        zVal[j] = zBoot[j]; 
    }
    */
    // for(j=0; j<10; j++) printf("\n %e", xVal[j]);
    
    return 0;
}


int SubVolAssign(double x, double y, double z, double CubeSize){
    int boxlabel;
    
    int xlow, ylow, zlow;
    
    xlow = ylow = zlow = 1;
    
    if(x<=CubeSize/2.) xlow = 0;
    if(y<=CubeSize/2.) ylow = 0;
    if(z<=CubeSize/2.) zlow = 0;
    
    boxlabel = xlow + 2*ylow + 2*2*zlow;    
    
    return boxlabel;
}
