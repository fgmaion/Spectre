int initialiseUnitMatrix(int n){
    unitMatrix = matrix(1, n, 1, n);

    for(i=1; i<n+1; i++){
        for(j=1; j<n+1; j++) unitMatrix[i][j] = 0.0;
    }

    for(j=1; j<n+1; j++)     unitMatrix[j][j] = 1.0;
    return 0;
}


double sdltNz_fitting(float array[], double z){
    // Numerical recipes indexing. 
    
    float alpha = array[1];
    float beta  = array[2];
    float z0    = array[3];
    float A     = array[4];

    // Base values
    // alpha = 8.603; 
    //  beta = 1.448;
    //    z0 = 0.191;
    //     A = 3.103;
    
    // Seemingly a 40% correction needs to be applied, eqn 2. from sdlt correlation fn. paper.
    return (100./40.)*zBinWidth*A*pow(z/z0, alpha)*exp(-1.*pow(z/z0, beta))*CSR(z);
}


// Chi squared. Numerical recipes indexing
float Chi2nz(float array[]){
    double Chi2 = 0.0;
    
    for(j=13; j<40; j++) weights[j] = fabs(1.0/sdltNz_fitting(array, xdata[j]));
    
    for(j=13; j<40; j++)  Chi2 += weights[j]*pow(ydata[j] - (*pt2Theoryfn)(array, xdata[j]), 2.);
    
    return (float) Chi2;
}


float redChi2(float array[]){
    return Chi2nz(array)/ (float) dof;
}


int quickMCfit(){
    mc_loopCount = 0;
    
    float          minChi2;
    
    for(j=13; j<40; j++) weights[j]  = fabs(1.0/sdltNz_fitting(nz_startParams, xdata[j]));
    
    minChi2                = redChi2(nz_startParams);
    
    float    fracPrecision = 10.;
    
    // while(minChi2 > 4.5){
    while(mc_loopCount < 1){
        if(mc_loopCount % 5000 == 0)  printf("\n\t\t %e \t %e \t %e \t %e \t %e", nz_startParams[1], nz_startParams[2], nz_startParams[3], nz_startParams[4], minChi2);

        // if(mc_loopCount % 5000 == 0)  printf("\n\t\t %e \t %e \t %e \t %e \t %e", stepArray[1], stepArray[2], stepArray[3], stepArray[4], minChi2);

        mc_jump(nz_startParams, &minChi2, 4, fracPrecision);
        
        mc_loopCount += 1;
    }
    
    return 0;
}


int mc_jump(float paramsArray[], float* minChi2, int paramNumber, float fracPrecision){
    for(i=1; i<paramNumber+1; i++){  
        for(k=0; k<paramNumber+1; k++) NewParams[k] = paramsArray[k];
    
        NewParams[i]      +=  (-0.5 + gsl_rng_uniform(gsl_rng_r))*stepArray[i];
        
        if(redChi2(NewParams) < *minChi2){
            *minChi2             = redChi2(NewParams);    
            
            // if(fabs(NewParams[i] - paramsArray[i]) > fabs(stepArray[i])){
            //    stepArray[i] *= 1.001;
            // }
        
        // *fracPrecision  = fabs(*minChi2/Chi2nz(NewParams) - 1.);
            paramsArray[i] = NewParams[i];    
        }
    }

    return 0;
}

int fitnz(){
    pt2Theoryfn    = &sdltNz_fitting;
    
    dof            = (39 - 13) - 4;
    
    stepArray      = malloc(5*sizeof(*stepArray));
    nz_startParams = malloc(5*sizeof(*nz_startParams));
    NewParams      = malloc(5*sizeof(*NewParams));
    
    stepArray[0]   = 0.0;
    stepArray[1]   = 0.1;
    stepArray[2]   = 0.1;
    stepArray[3]   = 0.05;
    stepArray[4]   = 4.0;
    
    // Base values:  alpha = 8.603, beta = 1.448, z0 = 0.191, A = 35.*3.103

/*
    // chi2: 5.18
    nz_startParams[0] =  0.0;
    nz_startParams[1] =  6.986100;
    nz_startParams[2] =  1.303534;
    nz_startParams[3] =  0.1792731;
    nz_startParams[4] =  257.1606;
*/

  // red chi2: 5.03
    nz_startParams[0] =  0.0;
    nz_startParams[1] =  7.041953;
    nz_startParams[2] =  1.253965;
    nz_startParams[3] =  0.1632144;
    nz_startParams[4] =  165.2622;

    printf("\n\nStarting params:\n%e \t %e \t %e \t %e \t %e\n\n", nz_startParams[0], nz_startParams[1], nz_startParams[2], nz_startParams[3], nz_startParams[4]);

    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/nz/HODMocks_MockAvg_nz_%.2f.dat", zBinWidth);
    inputfile     = fopen(filepath, "r");
    
    ch            = 0;
    len           = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  len += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    xdata          = malloc(len*sizeof(*xdata));
    ydata          = malloc(len*sizeof(*ydata));
    edata          = malloc(len*sizeof(*edata));
    weights        = malloc(len*sizeof(*weights));
    
    zdata          = malloc(len*sizeof(*zdata));
    Nzdata         = malloc(len*sizeof(*Nzdata));    
    
    for(j=0; j<len; j++){
        fscanf(inputfile, "%lg \t %*lg \t %d \t %*lg \t %*lg \t %*lg \t %*lg\n", &zdata[j], &Nzdata[j]);
        
	   xdata[j]    = zdata[j];
	   
	   // Unit solid angle conversion.
	   ydata[j]    = (double) Nzdata[j]/steradians2sqdegs(VIPERS_SolidAngle);
	   
	   edata[j]    = 0.0;
     weights[j]    = 1.0;
    }

    fclose(inputfile);

    // initialiseUnitMatrix(4);

    /* 
        Minimization of a function func of n variables.  Input consists of an initial starting point p[1..n];
        an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions (usually the n 
        unit vectors); and ftol, the fractional tolerance in the function value such that failure to decrease
        by more than this amount in one iteration signals completion.  On output, p is set to the best point 
        found, xi is the then-current direction set, fret is the returned function value at p, and iter is the
        number of iterations taken.  The routine linmin is used. 
    */
    
    // powell(nz_startParams, unitMatrix, 4, 0.0001, &maxiter, &returnval, &redChi2);
    
    gsl_rng_env_setup();
    
    gsl_rng_T = gsl_rng_default;
    gsl_rng_r = gsl_rng_alloc(gsl_rng_T);
    
    quickMCfit();
    
    gsl_rng_free(gsl_rng_r);
    
    // for(j=0; j<len; j++)  printf("\n%e \t %e \t %e \t %e", xdata[j], ydata[j], weights[j], sdltNz_fitting(nz_startParams, xdata[j]));
    
    printf("\n\nChi^2: %e", Chi2nz(nz_startParams));
    printf("\nreduced Chi^2: %e", redChi2(nz_startParams));

    printf("\n\nBest fitting N(z) parameters:\n%e \t %e \t %e \t %e \t %e\n\n", nz_startParams[0], nz_startParams[1], nz_startParams[2], nz_startParams[3], nz_startParams[4]);

    sprintf(filepath, "%s/Data/nz/HODMocks_MinChi2_MockAvgNz_%.2f.dat", root_dir, zBinWidth);
    output = fopen(filepath, "w");
    
    for(j=13; j<40; j++){
        fprintf(output, "%e \t %e \t %e \n", xdata[j], ydata[j], sdltNz_fitting(nz_startParams, xdata[j]));
    }
    
    fclose(output);
    return 0;
}
