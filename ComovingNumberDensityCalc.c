int ComovingNumberDensityCalc(){
    zBinWidth                =                                     0.1;
    zBinNumber               =  (int)     floor((1.4 - 0.4)/zBinWidth);
    
    // Numerical recipes indexing for spline. arrays begin at a[1].
    
    redshiftSlices           =  (float *) realloc(redshiftSlices,          (zBinNumber+1)*sizeof(*redshiftSlices));
    ChiSlices                =  (float *) realloc(ChiSlices,               (zBinNumber+1)*sizeof(*ChiSlices));
    NumberAtReshift          =  (float *) realloc(NumberAtReshift,          zBinNumber*sizeof(*NumberAtReshift));
    ComovingNumberDensity    =  (float *) realloc(ComovingNumberDensity,    zBinNumber*sizeof(*ComovingNumberDensity));
    ComovingVolumeAtZ        =  (float *) realloc(ComovingVolumeAtZ,        zBinNumber*sizeof(*ComovingVolumeAtZ));

    // Second derivatives. 
    ComovingNumberDensity2d  =  (float *) realloc(ComovingNumberDensity2d,  zBinNumber*sizeof(*ComovingNumberDensity2d));

    for(j=1; j<zBinNumber+1;   j++)    redshiftSlices[j]   = 0.4 + (j-1)*zBinWidth;
    for(j=1; j<zBinNumber; j++)        NumberAtReshift[j]  = 0.0;
    
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<zBinNumber; i++){
            if((redshiftSlices[i] < zUtilized[j]) && (zUtilized[j] < redshiftSlices[i+1]) && (M_B[j] < absMagCut)){
                NumberAtReshift[i] += 1;
            }
        }
    }

    for(j=1; j<zBinNumber; j++){
        ComovingVolumeAtZ[j]     = VIPERS_SolidAngle*pow(3., -1.)*(pow(interp_comovingDistance(redshiftSlices[j+1]), 3.) - pow(interp_comovingDistance(redshiftSlices[j]), 3.));
        ComovingNumberDensity[j] = NumberAtReshift[j]/ComovingVolumeAtZ[j];
    }

    sprintf(filepath, "%s/Data/nz/HODMocks_nz.dat", root_dir);
    output = fopen(filepath, "w");
    
    for(j=1; j<zBinNumber; j++){
        fprintf(output, "%g \t %g \t %g \t %g \n", redshiftSlices[j], NumberAtReshift[j], ComovingVolumeAtZ[j], ComovingNumberDensity[j]);
    }

    fclose(output);

    for(j=1; j<zBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(redshiftSlices[j]);

    spline(ChiSlices, ComovingNumberDensity, zBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);

    return 0;
}


float interp_nz(float x){
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d, zBinNumber-1, x, &Interim);
    return Interim;
}


