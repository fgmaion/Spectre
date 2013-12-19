int MockAverageComovingdensity(){
    prepNumberDensityCalc();

    for(loopCount=1; loopCount<CatalogNumber; loopCount++){
        if(loopCount < 10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
        else                sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

        CatalogueInput(filepath);
        
        // Choice of redshift from zcos, zpec, zphot, zobs.
        zUtilized       =   &zcos[0];
  
        assignAcceptance();

        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<zBinNumber; i++){
	            if((redshiftSlices[i] < zUtilized[j]) && (zUtilized[j] < redshiftSlices[i+1])){
	                NumberAtRedshift[i]  += 1;
	                
	                MeanSliceRedshift[i] += zUtilized[j];
	            }
           }
       }
    }
  
   
    for(j=1; j<zBinNumber; j++){  
        ComovingVolumeAtZ[j]      = VIPERS_SolidAngle*pow(3., -1.)*(pow(interp_comovingDistance(redshiftSlices[j+1]), 3.) - pow(interp_comovingDistance(redshiftSlices[j]), 3.));
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j] /= NumberAtRedshift[j];
        }
        
        else{
            MeanSliceRedshift[j] = 0.5*(redshiftSlices[j+1] + redshiftSlices[j]); 
        }
        
        NumberAtRedshift[j]      /= CatalogNumber;
        
        ComovingNumberDensity[j]  = (float) NumberAtRedshift[j]/ComovingVolumeAtZ[j];
    }

    for(j=1; j<zBinNumber+1; j++) ChiSlices[j] = (float) interp_comovingDistance(MeanSliceRedshift[j]);

    sprintf(filepath, "%s/Data/nz/HODMocks_MockAverage_nz_%.2f.dat", root_dir, zBinWidth);
    output = fopen(filepath, "w");
    
    for(j=1; j<zBinNumber; j++){
        fprintf(output, "%e \t %e \t %d \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j], ComovingVolumeAtZ[j], ComovingNumberDensity[j], sdltNz(MeanSliceRedshift[j]), (float) NumberAtRedshift[j]/steradians2sqdegs(VIPERS_SolidAngle));
    }

    fclose(output);
   
    return 0;
}


int splineMockAvg_nz(){
    prepNumberDensityCalc();

    sprintf(filepath, "%s/Data/nz/HODMocks_MockAverage_nz_%.2f.dat", root_dir, zBinWidth);

    inputfile     = fopen(filepath, "r"); 
    
    for(j=1; j<zBinNumber; j++){
        // assignment suppression character, %*e.
        fscanf(inputfile, "%e \t %e \t %d \t %e \t %e \t %*e \t %*e\n", &MeanSliceRedshift[j], &ChiSlices[j], &NumberAtRedshift[j], &ComovingVolumeAtZ[j], &ComovingNumberDensity[j]);
    }
    
    fclose(inputfile);
    
    spline(ChiSlices, ComovingNumberDensity, zBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);
    
    return 0;    
}


double MockAvg_nz(double z){
    float fz = (float) z;
    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d, zBinNumber-1, fz, &fInterim);
    
    return (double) fInterim;
}
