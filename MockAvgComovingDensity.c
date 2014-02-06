int MockAverageComovingdensity(){
    float TotalW1W4area = 21.47; // sq. degs.

    // Equal intervals in comoving distance, for both W1 and W4 fields.  For all mocks.
    prepNumberDensityCalc();
    
    for(j=1; j<chiBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(0.4) + chiBinWidth*(j-1);

    for(forCount=1; forCount<CatalogNumber; forCount++){
        if(forCount < 10)   sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else                sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);
        
        // Choice of redshift from zcos, zpec, zphot, zobs.
        zUtilized       =   &zcos[0];
  
        assignAcceptance();

        // W1 catalogue.
        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<chiBinNumber; i++){
	            if((ChiSlices[i]<interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
                    NumberAtRedshift[i]   += 1;
                
                    MeanSliceRedshift[i]  += zUtilized[j];
            
                    SelectedGalNumber_nz  += 1;
                }
            }
        }
 
        // Now the W4 field. 
        if(forCount<10)   sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);

        // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
        zUtilized       =   &zcos[0];

        assignAcceptance();

        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<chiBinNumber; i++){
	            if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
	                NumberAtRedshift[i]   += 1;

	                MeanSliceRedshift[i]  += zUtilized[j];

	                SelectedGalNumber_nz  += 1;
	            }
            }
        }
    }
    
    for(j=1; j<chiBinNumber; j++) NumberAtRedshift[j]        /= CatalogNumber;
  
    for(j=1; j<chiBinNumber; j++) filteredNumberAtRedshift[j] = (float) NumberAtRedshift[j];
    
    for(j=1; j<chiBinNumber; j++) filtered_divfuncln_Atz[j]   = (float) log((float) NumberAtRedshift[j]/NChi_dChi(ChiSlices[j], chiBinWidth));
    
    Gaussianfilter(filteredNumberAtRedshift, ChiSlices, chiBinNumber, nzSigma);
    Gaussianfilter(filtered_divfuncln_Atz,   ChiSlices, chiBinNumber, nzSigma);
    
    for(j=1; j<chiBinNumber; j++) filtered_divfuncln_Atz[j]   = (float) NChi_dChi(ChiSlices[j], chiBinWidth)*exp(filtered_divfuncln_Atz[j]);

    for(j=1; j<chiBinNumber; j++){      // W1 + W4
        ComovingVolumeAtZ[j]            = sqdegs2steradians(TotalW1W4area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));

        ComovingNumberDensity[j]        = (float) (filteredNumberAtRedshift[j])/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j]       /= NumberAtRedshift[j];
        }
    }

    sprintf(filepath, "%s/Data/nz/HODMocks_MockAverage_nz_%.2f.dat", root_dir, chiBinWidth);
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++){
        fprintf(output, "%e \t %e \t %d \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j], filteredNumberAtRedshift[j], ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth), filtered_divfuncln_Atz[j]);
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
