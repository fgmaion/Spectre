int MockAverageComovingdensity(){
    double  TotalObservedGalaxies =   0.0;

    // Equal intervals in comoving distance, for both W1 and W4 fields.  For all mocks.
    prepNumberDensityCalc();
    
    for(j=1; j<chiBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(0.4) + chiBinWidth*(j-1);

    for(forCount=1; forCount<CatalogNumber; forCount++){
        if(forCount < 10)   sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else                sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);
        
        // CoordinateCalc();
        
        // Choice of redshift from zcos, zpec, zphot, zobs.
        zUtilized       =   &zcos[0];
  
        assignAcceptance();

        // W1 catalogue.
        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<chiBinNumber; i++){
	            if((ChiSlices[i]<interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
                    NumberAtRedshift[i]   += 1.;
                
                    MeanSliceRedshift[i]  += zUtilized[j];
            
                    SelectedGalNumber_nz  += 1;
                }
            }
        }
 
        // Now the W4 field. 
        if(forCount<10)   sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);
        
        // CoordinateCalc();

        // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
        zUtilized       =   &zcos[0];

        assignAcceptance();

        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<chiBinNumber; i++){
	            if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
	                NumberAtRedshift[i]   += 1.;

	                MeanSliceRedshift[i]  += zUtilized[j];

	                SelectedGalNumber_nz  += 1;
	            }
            }
        }
    }
    
    // Binning finished. Assignment for filtering. 
    for(j=1; j<chiBinNumber; j++)  NumberAtRedshift[j]         /= CatalogNumber;
        
    // Preserve initial normalisation. 
    for(j=1; j<chiBinNumber; j++)  TotalObservedGalaxies       += (double) NumberAtRedshift[j];
    
    NormedGaussianfilter(&identitytransform, &identitytransform,    NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies,  filteredNumberAtRedshift);

    NormedGaussianfilter(&ln_factor_fittedNz, &factor_fittedNz_exp, NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies,  filtered_divfuncln_Atz);
    
    for(j=1; j<chiBinNumber; j++){   // W1 + W4
        ComovingVolumeAtZ[j]         = sqdegs2steradians(TotalW1W4area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));

        ComovingNumberDensity[j]     = (float) (filteredNumberAtRedshift[j])/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j]    /= NumberAtRedshift[j];
        }
    }

    sprintf(filepath, "%s/Data/nz/HODMocks_MockAverage_nz_%.2f.dat", root_dir, chiBinWidth);
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++){
        fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j]/TotalW1W4area, filteredNumberAtRedshift[j]/TotalW1W4area, ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth)/TotalW1W4area, filtered_divfuncln_Atz[j]/TotalW1W4area);
    }
    
    fclose(output);
    
    return 0;
}


int splineMockAvg_nz(){
    prepNumberDensityCalc();

    sprintf(filepath, "%s/Data/nz/HODMocks_MockAverage_nz_%.2f.dat", root_dir, chiBinWidth);

    inputfile     = fopen(filepath, "r"); 
    
    for(j=1; j<chiBinNumber; j++){
        // assignment suppression character, %*e.
        fscanf(inputfile, "%e \t %e \t %*e \t %*e \t %e \t %*e \t %*e\n", &MeanSliceRedshift[j], &ChiSlices[j], &ComovingNumberDensity[j]);
    }
    
    fclose(inputfile);
    
    spline(ChiSlices, ComovingNumberDensity, chiBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);
    
    return 0;    
}


double MockAvg_nz(double Chi){
    float fChi = (float) Chi;
    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d, chiBinNumber-1, fChi, &fInterim);
    
    return (double) fInterim;
}
