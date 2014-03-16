int ComovingNumberDensityCalc_Data(){
    double  TotalObservedGalaxies    =    0.0;
    double  dimmestAcceptedMagnitude =  -99.0;

    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    prepNumberDensityCalc();
    
    for(j=1; j<chiBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(0.4) + chiBinWidth*(j-1);

    // Now the W1 field, for the real data. 
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CatalogueInput(filepath);

    // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
    zUtilized       =   &zcos[0];

    assignAcceptance();
    
    // Reinitialise.
    dimmestAcceptedMagnitude =  -99.0;

    // W1 catalogue.
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i]<interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1]) && (M_B[j]<absMagCut)){
                NumberAtRedshift[i]   += 1.;
                
                MeanSliceRedshift[i]  += zUtilized[j];
                
                if(M_B[j] > dimmestAcceptedMagnitude){
                    dimmestAcceptedMagnitude = M_B[j];
                }
            }
        }
    }

    printf("\nDimmest magnitude accepted in W1:  %le", dimmestAcceptedMagnitude);

    // Now the W4 field. 
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CatalogueInput(filepath);

    // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
    zUtilized       =   &zcos[0];

    assignAcceptance();
    
    // Reinitialise.
    dimmestAcceptedMagnitude =  -99.0;

    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1]) && (M_B[j]<absMagCut)){
	            NumberAtRedshift[i]   += 1.;

	            MeanSliceRedshift[i]  += zUtilized[j];
	            
	            if(M_B[j] > dimmestAcceptedMagnitude){
                    dimmestAcceptedMagnitude = M_B[j];
                }
	        }
        }
    }
    
    printf("\nDimmest magnitude accepted in W4:  %le", dimmestAcceptedMagnitude);

    for(j=1; j<chiBinNumber; j++)  TotalObservedGalaxies       += (double) NumberAtRedshift[j];
    
    NormedGaussianfilter(&identitytransform,  &identitytransform,   NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies, filteredNumberAtRedshift);
    
    NormedGaussianfilter(&ln_factor_fittedNz, &factor_fittedNz_exp, NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies, filtered_divfuncln_Atz);

    for(j=1; j<chiBinNumber; j++){      // W1 + W4
        ComovingVolumeAtZ[j]            = sqdegs2steradians(TotalW1W4area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));
        
        ComovingNumberDensity[j]        = (float) (filteredNumberAtRedshift[j])/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j]     /= NumberAtRedshift[j];
        }
    }
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/nz/HODMocks_00%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma);
    else              sprintf(filepath, "%s/Data/nz/HODMocks_0%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma); 
    
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++) fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j]/TotalW1W4area, filteredNumberAtRedshift[j]/TotalW1W4area, ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth), filtered_divfuncln_Atz[j]/TotalW1W4area);
    
    fclose(output);
    
    return 0;
}
