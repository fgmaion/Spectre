int MockAverageComovingdensity(){
    // TotalObservedGalaxies    =  0.0;
    // dimmestAcceptedMagnitude = -99.;

    // Equal intervals in comoving distance, for both W1 and W4 fields.  For all mocks.
    prepNumberDensityCalc();
    
    for(j=0; j<chiBinNumber+1; j++){
       ChiSlices[j] = interp_comovingDistance(0.2) + chiBinWidth*j;
    }
    
    for(forCount=1; forCount<CatalogNumber+1; forCount++){
        if(forCount<10)     sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else                sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);
        
        // Choice of redshift from zcos, zpec, zphot, zobs.
        zUtilized       =   &zcos[0];
  
        // assignAcceptance();
        
        int BinnedGals = 0;

        // W1 catalogue.
        for(j=0; j<Vipers_Num; j++){
            for(i=0; i<chiBinNumber+1; i++){
	            if((ChiSlices[i]< interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) <= ChiSlices[i+1]) && (M_B[j]< absMagCut)){
                    NumberAtRedshift[i]   += 1.;
                    
                    BinnedGals            += 1;
                
                    MeanSliceRedshift[i]  += zUtilized[j];      
                }
           }
        }
        
        printf("\nVipers number: %d \t %d", Vipers_Num, BinnedGals);
        
        /*
        // Now the W4 field. 
        if(forCount<10)   sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, forCount);
        else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, forCount);

        CatalogueInput(filepath);

        // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
        zUtilized       =   &zcos[0];

        assignAcceptance();

        for(j=0; j<Vipers_Num; j++){
            for(i=1; i<chiBinNumber; i++){
	            if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1]) && (M_B[j]<absMagCut)){
	                NumberAtRedshift[i]   += 1.;

	                MeanSliceRedshift[i]  += zUtilized[j];    
	            }
            }
        }
        
        */
    }
    
    // Binning finished. Assignment for filtering. 
    
    for(j=0; j<chiBinNumber; j++){
        if(NumberAtRedshift[j] != 0){
             MeanSliceRedshift[j]        /= NumberAtRedshift[j];
        }
        
        else{
            MeanSliceRedshift[j]          = interp_inverseComovingDistance(ChiSlices[j] + 0.5*chiBinWidth); 
        } 
    }
    
    for(j=0; j<chiBinNumber; j++)  NumberAtRedshift[j]         /= CatalogNumber;
        
    // Preserve initial normalisation. 
    for(j=0; j<chiBinNumber; j++)  TotalObservedGalaxies       += (double) NumberAtRedshift[j];
    
    printf("\n\nTotal observed galaxies: %e", TotalObservedGalaxies);

    NormedGaussianfilter(&identitytransform, &identitytransform, NumberAtRedshift, ChiSlices, chiBinNumber, chiBinWidth, nzSigma, TotalObservedGalaxies,  filteredNumberAtRedshift);

    // NormedGaussianfilter(&ln_factor_fittedNz, &factor_fittedNz_exp, NumberAtRedshift, ChiSlices, chiBinNumber, chiBinWidth, nzSigma, TotalObservedGalaxies,  filtered_divfuncln_Atz);
    
    // Luminosity distance at redshift 0.9, (1+z)*R0*r. 
    // double Dl_z9;
    
    // Dl_z9        = (1. + 0.9)*interp_comovingDistance(0.9);
    
    for(j=0; j<chiBinNumber; j++){   // W1 + W4
        ComovingVolumeAtZ[j]                = sqdegs2steradians(W1area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));

        ComovingNumberDensity[j]            = (float) NumberAtRedshift[j]/ComovingVolumeAtZ[j];
        
        filteredComovingNumberDensity[j]    = filteredNumberAtRedshift[j]/ComovingVolumeAtZ[j];
        
        // LuminosityDistance[j]            = (1. + MeanSliceRedshift[j])*interp_comovingDistance(MeanSliceRedshift[j]);
    
        // Schechter_fn[j]                  = 0.0045*exp(-1.*(pow(LuminosityDistance[j]/Dl_z9, 2.) - 1.));
    }
    
    /*
    sprintf(filepath, "%s/Data/nz/HODMocks_MockAvgNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, chiBinWidth, absMagCut, nzSigma);
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++){
        fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j] + 0.5*chiBinWidth, NumberAtRedshift[j]/TotalW1W4area, filteredNumberAtRedshift[j]/TotalW1W4area, ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth)/TotalW1W4area, filtered_divfuncln_Atz[j]/TotalW1W4area);
    }
    
    fclose(output);
    */
    
    sprintf(filepath, "%s/Data/nz/HODMocks_MockAvg_nz_chiSliced_%.1f_W1_%.2f_%.1f.dat", root_dir, chiBinWidth, absMagCut, nzSigma);
    
    output = fopen(filepath, "w");

    for(j=0; j<chiBinNumber; j++){
        fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j] + 0.5*chiBinWidth, NumberAtRedshift[j], ComovingNumberDensity[j], filteredNumberAtRedshift[j], filteredComovingNumberDensity[j]);
    }
    
    fclose(output);
    
    return 0;
}


double MockAvg_nz(double Chi){
    float fChi = (float) Chi;
    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d, chiBinNumber-1, fChi, &fInterim);
    
    return (double) fInterim;
}


int splineMockAvg_nz(){
    prepNumberDensityCalc();

    // sprintf(filepath, "%s/Data/nz/HODMocks_MockAvgNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, chiBinWidth, absMagCut, nzSigma);
    sprintf(filepath, "%s/Data/nz/HODMocks_MockAvg_nz_chiSliced_%.1f_W1_%.2f_%.1f.dat", root_dir, chiBinWidth, absMagCut, nzSigma);
  
    inputfile     = fopen(filepath, "r"); 
    
    for(j=0; j<chiBinNumber; j++){
        // assignment suppression character, %*e.
        fscanf(inputfile, "%*e \t %e \t %*e \t %e %*e %*e \n", &ChiSlices[j], &ComovingNumberDensity[j]);
    }
    
    fclose(inputfile);
    
    spline(ChiSlices, ComovingNumberDensity, chiBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);
    
    
    sprintf(filepath, "%s/Data/splineChecks/HODMocks_MockAvg_nz_chiSliced_%.1f_W1_%.2f_%.1f.dat", root_dir, chiBinWidth, absMagCut, nzSigma);
    
    output = fopen(filepath, "w");
    
    double interimChi;
    
    for(j=0; j<2500; j++){
        interimChi = 500. + j*1.;
     
        fprintf(output, "%e \t %e \n", interimChi,  MockAvg_nz(interimChi));
    }
    
    fclose(output);
    
    return 0;    
}
