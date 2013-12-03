int PkCalc(){
    printf("\n\nAssigning FFT in array.");
    
    // Including the normalisation of the mask, TotalVolume/TotalSurveyedVolume, such that the filter fn. has unit mean."
    for(j=0; j<n0*n1*n2; j++) in[j][0] = densityArray[j]*FKPweights[j]*(TotalVolume/TotalSurveyedVolume)*booldensity[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");
    
    // 0: Subtract shot noise for a real survey, 1: Neglect shot noise subtraction for FFT of window function. 
    PkCorrections(0);
    
    PkBinningCalc(n0*n1*n2, PkArray);
    
    sprintf(filepath, "%s/Data/Del2k/midK_Del2k_%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], TotalVolume*binnedPk[j], modesPerBin[j], linearErrors[j]);
    fclose(output);
    
    return 0;
}


int PkCorrections(int WindowFuncParam){
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                GaussianFilter                     = exp(-1.*kSq*0.5*(pow(1., 2.)));

                WindowFunc                         = 1.;

                if(k_x != 0.){
		            WindowFunc                    *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);}
                
                if(k_y != 0.){
		            WindowFunc                    *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);}
                
                if(k_z != 0.){
		            WindowFunc                    *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);}

                H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
                H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];

                // The NGP window function is a real, scale dependent quantity.
                H_kReal                           /= WindowFunc;
                H_kImag                           /= WindowFunc;

                PkArray[Index][0]                  = pow(kSq, 0.5);
				
				// Subtract shot noise contribution for a real survey.  Neglect this for window function calculation.  FKP corrections. 
				if(WindowFuncParam == 0){
				    PkArray[Index][1]              = pow(H_kReal, 2.) + pow(H_kImag, 2.);
				    //PkArray[Index][1]           -= 1./TotalZADEWeight;
				    //PkArray[Index][1]           -= fkpShotNoiseCorr;
				    PkArray[Index][1]             *= pow(fkpWeightedVolume, 2.)*pow(fkpSqWeightsVolume, -1.)*pow(TotalVolume, -1.);
                }
                
                // PkCorrections called to correct NGP for survey window function.
                if(WindowFuncParam == 1){
                    PkArray[Index][1]              = pow(H_kReal, 2.) + pow(H_kImag, 2.);
	            }
	        }
        }
    }
    
    return 0;
}


int PkBinningCalc(int NumberModes, float Array[][2]){
    printf("\n\nPerforming binning calculation.");
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = j*kbinInterval;
    
    for(j=0; j<100; j++) printf("\n %f \t %f", Array[j][0], Array[j][1]);

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    printf("\nSorting mod k array.");
    qsort(&Array, NumberModes, sizeof(Array[0]), FirstColumnCompare);
    
    LowerBinIndex = 0;
    UpperBinIndex = 0;
    
    for(j=0; j<kBinNumb-1; j++){
        modesPerBin[j]  =              0;
        midKBin[j]      =            0.0;
        meanKBin[j]     =            0.0;
        binnedPk[j]     =            0.0;
        linearErrors[j] =            0.0;
    }
    
    for(j=0; j<kBinNumb-1; j++){
        for(i=LowerBinIndex; i<NumberModes; i++){
            if(Array[i][0] > kBinLimits[j+1]){
                UpperBinIndex = i;
                break;
            } 
        }

        for(i=LowerBinIndex; i<UpperBinIndex; i++){
            meanKBin[j]    += Array[i][0];
            binnedPk[j]    += Array[i][1];
            modesPerBin[j] += 1;
        }

	    if(modesPerBin[j]  != 0)  meanKBin[j]  /= modesPerBin[j];
        if(modesPerBin[j]  != 0)  binnedPk[j]  /= modesPerBin[j];

        // Peacock and Nicholson 1991, pg 313. above eqn (20).
        // Result of summing over a shell in k space containing m modes, should be a Gaussian random variable with variance 2.m/N^2  
              
        linearErrors[j]     = sqrt(2.*modesPerBin[j]/(TotalZADEWeight*TotalZADEWeight));

        del2[j]             = pow(meanKBin[j], 3.)*TotalVolume*binnedPk[j]*(4.*pi)/pow(2.*pi, 3.);
        
        LowerBinIndex       = UpperBinIndex;
        
        midKBin[j]          = 0.5*(kBinLimits[j] + kBinLimits[j+1]);
    }
    
    return 0;
}


int BinnedPkForMuInterval(float lowerMuLimit, float upperMuLimit, char filepath[], double Array[]){
    printf("\nPerforming P(k) binning for given mu interval, %f < mu < %f", lowerMuLimit, upperMuLimit);

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kmodulus                           = pow(pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.), 0.5);

                mu                                 = k_x/kmodulus;
                
                if(kmodulus == 0.)            mu   = 0.0;      

                if((lowerMuLimit < fabs(mu)) && (fabs(mu) < upperMuLimit)){
                    PkArray[Num_ModesInMuInterval][0]              = kmodulus;
                    PkArray[Num_ModesInMuInterval][1]              = Array[Index];
                    Num_ModesInMuInterval                         += 1;
                }
            }
        }
    }
    
    printf("\nNumber of modes in interval:  %d", Num_ModesInMuInterval);
    
    PkBinningCalc(Num_ModesInMuInterval, PkArray);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)     fprintf(output, "%d \t %g \t %g\n", modesPerBin[j], midKBin[j], binnedPk[j]);
    fclose(output);
    
    return 0;
}


int TwoDpkBinningCalc(){
    printf("\nPerforming 2d binning calc.");
    
    // Binned Pk assignment. 
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = j*kbinInterval;

    // Order by k_los and k_perp.
    printf("\nSorting 2d array.");
    qsort(TwoDpkArray, n0*n1*n2, sizeof(TwoDpkArray[0]), FirstColumnCompare);

    int losLowerBinIndex  = 0;
    int losUpperBinIndex  = 0;
    
    int perpLowerBinIndex = 0;
    int perpUpperBinIndex = 0;
    
    int m                 = 0;
    
    for(j=0; j<kBinNumb-1; j++){
        for(k=0; k<kBinNumb-1; k++){
            zSpaceBinnedPk[j][k]    = 0.0;
            mean_perpk[j][k]        = 0.0;
            mean_losk[j][k]         = 0.0;
            zSpacemodesPerBin[j][k] = 0;    
        }
    
    }
    
    for(j=0; j<kBinNumb-1; j++){
        for(i=losLowerBinIndex; i<n0*n1*n2; i++){
            if(TwoDpkArray[i][0] > kBinLimits[j+1]){
                losUpperBinIndex = i;
                break;
            } 
        }

        perpLowerBinIndex = losLowerBinIndex;

	qsort(&TwoDpkArray[losLowerBinIndex], losUpperBinIndex-losLowerBinIndex, sizeof(TwoDpkArray[0]), SecondColumnCompare);
       
        for(k=0; k<(kBinNumb-1); k++){
            for(m=perpLowerBinIndex; m<losUpperBinIndex; m++){
                if(TwoDpkArray[m][1] > kBinLimits[k+1]){
                    perpUpperBinIndex = m;
                    break;
                } 
            }
        
            for(m=perpLowerBinIndex; m<perpUpperBinIndex; m++){
                zSpaceBinnedPk[j][k]    += TwoDpkArray[m][2];
                mean_perpk[j][k]        += TwoDpkArray[m][1];
                mean_losk[j][k]         += TwoDpkArray[m][0];
                zSpacemodesPerBin[j][k] += 1;    
            }
            
            if(zSpacemodesPerBin[j][k] != 0)  zSpaceBinnedPk[j][k]  /= zSpacemodesPerBin[j][k];
            if(zSpacemodesPerBin[j][k] != 0)  mean_perpk[j][k]      /= zSpacemodesPerBin[j][k];
            if(zSpacemodesPerBin[j][k] != 0)  mean_losk[j][k]       /= zSpacemodesPerBin[j][k];
            
            perpLowerBinIndex = perpUpperBinIndex;
        }
        
        losLowerBinIndex = losUpperBinIndex;
    }

    return 0;
}
