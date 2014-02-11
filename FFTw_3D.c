int PkCalc(){
    printf("\n\nAssigning FFT in array.");
    
    // The true density field multiplied by a mask, W(x). 
    for(j=0; j<n0*n1*n2; j++) in[j][0] = densityArray[j]*Cell_AppliedWindowFn[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");
    
    // 0: Subtract shot noise for a real survey, 1: Neglect shot noise subtraction for FFT of window function. 
    PkCorrections(0);
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/Del2k/midK_Del2k_%s_00%d.dat", root_dir, surveyType, loopCount);
    else              sprintf(filepath, "%s/Data/Del2k/midK_Del2k_%s_0%d.dat", root_dir, surveyType, loopCount);
    Monopole(filepath);

    // 2D Power spectrum.
    
    // if(loopCount<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_00%d.dat", root_dir, surveyType, loopCount);
    // else              sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_0%d.dat", root_dir, surveyType, loopCount);
    // Cartesian2Dpk();
    
    // if(loopCount<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observedpolar2Dpk_%s_00%d.dat", root_dir, surveyType, loopCount);
    // else              sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observedpolar2Dpk_%s_0%d.dat", root_dir, surveyType, loopCount);
    // polar2DpkCalc(filepath);
    
    // Quadrupole.
    // MultipoleCalc(2, kBinNumb, kQuadrupole);
    
    return 0;
}


int Monopole(char filepath[]){
    // Peacock 2dF:                 beta = 0.43 +- 0.07
    
    PkBinningCalc(n0*n1*n2, PkArray);

    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %e \t %d \t %e \n", meanKBin[j], del2[j], TotalVolume*binnedPk[j], modesPerBin[j], linearErrors[j]);

    fclose(output);

    return 0;
}


int Cartesian2Dpk(char filepath[]){
    for(j=0;  j<loskBinNumb;  j++)   loskBinLimits[j]  = j*2.*kbinInterval;
    for(j=0; j<perpkBinNumb;  j++)   perpkBinLimits[j] = loskBinLimits[j];
    
    DualBinning(n0*n1*n2, loskBinNumb, loskBinLimits, perpkBinNumb, perpkBinLimits, TwoDpkArray, zSpaceBinnedPk, mean_losk, mean_perpk, zSpacemodesPerBin);
    
    output = fopen(filepath, "w");

    for(j=0; j<loskBinNumb-1; j++){
        for(k=0; k<perpkBinNumb-1; k++){
            fprintf(output, "%g \t %g \t %g \t %d \n", mean_losk[j][k], mean_perpk[j][k], TotalVolume*zSpaceBinnedPk[j][k], zSpacemodesPerBin[j][k]);
        }
    }

    fclose(output);

    return 0;
}


int polar2DpkBinning(char filepath[]){
    for(j=0; j<kBinNumb;  j++)        kBinLimits[j]  =               kbinInterval*j;
    for(j=0; j<muBinNumb; j++)       muBinLimits[j]  =   (1.0/(double) muBinNumb)*j;

    DualBinning(n0*n1*n2, muBinNumb, muBinLimits, kBinNumb, kBinLimits, polar2Dpk, polar2DBinnedPk, mean_mu, mean_modk, polar_modesPerBin);
    
    output = fopen(filepath, "w");

    for(j=0; j<muBinNumb-1; j++){
        for(k=0; k<kBinNumb-1; k++){
            fprintf(output, "%g \t %g \t %g \t %d \n", mean_mu[j][k], mean_modk[j][k], TotalVolume*polar2DBinnedPk[j][k], polar_modesPerBin[j][k]);
        }
    }

    fclose(output);
    
    return 0;
}


double legendrePolynomials(int Order, double mu){
    switch(Order){
        case 0: 
            return 1.;

        case 2:
            return 0.5*(3.*mu*mu -1.);
	       
        case 4:
            return (35.*pow(mu, 4.) - 30.*mu*mu + 3.)/8.;
    } 
}


int MultipoleCalc(int kBinNumb, double meankBin[], double kMonopole[], double kQuadrupole[], int ModesPerBin[], double** Array){
    // monopole     factor:        1. + (2./3.)*beta + 0.2*beta*beta;
    // quadrupole   factor:    (4./3.)*beta + (4./7.)*beta**2 = 0.68;
    // hexadecapole factor:   (8./35.)*beta*beta); 

    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    printf("\n\nPerforming multipole calculation.");
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = j*kbinInterval;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    printf("\nSorting mod k array.");
    
    qsort(Array, n0*n1*n2, sizeof(Array[0]), FirstColumnCompare);

    LowerBinIndex = 0;
    UpperBinIndex = 0;
    
    for(k=0; k<kBinNumb-1; k++){
        meankBin[k]    = 0.0;
        kMonopole[k]   = 0.0;
        kQuadrupole[k] = 0.0;
        ModesPerBin[k] =   0;
    }
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_%s_kbin_%.3f_00%d.dat", root_dir,  surveyType, kbinInterval, loopCount);
    else              sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_%s__kbin_%.3f_0%d.dat", root_dir, surveyType, kbinInterval, loopCount);

    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++){
        //  Linear regression against P_i(k, \mu_i)  
        //  L_i  = 0.5*(3.*\mu_i**2 -1.)
        //  P_i  = P_i(k, \mu_i)
    
        double Li          = 0.0;
        double Pi          = 0.0;
    
        double Sum_Li      = 0.0;
        double Sum_Li2     = 0.0;
    
        double Sum_Pi      = 0.0;
        double Sum_PiLi    = 0.0;
        
        // Find the range of indices corresponding to the modes in a given k interval. 
        for(i=LowerBinIndex; i<n0*n1*n2; i++){
            if(Array[i][0] > kBinLimits[j+1]){
                UpperBinIndex = i;
                break;
            } 
        }
        
        printf("\n %d \t %d", LowerBinIndex, UpperBinIndex);
        
        for(i=LowerBinIndex; i<UpperBinIndex; i++){
            meankBin[j]    += Array[i][0];
            modesPerBin[j] += 1;
            
                     // L_i = 0.5*(3.*mu**2 -1.)
                     
            Li              = 0.5*(3.*pow(Array[i][1], 2.) - 1.);     
            Pi              = Array[i][2];
            
            Sum_Li         += Li; 
            Sum_Li2        += Li*Li;
            
            Sum_Pi         += Pi;
            Sum_PiLi       += Pi*Li;
        }
    
        meankBin[j]        /= modesPerBin[j];
        	    
        // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
        // Required to invert AP  = B. 2x2 matrix inversion.
         	
        // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
        //            (c d)        

        // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li
     
        double detA;
        // det = ad - bc.
     
        detA = modesPerBin[j]*Sum_Li2 - Sum_Li*Sum_Li;
     
        // (P_0, P_2)^T = (A^-1)B = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
          kMonopole[j]  = (1./detA)*(Sum_Li2*Sum_Pi - Sum_Li*Sum_PiLi);
        kQuadrupole[j]  = (1./detA)*(-1.*Sum_Li*Sum_Pi + modesPerBin[j]*Sum_PiLi);
        
        fprintf(output, "%f \t %f \t %f \t %d \n", (float) meanKBin[j], (float) TotalVolume*kMonopole[j], (float) TotalVolume*kQuadrupole[j], modesPerBin[j]);   
    
        LowerBinIndex = UpperBinIndex;
    } 
    
    fclose(output);
    
    return 0;
}


float ShotNoise(){
    // Shot noise correction for a varying background number density. In the Peeble's convention for P(k).
    
    float numerator;
    float denominator;
    
    float flowChi;
    float fhiChi;
    
    flowChi = (float)  LowerChiLimit;
    fhiChi  = (float)  UpperChiLimit;
    
    numerator   =  qromb(&nbarChi2, flowChi, fhiChi);
    denominator =  qromb(&Chi2, flowChi, fhiChi);
    
    return numerator/denominator;
}


float nbarChi2(float Chi){
    return (float) pow((*pt2nz)(Chi), -1.)*Chi*Chi;
}


float Chi2(float Chi){
    return (float) Chi*Chi;
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
                
                kmodulus                           = pow(pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.), 0.5);
                
                mu                                 = k_x/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;      
                
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
				    
				    // The true density field is multiplied by a mask, correct for this with the measured P(k).
                    PkArray[Index][1]             *= TotalVolume*pow(fkpSqWeightsVolume, -1.);
                    
                    PkArray[Index][1]             -= (1./TotalVolume)*ShotNoise();
                 }
                
                // PkCorrections called to correct NGP for survey window function.
                if(WindowFuncParam == 1){
                    PkArray[Index][1]              = pow(H_kReal, 2.) + pow(H_kImag, 2.);
	            }
	            
	            polar2Dpk[Index][0]                = kmodulus;
	            polar2Dpk[Index][1]                = fabs(mu);
	            polar2Dpk[Index][2]                = PkArray[Index][1];
	            
	            TwoDpkArray[Index][0]              = fabs(k_x);                      // Line of sight wavevector. 
	            TwoDpkArray[Index][1]              = pow(k_y*k_y + k_z*k_z, 0.5);    // perpendicular wavevector.
                TwoDpkArray[Index][2]              = PkArray[Index][1];
	            
	        }
        }
    }
    
    printf("\nShot noise correction:  %e", ShotNoise());
    
    return 0;
}


int PkBinningCalc(int NumberModes, double** Array){
    printf("\n\nPerforming binning calculation.");
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = j*kbinInterval;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    printf("\nSorting mod k array.");
    
    qsort(Array, NumberModes, sizeof(Array[0]), FirstColumnCompare);

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
              
        // linearErrors[j]     = sqrt(2.*modesPerBin[j]/(TotalZADEWeight*TotalZADEWeight));

        del2[j]             = pow(meanKBin[j], 3.)*TotalVolume*binnedPk[j]*(4.*pi)/pow(2.*pi, 3.);
        
        LowerBinIndex       = UpperBinIndex;
        
        midKBin[j]          = 0.5*(kBinLimits[j] + kBinLimits[j+1]);
    }

    return 0;
}


int BinnedPkForMuInterval(double lowerMuLimit, double upperMuLimit, char filepath[], double Array[]){
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


int DualBinning(int NumberModes, int firstBinNumb, double firstBinLimits[], int secndBinNumb, double secondBinLimits[], double** DualParamArray, double** BinnedDualParamArray, double** mean_firstCol, double** mean_secndCol, int** modesPerBin){
    int m;
    
    printf("\nPerforming 2D binning calc.");
    
    qsort(DualParamArray, NumberModes, sizeof(DualParamArray[0]), FirstColumnCompare);

    // Order by first then second column.
    printf("\nDual param array sorted.");

    int firstColLowerBinIndex     = 0;
    int firstColUpperBinIndex     = 0;
    
    int secndColLowerBinIndex     = 0;
    int secndColUpperBinIndex     = 0;

    for(j=0; j<firstBinNumb-1; j++){
        for(k=0; k<secndBinNumb-1; k++){
            BinnedDualParamArray[j][k] = 0.0;
            mean_firstCol[j][k]        = 0.0;
            mean_secndCol[j][k]        = 0.0;
            modesPerBin[j][k]          =   0;    
        }
    }
    
    for(j=0; j<firstBinNumb-1; j++){
        for(i=firstColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][0] > firstBinLimits[j+1]){
                firstColUpperBinIndex = i;
                break;
            } 
        }

        secndColLowerBinIndex = firstColLowerBinIndex;

	qsort(&DualParamArray[firstColLowerBinIndex], firstColUpperBinIndex - firstColLowerBinIndex, sizeof(DualParamArray[0]), SecondColumnCompare);
       
        for(k=0; k<secndBinNumb-1; k++){
            for(m=secndColLowerBinIndex; m<firstColUpperBinIndex; m++){
                if(DualParamArray[m][1] > secondBinLimits[k+1]){
                    secndColUpperBinIndex = m;
                    break;
                } 
            }
        
            for(m=secndColLowerBinIndex; m<secndColUpperBinIndex; m++){
                BinnedDualParamArray[j][k]    += DualParamArray[m][2];
                mean_firstCol[j][k]           += DualParamArray[m][0];
                mean_secndCol[j][k]           += DualParamArray[m][1];
                modesPerBin[j][k]             += 1;    
            }
            
            if(modesPerBin[j][k] != 0)  BinnedDualParamArray[j][k] /= modesPerBin[j][k];
            if(modesPerBin[j][k] != 0)  mean_firstCol[j][k]        /= modesPerBin[j][k];
            if(modesPerBin[j][k] != 0)  mean_secndCol[j][k]        /= modesPerBin[j][k];
            
            secndColLowerBinIndex = secndColUpperBinIndex;
        }
        
        firstColLowerBinIndex = firstColUpperBinIndex;
    }

    return 0;
}
