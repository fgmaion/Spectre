int PkCalc(){
    printf("\n\nAssigning FFT in array.");
    
    // The true density field multiplied by a mask, W(x). 
    for(j=0; j<n0*n1*n2; j++) in[j][0] = densityArray[j]*Cell_AppliedWindowFn[j]*BootStrap_Wght[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");
    
    PkCorrections();
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/Del2k/midK_Del2k_%s_kInterval_%.2f_00%d.dat", root_dir, surveyType, kbinInterval, loopCount);
    else              sprintf(filepath, "%s/Data/Del2k/midK_Del2k_%s_kInterval_%.2f_0%d.dat", root_dir, surveyType, kbinInterval, loopCount);
    
    // Monopole(filepath);

    observedQuadrupole(polarPk_modeCount);
        
    // Cartesian2Dpk();
    
    // sprintf(filepath, "%s/Data/muSlicedPk/%s_muSlicedPk_0.0mu0.2.dat", root_dir, surveyType);
    // BinnedPkForMuInterval(0.0, 0.2, filepath, PkCube);

    // sprintf(filepath, "%s/Data/muSlicedPk/%s_muSlicedPk_0.3mu0.5.dat", root_dir, surveyType);
    // BinnedPkForMuInterval(0.3, 0.5, filepath, clippedPk);

    // sprintf(filepath, "%s/Data/muSlicedPk/%s_muSlicedPk_0.6mu0.8.dat", root_dir, surveyType);
    // BinnedPkForMuInterval(0.6, 0.8, filepath, clippedPk);
    
    return 0;
}


int Monopole(char filepath[]){
    // Peacock 2dF:                 beta = 0.43 +- 0.07
    
    MonopoleCalc(kBinNumb, meanKBin, kMonopole, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0.0, 1.0, 1);;

    output = fopen(filepath, "w");
    
    // Only half the modes are independent. 
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %e \n", meanKBin[j], binnedPk[j], (*pt2Pk)(meanKBin[j]));

    fclose(output);

    return 0;
}


int observedQuadrupole(int modeCount){
    // sprintf(filepath, "%s/Data/Multipoles/Polar2Dpk_%s_%.3f_%d.dat", root_dir, surveyType, kbinInterval, loopCount);

    // polar2DpkBinning(modeCount);
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/Multipoles/Multipoles_%s_kbin_%.3f_00%d.dat", root_dir,  surveyType, kbinInterval, loopCount);
    else              sprintf(filepath, "%s/Data/Multipoles/Multipoles_%s_kbin_%.3f_0%d.dat", root_dir, surveyType, kbinInterval, loopCount);

    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, modeCount, filepath, kbinInterval, 0.0, 1.0, 1);
  
    // HexadecapoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, kHexadecapole, modesPerBin, polar2Dpk);

  return 0;
}


int Cartesian2Dpk(){
    for(j=0;  j<loskBinNumb;  j++)   loskBinLimits[j]  = (j+1)*(10.*kbinInterval);
    for(j=0; j<perpkBinNumb;  j++)   perpkBinLimits[j] = (j+1)*( 4.*kbinInterval);
    
    DualBinning(n0*n1*n2, loskBinNumb, loskBinLimits, perpkBinNumb, perpkBinLimits, TwoDpkArray, zSpaceBinnedPk, mean_losk, mean_perpk, zSpacemodesPerBin);
    
    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMask_2Dpk.dat", root_dir);
    
    output = fopen(filepath, "w");

    for(j=0; j<loskBinNumb-1; j++){
        for(k=0; k<perpkBinNumb-1; k++){
            fprintf(output, "%e \t %e \t %e \n", mean_losk[j][k], mean_perpk[j][k], TotalVolume*zSpaceBinnedPk[j][k]);
        }
    }

    fclose(output);
    
    return 0;
}


int polar2DpkBinning(int modeCount){
    for(j=0; j<kBinNumb;  j++)        kBinLimits[j]  =                     kbinInterval*(j+1);
    for(j=0; j<muBinNumb; j++)       muBinLimits[j]  =   (1.0/(double) (muBinNumb - 1))*(j+1);

    DualBinning(modeCount, muBinNumb, muBinLimits, kBinNumb, kBinLimits, polar2Dpk, polar2DBinnedPk, mean_mu, mean_modk, polar_modesPerBin);

    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMask_2Dpk_W2k.dat", root_dir);
    
    output = fopen(filepath, "w");

    printf("\nNumber of mu bins: %d", muBinNumb);

    for(j=0; j<muBinNumb-1; j++){
        for(k=0; k<kBinNumb-1; k++){
                    fprintf(output, "%e \t %e \t %e \t %e \n", mean_mu[j][k], mean_modk[j][k], TotalVolume*polar2DBinnedPk[j][k], (*pt2Pk)(mean_modk[j][k])*pow(1. + beta*pow(mean_mu[j][k], 2.), 2.)/(1. + 0.5*pow(mean_modk[j][k]*mean_mu[j][k]*velDispersion, 2.)));
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


int MultipoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double** Array, int modeCount, char filepath[], double Interval, double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    printf("\n\nPerforming multipole calculation. (to Quadrupole order)");

    double      midBin[modBinNumb];
    double   BinLimits[modBinNumb];
    int    modesperbin[modBinNumb];

    for(j=0; j<modBinNumb;     j++)     BinLimits[j]  =     (j+1)*Interval;
    
    for(j=0; j<(modBinNumb-1); j++)     midBin[j]     = (j + 1.5)*Interval;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    LowerBinIndex = 0;
    UpperBinIndex = 0;
    
    for(k=0; k<modBinNumb-1; k++){
        mean_modBin[k]    = 0.0;
        Monopole[k]       = 0.0;
        Quadrupole[k]     = 0.0;
        modesperbin[k]    =   0;
    }
    
    for(i=0; i<modeCount; i++){
        if(Array[i][0] >= BinLimits[0]){
            LowerBinIndex = i; 
            break;
        }
    }
    
    for(j=0; j<modBinNumb-1; j++){
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
        for(i=LowerBinIndex; i<modeCount; i++){
            if(Array[i][0] > BinLimits[j+1]){
                UpperBinIndex = i;
                break;
            } 
        }
        
        for(i=LowerBinIndex; i<UpperBinIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){        
                mean_modBin[j] += Array[i][0];
                modesperbin[j] += 1;
            
                         // L_i = 0.5*(3.*mu**2 -1.)
                     
                Li              = 0.5*(3.*pow(Array[i][1], 2.) - 1.);     
                Pi              = Array[i][2];
            
                Sum_Li         += Li; 
                Sum_Li2        += Li*Li;
            
                Sum_Pi         += Pi;
                Sum_PiLi       += Pi*Li;
            }
         }
    
         mean_modBin[j]        /= modesperbin[j];
        	    
         // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
         // Required to invert AP  = B. 2x2 matrix inversion.
         	
         // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
         //            (c d)        

         // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li
     
         double detA;
         // det = ad - bc.
     
         detA                       = modesperbin[j]*Sum_Li2 - Sum_Li*Sum_Li;
     
         // (P_0, P_2)^T = (A^-1)B  = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
                     Monopole[j]    = (1./detA)*( Sum_Li2*Sum_Pi - Sum_Li*Sum_PiLi);
                   Quadrupole[j]    = (1./detA)*(-Sum_Li*Sum_Pi  + modesperbin[j]*Sum_PiLi);
                   
         // for(jj=0; jj<= j; jj++)  printf("\n%d \t %e \t %e \t %e \t %d", j, mean_modBin[jj], Monopole[jj], Quadrupole[jj], modesperbin[jj]);
         
         // printf("\n\n");
         
         /*
         for(k=LowerBinIndex; k<UpperBinIndex; k++){
           if((mu_lolimit<Array[k][1]) && (Array[k][1]<mu_hilimit)){
             Li                       = 0.5*(3.*pow(Array[k][1], 2.) - 1.);
         
             kMonopole_expError[j]   += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(Sum_Li2 - Li*Sum_Li), 2.);
           
             kQuadrupole_expError[j] += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(-1.*Sum_Li + + modesperbin[j]*Li), 2.);
           }
         }

           kMonopole_expError[j]   *= pow(detA, -2.); 
         kQuadrupole_expError[j]   *= pow(detA, -2.);
         */
         
         LowerBinIndex   = UpperBinIndex;
     } 
    
     if(fileOutput == 1){       
         output = fopen(filepath, "w");
    
         for(j=0; j<modBinNumb-1; j++){
             fprintf(output, "%e \t %e \t %e \t %d \n", midBin[j], Monopole[j], Quadrupole[j], modesperbin[j]);   
         }
    
         fclose(output);
     }
    
     return 0;
}


int HexadecapoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double Hexadecapole[], double** Array, int modeCount, char filepath[], double Interval, double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.) + P_4(k)*(35x^4 -30x^2 +3)/8.
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    printf("\n\nPerforming Monopole calculation (to Hexadecapole order).");
    
    double midBin[modBinNumb];
    double BinLimits[modBinNumb];
    int    modesperbin[modBinNumb];
    
    for(j=0; j<modBinNumb;     j++)  BinLimits[j]  = (j + 1)*Interval;
    
    for(j=0; j<(modBinNumb-1); j++)  midBin[j]     = (j + 1.5)*Interval;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    LowerBinIndex = 0;
    UpperBinIndex = 0;
    
    for(k=0; k<modBinNumb-1; k++){
        mean_modBin[k]      = 0.0;
        Monopole[k]         = 0.0;
        Quadrupole[k]       = 0.0;
        Hexadecapole[k]     = 0.0;
        modesperbin[k]      =   0;
    }
    
    for(i=0; i<modeCount; i++){
        if(Array[i][0] >= BinLimits[0]){
            LowerBinIndex = i; 
            break;
        }
    }
    
    for(j=0; j<modBinNumb-1; j++){
        //  Linear regression against P_i(k, \mu_i)  
        //  L_i  = 0.5*(3.*\mu_i**2 -1.)
        //  P_i  = P_i(k, \mu_i)
        
        double Li          = 0.0;
        double Pi          = 0.0;
        double Ji          = 0.0;
    
        double Sum_Li      = 0.0;
        double Sum_Li2     = 0.0;
        double Sum_Ji      = 0.0;
        double Sum_Ji2     = 0.0;
    
        double Sum_Pi      = 0.0;
        double Sum_PiLi    = 0.0;
        double Sum_PiJi    = 0.0;
        double Sum_LiJi    = 0.0;
        
        // Find the range of indices corresponding to the modes in a given k interval. 
        for(i=LowerBinIndex; i<modeCount; i++){
            if(Array[i][0] > BinLimits[j+1]){
                UpperBinIndex = i;
                break;
            } 
        }
        
        for(i=LowerBinIndex; i<UpperBinIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){        
                mean_modBin[j]    += Array[i][0];
                modesperbin[j]    += 1;
                    
                Pi                 = Array[i][2];
                Li                 = 0.5*(3.*pow(Array[i][1], 2.) - 1.);    
                Ji                 = (35.*pow(Array[i][1], 4.) - 30.*pow(Array[i][1], 2.) + 3.)/8.; 
            
                Sum_Li            += Li; 
                Sum_Li2           += Li*Li;
                Sum_Ji            += Ji;
                Sum_Ji2           += Ji*Ji;
            
                Sum_Pi            += Pi;
                Sum_PiLi          += Pi*Li;
                Sum_PiJi          += Pi*Ji;
                Sum_LiJi          += Li*Ji;
            }
        }
    
        mean_modBin[j]            /= modesperbin[j];
        	    
        // For a matrix A, paramater vector, (P_0, P_2, P_4)^T, P and vector B
        // Required to invert AP  = B. 3x3 matrix inversion.
         	
        // A reads as (a b c) for a = sum_Modes 1, b = sum L_i, c = sum J_i, d = sum_Modes L_i, e = sum_Modes Li**2, f = sum_Modes Ji*Li, g = sum_Modes Ji, h = sum_Modes Li*Ji, i = sum_Modes Ji*Ji 
        //            (d e f)
        //            (g h i)        

        // and B = (b1, b2, b3)^T for b1 = sum_Modes hat Pi, b2 = sum_Modes (hat Pi)*Li, b3 = sum_Modes (hat Pi)*Ji
     
        double detA;
        // det = ad - bc.
     
        detA = modesperbin[j]*(Sum_Li2*Sum_Ji2 - Sum_LiJi*Sum_LiJi) - Sum_Li*(Sum_Li*Sum_Ji2 - Sum_LiJi*Sum_Ji) + Sum_Ji*(Sum_Li*Sum_LiJi - Sum_Li2*Sum_Ji);
     
        // (P _0, P_2)^T = (A^-1)B = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
            Monopole[j]  = (1./detA)*(     (Sum_Li2*Sum_Ji2 - Sum_LiJi*Sum_LiJi)*Sum_Pi - (Sum_Li*Sum_Ji2          - Sum_Ji*Sum_LiJi)*Sum_PiLi + (Sum_Li*Sum_LiJi         - Sum_Ji*Sum_Li2)*Sum_PiJi);
          Quadrupole[j]  = (1./detA)*( -1.*(Sum_Li*Sum_Ji2  - Sum_LiJi*Sum_Ji  )*Sum_Pi + (modesperbin[j]*Sum_Ji2  - Sum_Ji*Sum_Ji  )*Sum_PiLi - (modesperbin[j]*Sum_LiJi - Sum_Ji*Sum_Li )*Sum_PiJi);
        Hexadecapole[j]  = (1./detA)*(     (Sum_Li*Sum_LiJi - Sum_Li2*Sum_Ji   )*Sum_Pi - (modesperbin[j]*Sum_LiJi - Sum_Li*Sum_Ji  )*Sum_PiLi + (modesperbin[j]*Sum_Li2  - Sum_Li*Sum_Li )*Sum_PiJi);
    
        LowerBinIndex = UpperBinIndex;
    } 
    
    if(fileOutput == 1){
        output = fopen(filepath, "w");
    
        for(j=0; j<modBinNumb-1; j++){    
            fprintf(output, "%e \t %e \t %e \t %e \t %d \n", midBin[j], Monopole[j], Quadrupole[j], Hexadecapole[j], modesperbin[j]);   
        }
    
        fclose(output);
    }
    
    return 0;
}


double CubeShot(double irrev){
    // Comoving distance for evaluation is irrelevant. 
    return pow(CubeMeanNumberDensity(1500.), -1.);
}


double lightconeShot(double irrev){
    // Shot noise correction for a varying background number density. In the Peeble's convention for P(k).
    
    float numerator;
    float denominator;
    
    float flowChi;
    float fhiChi;
    
    flowChi = (float)  LowerChiLimit;
    fhiChi  = (float)  UpperChiLimit;
    
    numerator   =  qromb(&nbarChi2, flowChi, fhiChi);
    denominator =  qromb(&Chi2, flowChi, fhiChi);
    
    return (double) numerator/denominator;
}


float nbarChi2(float Chi){
    return (float) pow((*pt2nz)(Chi), -1.)*Chi*Chi;
}


float Chi2(float Chi){
    return (float) Chi*Chi;
}


int PkCorrections(){
    polarPk_modeCount = 0;

    double expectation, Power, amplitude;

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
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;      
                
                GaussianFilter                     = exp(-1.*kSq*0.5*(pow(1., 2.)));

                WindowFunc                         = 1.;

                if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);
                
                if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);
                
                if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);

                H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
                H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];
                
                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                // H_kReal                        /= pow(WindowFunc, 2.);
                // H_kImag                        /= pow(WindowFunc, 2.);
                
                PkArray[Index][0]                  = kmodulus;
        
                PkArray[Index][1]                  = pow(H_kReal, 2.) + pow(H_kImag, 2.);
        
                // Rescale measured P(k) in the cube to the V=1 camb convention.
                PkArray[Index][1]                 *= TotalVolume; 
                    
                // Account for the affect of the window on the measured P(k).
                PkArray[Index][1]                 /= fkpSqWeightsVolume*pow(TotalVolume, -1.);
                    
                // Clipping corrected shot noise estimate. 
                // PkArray[Index][1]              -= (1./TotalVolume)*(*pt2shot)(1.)*(1. - clippedVolume/TotalVolume);
                 
                // PkCorrections called to correct NGP for survey window function.

	            if((kmodulus > 0.000001) && (k_z >= 0.)){
	                // Only half the modes are independent.
	            	            
		            // Issue with mu for a zeroth length vector being ill defined. 
		            polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		            polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		            polar2Dpk[polarPk_modeCount][2]    = PkArray[Index][1];
		            
		            // printf("\n%e \t %e \t %e", kmodulus, fabs(mu), PkArray[Index][1]);
		            
		            // TwoDpkArray[Index][0]              = fabs(k_z);                      // Line of sight wavevector. 
	                // TwoDpkArray[Index][1]              = pow(k_y*k_y + k_x*k_x, 0.5);    // perpendicular wavevector.
                    // TwoDpkArray[Index][2]              = PkArray[Index][1];
                    
                    polarPk_modeCount                 += 1;
	            }
	        }
        }
    }
    
    // sprintf(filepath, "%s/Data/SpectralDistortion/fullCube_muBinsPk_0.0.dat", root_dir);
    
    // BinnedPkForMuInterval(0.0, 0.2, filepath, n0*n1*n2);
    
    return 0;
}


int Gaussianfield(){
    int m0, m1, m2;
    
    double Power, amplitude, phase, expectation; 
    
    printf("\nTotal vol.: %e", TotalVolume);
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                m0 = k;
                m1 = j;
                m2 = i;

                if(m2>n2/2)  m2                   -= n2;
                if(m1>n1/2)  m1                   -= n1;
                if(m0>n0/2)  m0                   -= n0;
                
                k_x                                = kIntervalx*m2;
                k_y                                = kIntervaly*m1;
                k_z                                = kIntervalz*m0;
                
                Index                              = k*n1*n2 + j*n2 + i;
                
                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);

                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;
                                                                                     // shot noise contribution. 
                // expectation                     = (*pt2Pk)(kmodulus)/TotalVolume; //  + (1./TotalVolume)*(*pt2shot)(1.);

                // Monopole and Quadrupole components, with z taken as the line of sight direction. 
                // expectation                     = (kaiser_multipole(kmodulus, beta, 0) + kaiser_multipole(kmodulus, beta, 2)*0.5*(3.*mu*mu -1.))*Pk_powerlaw(kmodulus, 5., 1.8)/TotalVolume;

                expectation                        = (*pt2Pk)(kmodulus)/TotalVolume;
                
                expectation                       *= 1. + 0.5*pow(mu, 2.);
            
                // expectation                    *= pow(1. + beta*pow(mu, 2.), 2.);
                
                // expectation                    /= 1. + 0.5*pow(kmodulus*mu*velDispersion, 2.);

                Power                              = -log(gsl_rng_uniform(gsl_ran_r))*expectation;
                
                amplitude                          = sqrt(Power);
                // amplitude                          = sqrt(expectation);
                
                phase                              = 2.*pi*gsl_rng_uniform(gsl_ran_r);
                
                // purely real.
                if(k_x == NyquistWaveNumber) phase = 0.0;
                if(k_y == NyquistWaveNumber) phase = 0.0;
                if(k_z == NyquistWaveNumber) phase = 0.0;
                
                // printf("\n%le %le \t %le \t %le \t %le", kmodulus, (*pt2Pk)(kmodulus), Power, amplitude, phase);
                
                // Assuming Cic assignment scheme
                out[Index][0]                      = amplitude*cos(phase);
                out[Index][1]                      = amplitude*sin(phase);
                
                WindowFunc                         = 1.;

                if(k_x != 0.){
		            WindowFunc                    *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);}
                
                if(k_y != 0.){
		            WindowFunc                    *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);}
                
                if(k_z != 0.){
		            WindowFunc                    *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);}		      
	        
	            // out[Index][0]                  *= pow(WindowFunc, 2.);
	        	// out[Index][1]                  *= pow(WindowFunc, 2.);
	        }
        }
    }
    
    // Zero mean. 
    out[0][0] = 0.0;
    out[0][1] = 0.0;
    
    int negkIndex;
            
    // Hermitian condition. 
    for(k=0; k<n0/2 +1; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index     = k*n1*n2 + j*n2 + i;
                       
                negkIndex = 0;
                     
                if(i != 0) negkIndex += (n2-i); 
                if(j != 0) negkIndex += (n1-j)*n2; 
                if(k != 0) negkIndex += (n0-k)*n1*n2; 

                out[negkIndex][0] =     out[Index][0];
                out[negkIndex][1] = -1.*out[Index][1];
            }
        }
    }
    
    fftw_execute(iplan);

    for(j=0; j<n0*n2*n1; j++)  densityArray[j] = in[j][0];
    
    double GRF_var = 0.0, GRF_mean = 0.0;
    
    apparent_mean = 0.0;
    
    for(j=0; j<n0*n1*n2; j++){ 
        GRF_var        += pow(densityArray[j], 2.);
        GRF_mean       +=     densityArray[j];
    }
    
    GRF_var     /= n0*n1*n2;
    GRF_mean    /= n0*n1*n2;                           
    
    printf("\n\nMean of the Gaussian realisation: %e", GRF_mean);
    printf("\n\nVariance of Gaussian realisation: %e", GRF_var);
    
    for(j=0; j<n0*n2*n1; j++)  apparent_mean       += densityArray[j]*Cell_AppliedWindowFn[j];
    
    // Number of non-empty cells in mask. 
    apparent_mean    /= 123004.;
    
    // for(j=0; j<n0*n1*n2; j++)  densityArray[j]     -= apparent_mean;
    
    printf("\n\nApparent (weighted) mean of the Gaussian realisation: %e", apparent_mean);
    
    return 0;
}


int BinnedPkForMuInterval(double lowerMuLimit, double upperMuLimit, char filepath[], int modeCount){
    printf("\nPerforming P(k) binning for given mu interval, %f < mu < %f", lowerMuLimit, upperMuLimit);

    int muInterval_modeCount = 0; 

    for(k=0; k<modeCount; k++){
        if((lowerMuLimit<polar2Dpk[k][1]) && (polar2Dpk[k][1] < upperMuLimit)){
            muIntervalPk[k][0]    = polar2Dpk[k][0];
            muIntervalPk[k][1]    = polar2Dpk[k][2];  
        
            muInterval_modeCount += 1;
        }
    }
    
    printf("\nNumber of modes in interval:  %d", muInterval_modeCount);
    
    MonopoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0.0, 1.0, 1);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++)     fprintf(output, "%e \t %e \t %e \n", meanKBin[j], TotalVolume*binnedPk[j], splintConvMono(midKBin[j]) + splintConvQuad(midKBin[j])*LegendrePolynomials(0.5*(upperMuLimit - lowerMuLimit), 2));
    
    fclose(output);
    
    return 0;
}


int MonopoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double** Array, int modeCount, char filepath[], double Interval, double mu_lolimit, double mu_hilimit, int fileOutput){
    printf("\n\nPerforming multipole calculation. (to Monopole order)");
    
    double      midBin[modBinNumb];
    double   BinLimits[modBinNumb];
    int    modesperbin[modBinNumb];

    for(j=0; j<modBinNumb;     j++)     BinLimits[j]  =     (j+1)*Interval;
    
    for(j=0; j<(modBinNumb-1); j++)     midBin[j]     = (j + 1.5)*Interval;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    LowerBinIndex = 0;
    UpperBinIndex = 0;
    
    for(j=0; j<modBinNumb-1; j++){
        modesperbin[j]  =              0;
        mean_modBin[j]  =            0.0;
        Monopole[j]     =            0.0;
        // linearErrors[j] =            0.0;
    }
    
    for(j=0; j<modeCount; j++){
        if(Array[j][0]   >= BinLimits[0]){
            LowerBinIndex = j; 
            break;
        }
    }
    
    for(j=0; j<modBinNumb-1; j++){
        for(i=LowerBinIndex; i<modeCount; i++){
            if(Array[i][0] > BinLimits[j+1]){
                UpperBinIndex = i;
                break;
            } 
        }
        
        for(i=LowerBinIndex; i<UpperBinIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){      
                // printf("\n%d \t %e \t %e \t %e", j, Array[i][0], Array[i][1], Array[i][2]);
              
                mean_modBin[j] += Array[i][0];
                   Monopole[j] += Array[i][2];
                modesperbin[j] += 1;
            }
        }
        	    
        if(modesperbin[j]  != 0)  mean_modBin[j]  /= modesperbin[j];
        if(modesperbin[j]  != 0)     Monopole[j]  /= modesperbin[j];

        // Peacock and Nicholson 1991, pg 313. above eqn (20).
        // Result of summing over a shell in k space containing m modes, should be a Gaussian random variable with variance 2.m/N^2  
              
        // linearErrors[j]     = sqrt(2.*modesPerBin[j]/(TotalZADEWeight*TotalZADEWeight));

        // del2[j]             = pow(mean_modBin[j], 3.)*TotalVolume*binnedPk[j]*(4.*pi)/pow(2.*pi, 3.);
        
        LowerBinIndex       = UpperBinIndex;
    }

    if(fileOutput == 1){
        output = fopen(filepath, "w");
        
        fprintf(output, "%e \t %e \t %d \n", Array[0][0], Array[0][2], 1);
    
        for(j=0; j<modBinNumb-1; j++){    
            fprintf(output, "%e \t %e \t %d \n", mean_modBin[j], Monopole[j], modesperbin[j]);   
        }
    
        fclose(output);
    }

    return 0;
}


int DualBinning(int NumberModes, int firstBinNumb, double firstBinLimits[], int secndBinNumb, double secondBinLimits[], double** DualParamArray, double** BinnedDualParamArray, double** mean_firstCol, double** mean_secndCol, int** modesPerBin){
    int m;
    
    printf("\nPerforming 2D binning calc.");
    
    qsort(DualParamArray, NumberModes, sizeof(DualParamArray[0]), FirstColumnCompare);

    // for(j=0; j<NumberModes; j++) printf("\n%e \t %e \t %e", DualParamArray[j][0], DualParamArray[j][1], DualParamArray[j][2]);  

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
    
    for(j=0; j<NumberModes; j++){
        if(DualParamArray[j][0]   >= firstBinLimits[0]){
            firstColLowerBinIndex = j; 
        
            break;
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
    
        // meets the lower bin limit
        for(i=secndColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][1]   >= secondBinLimits[0]){
                secndColLowerBinIndex = i; 
        
                break;
            }
        }      
       
        // meets the upper bin limit. 
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
        
        // printf("\n %d \t %d", firstColLowerBinIndex, firstColUpperBinIndex);
        
        firstColLowerBinIndex = firstColUpperBinIndex;
    }

    return 0;
}
