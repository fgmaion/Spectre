int WfSlice1D(int n0, double iCoors[], int axislabel, char axis){
    prepFFTw1D(n0);    

    double   tinyCellSize;

    tinyCellSize   = 3000./n0; 

    for(j=0; j<n0; j++){
        in1D[j][0] = 0.0;
        in1D[j][1] = 0.0;
    }
    
    for(j=0; j<rand_number; j++){
                                                      // lower limit for the i'th axis. 
        boxlabel           = (int) floor((iCoors[j] - AxisLimsArray[0][axislabel])/tinyCellSize);        
    
        in1D[boxlabel][0] += 1;
    }   
    
    fftw_execute(p1D);
    
    sprintf(Windowfunc_xSlice, "%s/%s_1D_%d_%cSlice_NGPuncorr.dat", WindowfuncSlices_dir, surveyType, n0, axis);
    
    output = fopen(Windowfunc_xSlice, "w");
    
    for(j=0; j<n0/2 + 1; j++)  fprintf(output, "%.9e \t %.9e\n", log10(j*2.*pi/3000.), log10(pow(out1D[j][0]/out1D[0][0], 2.) + pow(out1D[j][1]/out1D[0][0], 2.)));

    fclose(output);

    return 0;
}


int polar2dWf(){
    int polarPk_modeCount = 0;

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

	            if(kmodulus > 0.000001){	            
		            // Issue with mu for a zeroth length vector being ill defined. 
		            polar2Dpk[polarPk_modeCount][0]    = fabs(mu);
		            polar2Dpk[polarPk_modeCount][1]    = kmodulus;
		            polar2Dpk[polarPk_modeCount][2]    = W2_veck[Index];
		            
    	            polarPk_modeCount                 += 1;
                }		            
            }
        }
    }

    polar2DpkBinning(polarPk_modeCount);

    return 0;
}


int wfPkCalc(){    
    printf("\n\nBeginning window function calculation.");
    
    // The true density field is multiplied by a mask, Cell_AppliedWindowFn [W(x)].
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nWindow function FFT complete.");
    
    for(j=0; j<n0*n1*n2; j++) W2_veck[j] = pow(n0*n1*n2, -2.0)*(pow(out[j][0], 2.) + pow(out[j][1], 2.));
    
    // polar2dWf();
    
    // W2k_Multipoles();
    
    // 3D FFT of W^2(vec k).
    for(j=0; j<n0*n1*n2; j++)  out[j][0]  = W2_veck[j];
    for(j=0; j<n0*n1*n2; j++)  out[j][1]  = 0.0;
    
    fftw_execute(iplan);
    
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_re[j] =  in[j][0];
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_im[j] =  in[j][1];
    
    for(j=0; j<n0*n1*n2; j++) W2_vecr[j] = pow(FFTW2_vecr_re[j], 2.) + pow(FFTW2_vecr_im[j], 2.);
    
    /*
    int    m0, m1, m2;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               m0 = k;
               m1 = j;
               m2 = i;

               if(m2>n2/2)  m2 -= n2;
               if(m1>n1/2)  m1 -= n1;
               if(m0>n0/2)  m0 -= n0;

               Index         = k*n1*n2 + j*n2 + i;
               
               rmodulus_vec[Index][0] =           CellSize*sqrt(m2*m2 + m1*m1 + m0*m0);
               rmodulus_vec[Index][1] = (double)                                 Index;
            }
        }
    }
    
    //append an 'end of array' element to exit loop successfully. 
    rmodulus_vec[n0*n1*n2][0] =  9999.;
    rmodulus_vec[n0*n1*n2][1] =  9999.;
    
    qsort(rmodulus_vec, n0*n1*n2, sizeof(rmodulus_vec[0]), FirstColumnCompare);
    
    printf("\n\nBeginning corr. fn. reassignment.");
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = (j+1)*kbinInterval;
    
    j=0; 
    
    // j is the position in the sorted array. 
    while(j<n0*n1*n2){
        Interim = 0.0;
    
        for(k=0; k<(kBinNumb-1); k++)  Interim  += unitTheory(rmodulus_vec[j][0], kBinLimits[k], 0.5*kbinInterval, 0)*TotalVolume*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, 0);
    
        for(k=j; k<=n0*n1*n2; k++){
            i                                    = (int) ceil(rmodulus_vec[k][1]);
              
            Corrfn[i]                            = Interim;
        
            if(rmodulus_vec[k][0] > rmodulus_vec[j][0]){
                j = k;
                
                break;
            }
        }
    }
    
    Corrfn[0]   = 0.0;
    
    printf("\n\nCorrelaton fn. reassigned in the stepwise P(k) assumption.");
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
                
                Index                              = k*n1*n2 + j*n2 + i;
                
                kmodulus_vec[Index]                = kmodulus;
                mu_vec[Index]                      = fabs(mu);
            }
        }
    }
    
    // Carry out the convolution by inverse FFT of the correlation fn. x FFT of W^2(k).
    for(j=0; j<n0*n1*n2; j++){ 
        // Cannot use kMonopole values due to aliased measurement. 
        in[j][0]        =  Corrfn[j]*FFTW2_vecr_re[j];
        in[j][1]        =  Corrfn[j]*FFTW2_vecr_im[j];
    }
    
    fftw_execute(p);
        
    printf("\nWindow fn. convolution complete.");

    free2DBinning();

    assign2DPkMemory(muBinNumb, kBinNumb);
    
    // W2k_Multipoles();
    
    polarPk_modeCount = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
                
                Index                              = k*n1*n2 + j*n2 + i;
                
                kmodulus_vec[Index]                = kmodulus;
                mu_vec[Index]                      = fabs(mu);

                // Binning calculation following inverse transform. 
                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		        polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		        
		        // Convolved P(k), amplitude corrected for the effects of the convolution. 
		        polar2Dpk[polarPk_modeCount][2]    = pow(n0*n1*n2, -1.)*out[Index][0];
		        
		        polarPk_modeCount                 += 1;
            }
        }
    }
    
    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0);
    
    sprintf(filepath, "%s/Data/Multipoles/ConvolvedPk_Multipoles_kbin_%.5f_kmax_%.2f_CellSize_%.2f_MonoInput.dat", root_dir,  surveyType, kbinInterval, modkMax, CellSize);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<(kBinNumb-1); j++) fprintf(output, "%e \t %e \t %e \n", meanKBin[j], kMonopole[j], kQuadrupole[j]);
    
    fclose(output);
    
    // sprintf(filepath, "%s/Data/Multipoles/GaussianWf_rMultipoles_%s_kbin_%.5f_00%d.dat", root_dir,  surveyType, kbinInterval, loopCount);
    
    // corrfn_multipoles(FFTW2_vecr_re, filepath);
    */
    
    // W2r_Multipoles();
    
    newApproach_W2rMultipoles();
    
    return 0;
}


int newApproach_W2rMultipoles(){
    int    rbinNumb   =  200;
    
    double rbinlength =  3.0;
    
    windowfn_rspaceMultipoles(rbinNumb);
    
    for(j=0; j<rbinNumb;  j++){  
        windfn_rvals[j]         =         j*rbinlength;
    
        windfn_rMonopole[j]     =                  0.0;
        
        windfn_rQuadrupole[j]   =                  0.0;
                
        windfn_rHexadecapole[j] =                  0.0;
    }
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;
    
                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                    
                kmodulus                           = pow(kSq, 0.5);
                    
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
                
                Index                              = k*n1*n2 + j*n2 + i;
                    
                mu_vec[Index]                      = mu;
                    
                kmodulus_vec[Index]                = kmodulus;
            }
        }
    }
    
    printf("\n\nCalculating window autocorrelation monopole");
    
    for(j=0; j<rbinNumb; j++){
        printf("\n %d", j);
    
        for(k=0; k<n0*n1*n2; k++){
            windfn_rMonopole[j]                   +=  W2_veck[k]*gsl_sf_bessel_j0(kmodulus_vec[k]*windfn_rvals[j]);
            
            windfn_rQuadrupole[j]                 +=  W2_veck[k]*gsl_sf_bessel_j2(kmodulus_vec[k]*windfn_rvals[j])*LegendrePolynomials(mu_vec[k], 2);
            
            windfn_rHexadecapole[j]               +=  W2_veck[k]*gsl_sf_bessel_jl(4, kmodulus_vec[k]*windfn_rvals[j])*LegendrePolynomials(mu_vec[k], 4);
        }
    }
    
    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMaskWf_rSpaceMonopole_newApproach.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<rbinNumb; j++)  fprintf(output, "%e \t %e \t %e \t %e \n", windfn_rvals[j], windfn_rMonopole[j], windfn_rQuadrupole[j], windfn_rHexadecapole[j]);
    
    fclose(output);

    return 0;
}


/*
int wfPkCalc(){
    // Delivers the spherical average of W^2(q).
    
    printf("\n\nBeginning window function calculation.");
    
    // The true density field is multiplied by a mask, Cell_AppliedWindowFn [W(x)].
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nWindow function FFT complete.");
    
    for(j=0; j<n0*n1*n2; j++) W2_veck[j] = pow(n0*n1*n2, -2.0)*(pow(out[j][0], 2.) + pow(out[j][1], 2.));
    
    // 3D FFT of W^2(vec k).
    for(j=0; j<n0*n1*n2; j++)  in[j][0]  = W2_veck[j];
    for(j=0; j<n0*n1*n2; j++)  in[j][1]  = 0.0;
    
    fftw_execute(p);
    
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_re[j] =           pow(n0*n1*n2, -1.0)*out[j][0];
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_im[j] =           pow(n0*n1*n2, -1.0)*out[j][1];
    
    int    m0, m1, m2;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               m0 = k;
               m1 = j;
               m2 = i;

               if(m2>n2/2)  m2 -= n2;
               if(m1>n1/2)  m1 -= n1;
               if(m0>n0/2)  m0 -= n0;

               Index         = k*n1*n2 + j*n2 + i;
               
               rmodulus_vec[Index][0] =           CellSize*sqrt(m2*m2 + m1*m1 + m0*m0);
               rmodulus_vec[Index][1] = (double)                                 Index;
            }
        }
    }
    
    //append an 'end of array' element to exit loop successfully. 
    rmodulus_vec[n0*n1*n2][0] =  9999.;
    rmodulus_vec[n0*n1*n2][1] =  9999.;
    
    qsort(rmodulus_vec, n0*n1*n2, sizeof(rmodulus_vec[0]), FirstColumnCompare);
    
    printf("\n\nBeginning corr. fn. reassignment.");
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = (j+1)*kbinInterval;
    
    j=0; 
    
    // j is the position in the sorted array. 
    while(j<n0*n1*n2){
        Interim = 0.0;
    
        for(k=0; k<(kBinNumb-1); k++)  Interim  += unitTheory(rmodulus_vec[j][0], kBinLimits[k], 0.5*kbinInterval, 0)*(*pt2Pk)(kBinLimits[k]);
    
        for(k=j; k<=n0*n1*n2; k++){
            i                                    = (int) ceil(rmodulus_vec[k][1]);
              
            Corrfn[i]                            = Interim;
        
            if(rmodulus_vec[k][0] > rmodulus_vec[j][0]){
                j = k;
                
                break;
            }
        }
    }
    
    Corrfn[0]   = 0.0;
    
    printf("\n\nCorrelaton fn. reassigned in the stepwise P(k) assumption.");
    
    // Carry out the convolution by inverse FFT of the correlation fn. x FFT of W^2(k).
    for(j=0; j<n0*n1*n2; j++){ 
        // Cannot use kMonopole values due to aliased measurement. 
        out[j][0]        =              Corrfn[j]*FFTW2_vecr_re[j];
        out[j][1]        =              Corrfn[j]*FFTW2_vecr_im[j];
    }
    
    fftw_execute(iplan);
        
    printf("\nWindow fn. convolution complete.");
    
    // printWindowfuncSlices();
    
    // char xaxis = 'x';
    // char yaxis = 'y';
    // char zaxis = 'z';
    
    // Currently not sensitive to FKP weights. 
    // WfSlice1D(5000, rand_x, 0, xaxis);
        
    // WfSlice1D(5000, rand_y, 1, yaxis);
    
    // WfSlice1D(5000, rand_z, 2, zaxis);

    // 0: Subtract shot noise for a real survey, 1: Neglect shot noise subtraction for FFT of window function. 
    // PkCorrections(1);
    
    free2DBinning();

    assign2DPkMemory(muBinNumb, kBinNumb);
    
    W2k_Multipoles();
    
    polarPk_modeCount = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
                
                Index                              = k*n1*n2 + j*n2 + i;
                
                kmodulus_vec[Index]                = kmodulus;
                mu_vec[Index]                      = fabs(mu);

                // Binning calculation following inverse transform. 
                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		        polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		        
		        // Convolved P(k), amplitude corrected for the effects of the convolution. 
		        polar2Dpk[polarPk_modeCount][2]    = TotalVolume*in[Index][0]/fkpSqWeightsVolume;
		        
		        polarPk_modeCount                 += 1;
            }
        }
    }
    
    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0);
    
    sprintf(filepath, "%s/Data/Multipoles/ConvolvedPk_Multipoles_%s_kbin_%.5f_kmax_%.2f_CellSize_%.2f.dat", root_dir,  surveyType, kbinInterval, modkMax, CellSize);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<(kBinNumb-1); j++) fprintf(output, "%e \t %e \t %e \n", midKBin[j], kMonopole[j], (*pt2Pk)(midKBin[j]));
    
    fclose(output);
    
    // sprintf(filepath, "%s/Data/Multipoles/GaussianWf_rMultipoles_%s_kbin_%.5f_00%d.dat", root_dir,  surveyType, kbinInterval, loopCount);
    
    // corrfn_multipoles(FFTW2_vecr_re, filepath);
    
    return 0;
}
*/

int W2k_Multipoles(){
    polarPk_modeCount = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
                
                Index                              = k*n1*n2 + j*n2 + i;

                // Binning calculation following inverse transform. 
                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		        polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		        polar2Dpk[polarPk_modeCount][2]    = W2_veck[Index];
		      
		        polarPk_modeCount                 += 1;
            }
        }
    }
            
    sprintf(filepath, "%s/Data/Multipoles/VipersMaskWf_Multipoles_%s_kbin_%.5f.dat", root_dir,  surveyType, kbinInterval);
    
    MultipoleCalc(kBinNumb, meanKBin, wfMonopole, wfQuadrupole, polar2Dpk, polarPk_modeCount, filepath, 0.0, 1.0, kbinInterval, 0);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<(kBinNumb-1); j++) fprintf(output, "%e \t %e \t %e \n", meanKBin[j], wfMonopole[j], wfQuadrupole[j]);
    
    fclose(output);
    
    return 0;
}


int W2r_Multipoles(){
    int    rbinNumb   =  120;
    
    double rbinlength =  4.0;
    
    int    m0, m1, m2;

    double rbinLimits[rbinNumb];
    
    double xCell, yCell, zCell, rmodulus;
    
    free2DBinning();

    assign2DPkMemory(muBinNumb, rbinNumb);
    
    for(j=0; j<rbinNumb;  j++)    rbinLimits[j]    =    j*rbinlength;
    for(j=0; j<muBinNumb; j++)  muBinLimits[j]     =   (1.0/(double) (muBinNumb - 1))*(j+1);
    
    polarPk_modeCount        = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               m0 = k;
               m1 = j;
               m2 = i;

               if(m2>n2/2)  m2 -= n2;
               if(m1>n1/2)  m1 -= n1;
               if(m0>n0/2)  m0 -= n0;
            
	           xCell         =  CellSize*m2;
	           yCell         =  CellSize*m1;
	           zCell         =  CellSize*m0;

               rmodulus      = CellSize*sqrt(m2*m2 + m1*m1 + m0*m0);

               Index         = k*n1*n2 + j*n2 + i;
               
               rmodulus_vec[Index][0]             =        rmodulus;
               rmodulus_vec[Index][1]             = (double)  Index;

               mu                                 = m0/sqrt(m2*m2 + m1*m1 + m0*m0);
               if(rmodulus < 0.000001)       mu   = 0.0;      

		       polar2Dpk[polarPk_modeCount][0]    = rmodulus;
               polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		       polar2Dpk[polarPk_modeCount][2]    = FFTW2_vecr_re[Index];
                    
               polarPk_modeCount                 += 1;
            }
        }
    }
    
    printf("\n\nzero point of window fn. in delta space: %e", W2_vecr[0]);
    
    rmodulus_vec[n0*n1*n2][0] =           9999.;
    rmodulus_vec[n0*n1*n2][1] = (double)   9999;
    
    // DualBinning(polarPk_modeCount, muBinNumb, muBinLimits, rbinNumb, rbinLimits, polar2Dpk, polar2DBinnedPk, mean_mu, mean_modk, polar_modesPerBin);

    windowfn_rspaceMultipoles(rbinNumb);

    // sprintf(filepath, "%s/Data/SpectralDistortion/VipersBinaryMask_rSpaceMultipoles_Hex.dat", root_dir);

    // MultipoleCalc(rbinNumb, windfn_rvals, windfn_rMonopole, windfn_rQuadrupole, polar2Dpk, polarPk_modeCount, filepath, rbinlength, 0.0, 1.0, 1);
    
    // HexadecapoleCalc(rbinNumb, windfn_rvals, windfn_rMonopole, windfn_rQuadrupole, windfn_rHexadecapole, polar2Dpk, polarPk_modeCount, filepath, rbinlength, 0.0, 1.0, 1);
       
    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMaskWf_rSpaceMonopole.dat", root_dir);
    
    MonopoleCalc(rbinNumb, windfn_rvals, windfn_rMonopole, polar2Dpk, polarPk_modeCount, filepath, rbinlength, 0.0, 1.0, 1);
       
    return 0;
}


double splint_windfn_rSpaceMono(double r){
    double Interim;
    
    splint(windfn_rvals, windfn_rMonopole, windfn_rMonopole2d, rbinNumb, r, &Interim);
    
    // if(r>windfn_rvals[rbinNumb-1])  Interim = 0.0;
    
    if(r>windfn_rvals[rbinNumb-1])  Interim = 0.0;
    
    return Interim;
}


double splint_windfn_rSpaceQuad(double r){
    double Interim;
    
    splint(windfn_rvals, windfn_rQuadrupole, windfn_rQuadrupole2d, rbinNumb, r, &Interim);
    
    if(r>windfn_rvals[rbinNumb-1])  Interim = 0.0;
    
    return Interim;
}


double splint_windfn_rSpaceHex(double r){
    double Interim;
    
    splint(windfn_rvals, windfn_rHexadecapole, windfn_rHexadecapole2d, rbinNumb, r, &Interim);
    
    if(r>windfn_rvals[rbinNumb-1])  Interim = 0.0;
    
    return Interim;
}


int spline_wfMultipoles_deltaSpace(){
    rbinNumb   =  200;

    windowfn_rspaceMultipoles(rbinNumb);
    
    //  sprintf(filepath, "%s/Data/SpectralDistortion/VipersBinaryMask_rSpaceMultipoles_Hex.dat", root_dir);

    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMaskWf_rSpaceMonopole_newApproach.dat", root_dir);

    inputfile = fopen(filepath, "r");
    
    // for(j=0; j<rbinNumb; j++)  fscanf(inputfile, "%le \t %le \t %le \t \t %le %*d \n", &windfn_rvals[j], &windfn_rMonopole[j], &windfn_rQuadrupole[j], &windfn_rHexadecapole[j]);
    
    for(j=0; j<rbinNumb; j++){  
        fscanf(inputfile, "%le \t %le \t %le \t %le \n", &windfn_rvals[j], &windfn_rMonopole[j], &windfn_rQuadrupole[j], &windfn_rHexadecapole[j]);
    
        printf("\n%e \t %e \t %e \t %e", windfn_rvals[j], windfn_rMonopole[j], windfn_rQuadrupole[j], windfn_rHexadecapole[j]);
    }
    
    fclose(inputfile);
    
    spline(windfn_rvals, windfn_rMonopole,     rbinNumb, 1.0e31, 1.0e31, windfn_rMonopole2d);
        
    spline(windfn_rvals, windfn_rQuadrupole,   rbinNumb, 1.0e31, 1.0e31, windfn_rQuadrupole2d);
    
    spline(windfn_rvals, windfn_rHexadecapole, rbinNumb, 1.0e31, 1.0e31, windfn_rHexadecapole2d);
    
    printf("\n\nlimit on splint taken to be: %e", windfn_rvals[rbinNumb-1]);
    
    sprintf(filepath, "%s/Data/SpectralDistortion/VipersBinaryMask_rSpaceMultipoles_MonoSplint.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double delta;
    
    for(j=0; j<6000; j++){    
        delta = j/10.;
    
        fprintf(output, "%e \t %e \t %e \n", delta, splint_windfn_rSpaceMono(delta), splint_windfn_rSpaceQuad(delta));
    }
    
    fclose(output);
    
    return 0;
}
