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
    
    //append a 'end of array' element to exit loop successfully. 
    rmodulus_vec[n0*n1*n2][0] =           9999.;
    rmodulus_vec[n0*n1*n2][1] = (double)   9999;
    
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
        // for(k=0; k<(kBinNumb-1); k++)   Corrfn[j] += unitTheory(rmodulus_vec[j], kBinLimits[k], 0.5*kbinInterval, 0)*(*pt2Pk)(kBinLimits[k]);
    
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
                H_kReal                            = pow(n0*n1*n2, -1.0)*in[Index][0];
                H_kImag                            = pow(n0*n1*n2, -1.0)*in[Index][1];
				    
                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		        polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		        polar2Dpk[polarPk_modeCount][2]    =  H_kReal;
		        
		        polarPk_modeCount                 += 1;
            }
        }
    }
            
    sprintf(filepath, "%s/Data/Multipoles/ConvolvedPk_Multipoles_%s_kbin_%.5f_kmax_%.2f_CellSize_%.2f.dat", root_dir,  surveyType, kbinInterval, modkMax, CellSize);
    
    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, modesPerBin, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 1);
    
    // sprintf(filepath, "%s/Data/Multipoles/GaussianWf_rMultipoles_%s_kbin_%.5f_00%d.dat", root_dir,  surveyType, kbinInterval, loopCount);
    
    // corrfn_multipoles(FFTW2_vecr_re, filepath);
    
    /*
    printf("\nSum over states of W2k: %e", SumOverStates_W2k);

    PkBinningCalc(n0*n1*n2, PkArray, kbinInterval);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2q_%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%e \t %e\n", meanKBin[j], binnedPk[j]);

    // printf("\nWindow function P(k) calculation complete.");
    
    fclose(output);
    */
    
    return 0;
}


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
		        polar2Dpk[polarPk_modeCount][2]    = W2_veck[Index]*pow(fkpWeightedVolume, -2.);
		        
		        polarPk_modeCount                 += 1;
            }
        }
    }
            
    sprintf(filepath, "%s/Data/Multipoles/VipersMaskWf_Multipoles_%s_kbin_%.5f.dat", root_dir,  surveyType, kbinInterval);
    
    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, modesPerBin, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 1);
    
    return 0;
}
