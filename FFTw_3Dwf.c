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
    
    printf("\nFFT complete.");
    
    // printWindowfuncSlices();
    
    // char xaxis = 'x';
    // char yaxis = 'y';
    // char zaxis = 'z';
    
    // Currently not sensitive to FKP weights. 
    // WfSlice1D(5000, rand_x, 0, xaxis);
        
    // WfSlice1D(5000, rand_y, 1, yaxis);
    
    // WfSlice1D(5000, rand_z, 2, zaxis);

    // 0: Subtract shot noise for a real survey, 1: Neglect shot noise subtraction for FFT of window function. 
    PkCorrections(1);
    
    /*
    for(j=0; j<n0*n1*n2; j++){  
        k_x = kIntervalx*i;
        k_y = kIntervaly*j;
        k_z = kIntervalz*k;

        if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
        if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
        if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

        kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
        kmodulus                           = pow(pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.), 0.5);

        if(kmodulus < 0.25*NyquistWaveNumber){
            SumOverStates_W2k += PkArray[j][1];
        }
    }
    */
    
    printf("\nSum over states of W2k: %e", SumOverStates_W2k);

    PkBinningCalc(n0*n1*n2, PkArray, kbinInterval);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2q_%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%e \t %e\n", meanKBin[j], binnedPk[j]);

    // printf("\nWindow function P(k) calculation complete.");
    
    fclose(output);
    
    return 0;
}
