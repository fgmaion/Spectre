int wfPkCalc(){
    printf("\nBeginning window function calculation.");
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = (double) (TotalVolume/TotalSurveyedVolume)*booldensity[j]*FKPweights[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = (double) 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");

    PkCorrections(1);

    PkBinningCalc(n0*n1*n2);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2k_%s.dat", root_dir, surveyType);
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%f \t %f\n", meanKBin[j], binnedPk[j]);

    printf("\nWindow function P(k) calculation complete.");
    
    fclose(output);
    
    return 0;
}


int WindowfuncSlice(float kintervali, int ni, int x0, int y0, int z0, char filepath[]){
    // ints x0, y0, z0 specify whether it is true or false that the x,y and z co-ordinates have been taken to be zero for that slice.  
    
    float k_i                 = 0.0;
    int   Num_ModesInInterval = 0;

    for(j=0; j<ni; j++){

        k_i = kintervali*j;

        if(k_i>NyquistWaveNumber)  k_i    -= ni*kintervali;

        Index                              = ((float) z0)*n1*n2*j + ((float) y0)*n2*j + ((float) x0)*j;

        WindowFunc                         = 1.;
        
        if(k_i != 0.){
		  WindowFunc                      *= sin(pi*k_i*0.5/NyquistWaveNumber)/(pi*k_i*0.5/NyquistWaveNumber);}

        H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0]/WindowFunc;
        H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1]/WindowFunc;

        PkArray[Num_ModesInInterval][0]    = fabs(k_i);
        PkArray[Num_ModesInInterval][1]    = TotalVolume*(pow(H_kReal, 2.) + pow(H_kImag, 2.));
        Num_ModesInInterval               += 1;
    }
         
    printf("\n\nTotal number of modes for %d %d %d window fn. slice:  %d", x0, y0, z0, Num_ModesInInterval);          
         
    PkBinningCalc(Num_ModesInInterval);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)     fprintf(output, "%d \t %g \t %g\n", modesPerBin[j], meanKBin[j], binnedPk[j]);
    fclose(output);
    
    return 0;
}


int printWindowfuncSlices(){
    printf("\nPrinting window func. slices.");
    
    sprintf(WindowfuncSlices_dir, "%s/Data/WindowfuncSlices", root_dir);

    sprintf(Windowfunc_xSlice, "%s/%s_xSlice.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervalx, n2, 1, 0, 0, Windowfunc_xSlice);
    
    sprintf(Windowfunc_ySlice, "%s/%s_ySlice.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervaly, n1, 0, 1, 0, Windowfunc_ySlice);
    
    sprintf(Windowfunc_zSlice, "%s/%s_zSlice.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervalz, n0, 0, 0, 1, Windowfunc_zSlice);

    return 0;
}
