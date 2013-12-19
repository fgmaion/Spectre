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


int WindowfuncSlice(double kintervali, int ni, int x0, int y0, int z0, char filepath[]){
    double k_i                 = 0.0;
    int   Num_ModesInInterval  = 0;

    for(j=0; j<ni; j++){
        k_i = kintervali*j;
        
        if(k_i>NyquistWaveNumber)  k_i    -= ni*kintervali;

        Index                              = n1*n2*z0*j + n2*y0*j + x0*j;

        // Cell window fn. uncorrected. 
        H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
        H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];

        PkArray[Num_ModesInInterval][0]    = fabs(k_i);
        PkArray[Num_ModesInInterval][1]    = pow(H_kReal, 2.) + pow(H_kImag, 2.);
        Num_ModesInInterval               += 1;
    }
         
    printf("\n\nTotal number of modes for %d %d %d window fn. slice:  %d", x0, y0, z0, Num_ModesInInterval);          
         
    PkBinningCalc(Num_ModesInInterval, PkArray);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)     fprintf(output, "%d \t %e \t %e\n", modesPerBin[j], meanKBin[j], binnedPk[j]);
    fclose(output);
    
    return 0;
}
