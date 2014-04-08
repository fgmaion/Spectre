int printWindowfuncSlices(){
    printf("\nPrinting window func. slices.");
    
    sprintf(WindowfuncSlices_dir, "%s/Data/WindowfuncSlices", root_dir);

    sprintf(Windowfunc_xSlice, "%s/%s_xSlice_NGPcorrected.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervalx, n2, 1, 0, 0, Windowfunc_xSlice);
    
    sprintf(Windowfunc_ySlice, "%s/%s_ySlice_NGPcorrected.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervaly, n1, 0, 1, 0, Windowfunc_ySlice);
    
    sprintf(Windowfunc_zSlice, "%s/%s_zSlice_NGPcorrected.dat", WindowfuncSlices_dir, surveyType);
    WindowfuncSlice(kIntervalz, n0, 0, 0, 1, Windowfunc_zSlice);

    return 0;
}


int WindowfuncSlice(double kintervali, int ni, int x0, int y0, int z0, char filepath[]){
    double k_i                 = 0.0;
    double Pk                  = 0.0;

    output = fopen(filepath, "w");

    for(j=0; j<ni/2 + 1; j++){
        k_i = kintervali*j;
        
        Index                              = n1*n2*z0*j + n2*y0*j + x0*j;

        H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
        H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];

        WindowFunc                         = sin(pi*k_i*0.5/NyquistWaveNumber)/(pi*k_i*0.5/NyquistWaveNumber);
        
        // NGP correction.
        H_kReal                           /= WindowFunc;
        H_kImag                           /= WindowFunc;

	    Pk                                 = pow(H_kReal, 2.) + pow(H_kImag, 2.);

        // Shot noise subtraction. Not binned estimate, not necessarily expectation. 
        // Pk                                -= 1./rand_number;

        // Correct for window fn. normalisation. 
        Pk                                /= pow(fkpWeightedVolume/TotalVolume, 2.);
        
        fprintf(output, "%.9e \t %.9e\n", log10(k_i), log10(Pk));
    }
         
    fclose(output);
    
    return 0;
}
