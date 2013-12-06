int wfPkCalc(){
    printf("\nBeginning window function calculation.");
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = (double) (TotalVolume/TotalSurveyedVolume)*booldensity[j]*FKPweights[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = (double) 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");

    PkCorrections(1);

    PkBinningCalc(n0*n1*n2, PkArray);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2k_%s.dat", root_dir, surveyType);
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%f \t %f\n", meanKBin[j], binnedPk[j]);

    printf("\nWindow function P(k) calculation complete.");
    
    fclose(output);
    
    return 0;
}