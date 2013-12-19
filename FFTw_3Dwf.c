int wfPkCalc(){
    // Delivers the spherical average of W^2(q).

    printf("\nBeginning window function calculation.");
    
    // The true density field is multiplied by a mask, Cell_AppliedWindowFn [W(x)].
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");

    PkCorrections(1);

    PkBinningCalc(n0*n1*n2, PkArray);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2q_%s.dat", root_dir, surveyType);
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%e \t %e\n", meanKBin[j], binnedPk[j]);

    printf("\nWindow function P(k) calculation complete.");
    
    fclose(output);
    
    return 0;
}
