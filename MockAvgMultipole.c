int MockAvgMultipole(int mockNumber){
    double* inputMultipole;
    
    double* mockAvgMonopole;
    // double* mockAvgQuadrupole;
    
    inputMultipole                                 = (double *) malloc((kBinNumb-1)*sizeof(double));
    
    mockAvgMonopole                                = (double *) malloc((kBinNumb-1)*sizeof(double));
    // mockAvgQuadrupole                              = (double *) malloc((kBinNumb-1)*sizeof(double));
        
    for(i=1; i<mockNumber+1; i++){
        printf("\n %d", i);
    
        if(i<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles%s_kbin_%.2f_00%d.dat", root_dir, surveyType, kbinInterval, i);
        else      sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles%s_kbin_%.2f_0%d.dat",  root_dir, surveyType, kbinInterval, i);
    
        inputfile = fopen(filepath, "r");

        for(j=0; j<kBinNumb-1; j++){
            fscanf(inputfile, "%*g \t %*g \t %le \t %*d \n", &inputMultipole[j]);
    
            mockAvgMonopole[j] += (1./mockNumber)*inputMultipole[j];
                
        }

        fclose(inputfile);    
    } 
    
    
    sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/MockAvgObservedMonopole_%s_kbin_%.2f.dat", root_dir, surveyType. kbinInterval);
    
    output = fopen(filepath, "w");

    for(j=0; j<kBinNumb-1; j++){
        fprintf(output, "%le \n", mockAvgMonopole[j]);
    }

    fclose(output);
    
    return 0;
}
