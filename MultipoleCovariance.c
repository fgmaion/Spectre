int CovarianceMatrix(int mockNumber){
  
    double** Covariance;
    
    Covariance                                          = (double **) malloc((kBinNumb-1)*sizeof(double*));
    
    for(j=0; j<kBinNumb-1; j++) Covariance[j]           = (double  *) malloc((kBinNumb-1)*sizeof(double));

    for(j=0; j<kBinNumb-1; j++){
        for(k=0; k<kBinNumb-1; k++){
            Covariance[j][k] = 0.0;
        }
    }

    for(i=1; i<mockNumber+1; i++){
        if(i<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_00%d.dat", root_dir, surveyType, i);
        else      sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_0%d.dat",  root_dir, surveyType, i);
    
        inputfile = fopen(filepath, "r");

        for(j=0; j<loskBinNumb-1; j++){
            for(k=0; k<perpkBinNumb-1; k++){
                fscanf(inputfile, "%*g \t %*g \t %g \t %*d \n", input2Dpk[j][k]);
    
                (1./(double) mockNumber)*(input2Dpk[j][k] - zSpaceBinnedPk[j][k]);
            }
        }

        fclose(inputfile);    
    }

    return 0;
}
