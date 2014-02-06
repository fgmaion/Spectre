int MockAvg2Dpk(int mockNumber){
    clean2Dpk();
    
    double** input2Dpk;
    
    input2Dpk                                      = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    
    for(j=0; j<loskBinNumb-1; j++) input2Dpk[j]    = (double  *) malloc((perpkBinNumb-1)*sizeof(double));
    
    for(i=1; i<mockNumber+1; i++){
        printf("\n %d", i);
    
        if(i<10)  sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_00%d.dat", root_dir, surveyType, i);
        else      sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed2Dpk_%s_0%d.dat",  root_dir, surveyType, i);
    
        inputfile = fopen(filepath, "r");

        for(j=0; j<loskBinNumb-1; j++){
            for(k=0; k<perpkBinNumb-1; k++){
                fscanf(inputfile, "%*g \t %*g \t %le \t %*d \n", &input2Dpk[j][k]);
                
                printf("\n %g", input2Dpk[j][k]);
    
                zSpaceBinnedPk[j][k] += (1./mockNumber)*input2Dpk[j][k];
            }
        }

        fclose(inputfile);    
    }
    
    sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/MockAvgObserved2Dpk_%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");

    for(j=0; j<loskBinNumb-1; j++){
        for(k=0; k<perpkBinNumb-1; k++){
            fprintf(output, "%le \n", zSpaceBinnedPk[j][k]);
        }
    }

    fclose(output);
    
    return 0;
}
