/*
int IntegralConstraintCorrection(){    
    ConvPkZeroPoint     = SphConvZeroPoint();
    printf("\nConvolved P(k) zero point calculated to be: %e", ConvPkZeroPoint);
    
    // sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    // output = fopen(filepath, "w");    
    // fprintf(output, "%f \t %f \n", 0.0, ConvolvedPkZeroPoint);  
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    // fclose(output);
    
    
    // Window fn. evaluated on the same regular grid in k on which the convolved P(k) is estimated, normalised such that W^2 = 1.0 at k=0.0
    // for(j=0; j<kBinNumb-1; j++) ConvolvedPk[j] -= splintWindowfunc(midKBin[j])*ConvolvedPkZeroPoint;

    // sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Corrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    // output = fopen(filepath, "w");    
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    // fclose(output);

    return 0;
}
*/
