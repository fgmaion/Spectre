int ConvolveTheory(){
    printf("\n\nConvolving HOD theory P(k).");
    
    // HOD theory P(k).
    inputPK();
    
    // Calculate the convolution of HOD P(k) with the upadded measurement of the window fn.
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2k_%s.dat", root_dir, surveyType);
    splineInputWindowfn(filepath);

    // EvaluateConvolution(splintMatterPk, splintWindowfunc, "splint_unpadded");
    // freeConvolutionMemory();

    // Calculate the convolution of the HOD P(k) with the padded measurement of the window fn.
    // sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_pad3W2k_%s.dat", root_dir, surveyType);
    // splineInputWindowfn(filepath);

    // EvaluateConvolution(splintMatterPk, splintWindowfunc, "splint_padded");
    
    // Now calculate the integral constraint correction.
    // IntegralConstraintCorrection();
    
    // Now the analytic "truth" for the spherical case.
    // EvaluateConvolution(AnalyticSpherical, "splintanalytic");
    
    // printInterpPk();
    
    pt2Pk                       = &Analytic2powerlaw;
    pt2Windowfn                 = &AnalyticGaussian;
    
    EvaluateConvolution("powerlaw2_Gauss");
    
    return 0;
}


float EvaluatefilterNormalisation(){
    // Interpolated Window fn. evaluation (spherical average of).
    return qromb(NormalisationKernel, 0.0, 0.9);
}


float NormalisationKernel(float q){
    return 4.*pi*q*q*(*pt2Windowfn)(q);
}


int EvaluateConvolution(char type[]){
    // splint the measured window fn.
    FilterNormalisation = EvaluatefilterNormalisation();
    printf("\nfilter normalisation:  %f", FilterNormalisation);
    
    kConvScope = 1.0;
    
    printf("\nq kernel check:  %f\n", qConvKernel(0.2));

    /*
    // For each k mode interval in convolved P(k).
    for(k=0; k<kBinNumb-1; k++)  ConvolvedPk[k] = ConvolvedPkQrombCalc(midKBin[k]);

    sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_%s_%s.dat", root_dir, surveyType, type);

    output = fopen(filepath, "w");    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    fclose(output);
    */
    return 0;
}


float ConvolvedPkQrombCalc(float kval){
    kConvScope  = kval;
    return qromb(qConvKernel, 0.0, 1.5);
}


float qConvKernel(float q){
    qConvScope = q;
    return qromb(muConvKernel, -1., 1.);
}


float muConvKernel(float mu){
    return  (2.0*pi/FilterNormalisation)*qConvScope*qConvScope*(*pt2Windowfn)(qConvScope)*(*pt2Pk)(pow(pow(kConvScope, 2.) + pow(qConvScope, 2.) - 2.0*kConvScope*qConvScope*mu, 0.5));
}


float ConvolvedPkZeroPointCalc(){
    return ConvolvedPkQrombCalc(0.0);
}


int IntegralConstraintCorrection(){    
    WindowfnZeroPointEval    = WindowFuncNR[1];
    printf("\nBinned   window fn. zero point:  %e", WindowFuncNR[1]);
    
    ConvolvedPkZeroPoint     = ConvolvedPkZeroPointCalc();
    printf("\nConvolved P(k) zero point calculated to be: %e", ConvolvedPkZeroPoint);
    
    sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    output = fopen(filepath, "w");    
    fprintf(output, "%f \t %f \n", 0.0, ConvolvedPkZeroPoint);
    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    fclose(output);
    
    // Window fn. evaluated on the same regular grid in k on which the convolved P(k) is estimated, normalised such that W^2 = 1.0 at k=0.0
    for(j=0; j<kBinNumb-1; j++) ConvolvedPk[j] -= splintWindowfunc(midKBin[j])*ConvolvedPkZeroPoint;

    sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Corrected/midK_Pk_%ssplint_padded.dat", root_dir, surveyType);
    output = fopen(filepath, "w");    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    fclose(output);

    return 0;
}


int splineInputWindowfn(char filepath[]){
    inputfile             = fopen(filepath, "r");
    
    ch                    = 0;
    inputWindowfnBinNumb  = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  inputWindowfnBinNumb += 1;
    } while (ch != EOF);

    prepConvolution(inputWindowfnBinNumb+1, kBinNumb);

    printf("\nWindow fn. line number:  %d", inputWindowfnBinNumb);

    rewind(inputfile);
    
    for(j=1; j<inputWindowfnBinNumb+1; j++) fscanf(inputfile, "%d \t %f \t %f \n", &midKmodesperbin[j], &midKBinNR[j], &WindowFuncNR[j]);
    fclose(inputfile);
                
    spline(midKBinNR, WindowFuncNR, inputWindowfnBinNumb, 1.0e31, 1.0e31, WindowFunc2d);

    return 0;
}


float splintWindowfunc(float EvalPoint){
    // Interpolated Window fn. evaluation (spherical average of).
    
    float Interim;
    
    splint(midKBinNR, WindowFuncNR, WindowFunc2d, inputWindowfnBinNumb, EvalPoint, &Interim);
    return Interim;
}


float splintMatterPk(float k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 

    float Interim;
    
    splint(sdltk, sdltPk, sdlt2d, 293, k, &Interim);
    return Interim;
}


float AnalyticSpherical(float k){
    float y = k*250.;

    return 3.*pow(y, -3.)*(sin(y) - y*cos(y));
}


float AnalyticGaussian(float q){
    float r0 = 10.;

    return expf(-q*q*r0*r0);
}


float Analytic2powerlaw(float k){
    float p0 = 150.;
    float  n =   2.;

    return p0*powf(k, n);
}


int printInterpPk(){
    printf("\nWriting regular spaced, interpolated HOD theory P(k).");
    
    sprintf(filepath, "%s/Data/InterpTheoryPk/regInterpCambExtendedPk_%.4f_hod_20.0.dat", root_dir, kbinInterval);
    output = fopen(filepath, "w");
    for(j=0; j<InterpK_binNumber; j++)  fprintf(output, "%f \t %f\n", kVals[j], interpolatedPk[j]);
    fclose(output);
    
    return 0;
}
