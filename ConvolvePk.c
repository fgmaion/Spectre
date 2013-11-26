int ConvolveTheory(){
    printf("\n\nConvolving HOD theory P(k).");
    
    // HOD theory P(k).
    inputPK();

    muInterval   = (1.0 - (-1.0))/MuIntegralPrecision;
    
    qSpaceVolume = (4./3.)*pi*pow(0.2, 3.);
    printf("\nq space volume:  %f\n", qSpaceVolume);
    
    // Calculate the convolution of HOD P(k) with the upadded measurement of the window fn.
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2k_%s.dat", root_dir, surveyType);
    splineInputWindowfn(filepath);

    EvaluateConvolution(splintWindowfunc, "splint_unpadded");
    freeConvolutionMemory();


    // Calculate the convolution of the HOD P(k) with the padded measuremen of the window fn.
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_pad3W2k_%s.dat", root_dir, surveyType);
    splineInputWindowfn(filepath);

    EvaluateConvolution(splintWindowfunc, "splint_padded");
    
    // Now calculate the integral constraint correction.
    IntegralConstraintCorrection();
    
    // Now the analytic "truth" for the spherical case.
    EvaluateConvolution(AnalyticSpherical, "splintanalytic");
    
    // printInterpPk();
    
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


int EvaluateConvolution(float (*EvaluateWindowfunc)(float), char type[]){
    // splint the measured window fn.
    FilterNormalisation = EvaluatefilterNormalisation((*EvaluateWindowfunc));

    // For each k mode interval in convolved P(k).
    for(k=0; k<kBinNumb-1; k++)  ConvolvedPk[k] = ConvolvedPkCalculation(midKBin[k], splintMatterPk, (*EvaluateWindowfunc));

    sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Del2k_%s_%s.dat", root_dir, surveyType, type);

    output = fopen(filepath, "w");    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], pow(2.*pi*pi, -1.)*ConvolvedPk[j]*pow(midKBin[j], 3.));
    fclose(output);

    return 0;
}


float EvaluatefilterNormalisation(float (*EvaluateWindowfunc)(float)){
    // Interpolated Window fn. evaluation (spherical average of).
    
    float Interim = 0.0;
    
    for(j=0; j<100; j++){
        WindowFuncEval       = (*EvaluateWindowfunc)(kVals[j]);
        Interim             += WindowFuncEval*kVals[j]*kVals[j]*(4.*pi)*(1./qSpaceVolume)*kbinInterval;
    }

    return Interim;
}


double ConvolvedPkCalculation(float kval, float (*EvaluateMatterPk)(float), float (*EvaluateWindowfunc)(float)){
    double result     =  0.0;
    muVal             = -1.0;

    // For each interval in delta mu_q.
    for(j=0; j<MuIntegralPrecision; j++){  
        
        // For each q mode interval.  
        for(q=0; q<100; q++){
            // Window fn. evaluation point. np.sqrt(k^2 + q^2 -2*k*q*\mu)
            PkEvalPoint    = pow(pow(kval, 2.) + pow(kVals[q], 2.) - 2.0*kval*kVals[q]*muVal, 0.5); 
            
            matterPkEval   = (*EvaluateMatterPk)(PkEvalPoint);
                
            WindowFuncEval = (*EvaluateWindowfunc)(kVals[q]);
                
            result += matterPkEval*pow(kVals[q], 2.0)*WindowFuncEval;
        }       
    
        muVal += muInterval;
    }

    result    *= 2.0*pi*muInterval*kbinInterval/(qSpaceVolume*FilterNormalisation);

    return result;
}


double ConvolvedPkZeroPointCalc(float (*EvaluateMatterPk)(float), float (*EvaluateWindowfunc)(float)){
    double result   = 0.0;
    
    // For each q mode interval.  
    for(q=0; q<100; q++){
        // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
        matterPkEval   = (*EvaluateMatterPk)(kVals[q]);
                
        // Interpolated Window fn. evaluation (spherical average of).
        WindowFuncEval = (*EvaluateWindowfunc)(kVals[q]);
                
        result += matterPkEval*pow(kVals[q], 2.0)*WindowFuncEval;
    }

    result     *= 4.*pi*kbinInterval/(qSpaceVolume*FilterNormalisation);

    return result;
}


int IntegralConstraintCorrection(){    
    WindowfnZeroPointEval    = WindowFuncNR[1];
    printf("\nBinned   window fn. zero point:  %e", WindowFuncNR[1]);
    
    ConvolvedPkZeroPoint     = ConvolvedPkZeroPointCalc(splintMatterPk, splintWindowfunc);
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


float splintMatterPk(float EvalPoint){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 

    float Interim;
    
    splint(sdltk, sdltPk, sdlt2d, 293, EvalPoint, &Interim);
    return Interim;
}


float AnalyticSpherical(float k){
    float Interim;
    float y = k*250.;

    Interim = 3.*pow(y, -3.)*(sin(y) - y*cos(y));
    
    return Interim;
}


int printInterpPk(){
    printf("\nWriting regular spaced, interpolated HOD theory P(k).");
    
    sprintf(filepath, "%s/Data/InterpTheoryPk/regInterpCambExtendedPk_%.4f_hod_20.0.dat", root_dir, kbinInterval);
    output = fopen(filepath, "w");
    for(j=0; j<InterpK_binNumber; j++)  fprintf(output, "%f \t %f\n", kVals[j], interpolatedPk[j]);
    fclose(output);
    
    return 0;
}
