int ConvolveSphericalSymm(){
    printf("\n\nConvolving HOD theory P(k).");
    
    // HOD theory P(k).
    inputPK();
    
    // Calculate the convolution of HOD P(k) with the upadded measurement of the window fn.
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_W2q_%s.dat", root_dir, surveyType);
    splineInputWindowfn(filepath);
    
    pt2Pk                       = &Analytic2powerlaw;
    pt2Windowfn                 = &AnalyticGaussian;

    // EvaluateSphericalConv("splint_unpadded");
    // freeConvolutionMemory();

    // Calculate the convolution of the HOD P(k) with the padded measurement of the window fn.
    // sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_pad3W2k_%s.dat", root_dir, surveyType);
    // splineInputWindowfn(filepath);

    // EvaluateSphericalConv("splint_padded");
    
    // Now calculate the integral constraint correction.
    // IntegralConstraintCorrection();
    
    // Now the analytic "truth" for the spherical case.
    // EvaluateSphericalConv("splintanalytic");
    
    // printInterpPk();
    
    EvaluateSphericalConv("powerlaw2_Gauss");
    
    kConvScope = 0.1;
    
    printf("\n %e", fqConvKernel(0.2));

    return 0;
}


int EvaluateSphericalConv(char Convtype[]){
    ConvNorm = SphericalW2NormEval();
    printf("\nW^2(q) norm (spherical symmetry):  %e", ConvNorm);

    // For each k mode interval in convolved P(k).
    // for(k=0; k<kBinNumb-1; k++)  ConvolvedPk[k] = SphConvPk(midKBin[k]);

    // sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Uncorrected/midK_Pk_%s_%s.dat", root_dir, surveyType, Convtype);

    // output = fopen(filepath, "w");    
    
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%f \t %e \n", midKBin[j], ConvolvedPk[j]);
    
    // fclose(output);
    return 0;
}


double NormKernel(double q){
    return 4.*pi*q*q*(*pt2Windowfn)(q);
}


float fNormKernel(float q){
    double qdub = (double) q;
    
    return (float) NormKernel(qdub);
}


double SphericalW2NormEval(){
    // Interpolated Window fn. evaluation (spherical average of).
    return (double) qromb(fNormKernel, 0.0, 0.9);
}


double SphConvPk(double kval){
    kConvScope  = (float) kval;
    
    return (double) qromb(fqConvKernel, 0.0, 4.0);
}

// Evaluated correctly.
float fqConvKernel(float q){
    double qdub = (double) q;

    return (float) qConvKernel(qdub);
}


double qConvKernel(double q){
    qConvScope = (float) q;
    
    return (double) qromb(fmuConvKernel, -1., 1.);
}


float fmuConvKernel(float q){
    double qdub = (double) q;
    
    return (float)  muConvKernel(qdub);
}


double muConvKernel(double mu){
    return  (2.0*pi/ConvNorm)*qConvScope*qConvScope*(*pt2Windowfn)(qConvScope)*(*pt2Pk)(pow(pow(kConvScope, 2.) + pow(qConvScope, 2.) - 2.0*kConvScope*qConvScope*mu, 0.5));
}


double SphConvZeroPoint(){
    return (double) SphConvPk(0.0);
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

    prepSphericalConv(inputWindowfnBinNumb+1, kBinNumb);

    printf("\nWindow fn. line number:  %d", inputWindowfnBinNumb);

    rewind(inputfile);
    
    for(j=1; j<inputWindowfnBinNumb+1; j++) fscanf(inputfile, "%d \t %f \t %f \n", &midKmodesperbin[j], &midKBinNR[j], &WindowFuncNR[j]);
    fclose(inputfile);
                
    spline(midKBinNR, WindowFuncNR, inputWindowfnBinNumb, 1.0e31, 1.0e31, WindowFunc2d);

    return 0;
}


double splintWindowfunc(double EvalPoint){
    // Interpolated Window fn. evaluation (spherical average of).
    
    float Interim;
    
    splint(midKBinNR, WindowFuncNR, WindowFunc2d, inputWindowfnBinNumb, (float) EvalPoint, &Interim);
    return (double) Interim;
}


double splintMatterPk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 

    float Interim;
    
    splint(sdltk, sdltPk, sdlt2d, 293, (float) k, &Interim);
    return (double) Interim;
}


int printInterpPk(){
    printf("\nWriting regular spaced, interpolated HOD theory P(k).");
    
    sprintf(filepath, "%s/Data/InterpTheoryPk/regInterpCambExtendedPk_%.4f_hod_20.0.dat", root_dir, kbinInterval);
    output = fopen(filepath, "w");
    for(j=0; j<InterpK_binNumber; j++)  fprintf(output, "%f \t %f\n", kVals[j], interpolatedPk[j]);
    fclose(output);
    
    return 0;
}
