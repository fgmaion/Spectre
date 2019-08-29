int prepNGP(){
    // densityArray               =  (double *)  realloc(densityArray,     n0*n1*n2*sizeof(*densityArray));
    // meanCellRedshift           =  (double *)  realloc(meanCellRedshift, n0*n1*n2*sizeof(*meanCellRedshift));  

    // FKPweights                 =  (double *)  realloc(FKPweights,       n0*n1*n2*sizeof(*FKPweights));

    // Cell_rotatedXvals       =  (double *)  realloc(Cell_rotatedXvals,       n0*n1*n2*sizeof(*Cell_rotatedXvals));
    // Cell_rotatedYvals       =  (double *)  realloc(Cell_rotatedYvals,       n0*n1*n2*sizeof(*Cell_rotatedYvals));
    // Cell_rotatedZvals       =  (double *)  realloc(Cell_rotatedZvals,       n0*n1*n2*sizeof(*Cell_rotatedZvals));

    // Cell_raVIPERSsystem     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_raVIPERSsystem));
    // Cell_decVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_decVIPERSsystem));
    // Cell_chiVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_chiVIPERSsystem));

    // Cell_VIPERSweights      =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_VIPERSweights));
    // Cell_VIPERSbools        =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_VIPERSbools));
    // Cell_SurveyEdge         =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_SurveyEdge));
    // Cell_ApodiseWeights     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_ApodiseWeights));
    // Cell_ShortDist2edge     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_ShortDist2edge));

    // Allocate memory for NGP density arrays.
    // for(j=0; j<n0*n1*n2; j++)             FKPweights[j] =                  1.0;
    // for(j=0; j<n0*n1*n2; j++)              densityArray[j] =                  0.0;
    // for(j=0; j<n0*n1*n2; j++)        Cell_SurveyEdge[j] =                  0.0;
    // for(j=0; j<n0*n1*n2; j++)          Cell_VIPERSbools[j] =                  0.0;
    // for(j=0; j<n0*n1*n2; j++)       meanCellRedshift[j] =                  0.0;
    // for(j=0; j<n0*n1*n2; j++)    Cell_ShortDist2edge[j] =  10.0*GibbsSkinDepth;
    // for(j=0; j<n0*n1*n2; j++)    Cell_ApodiseWeights[j] =                  1.0;
    
    Cell_SurveyLimitsMask      =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_SurveyLimitsMask));
    
    for(j=0; j<n0*n1*n2; j++)     Cell_SurveyLimitsMask[j] = 0.0;
    
    return 0;
}


/*
int prepSphericalConv(int inputWindowfnBinNumb, int evalConvPkNumb){    
    // Define arrays for Window function spline, with the usual NR offset. 
    midKmodesperbin  = (int   *)   malloc(inputWindowfnBinNumb*sizeof(int));
    midKBinNR        = (float *)   malloc(inputWindowfnBinNumb*sizeof(float)); 
    WindowFuncNR     = (float *)   malloc(inputWindowfnBinNumb*sizeof(float)); 
    WindowFunc2d     = (float *)   malloc(inputWindowfnBinNumb*sizeof(float)); 
    
    ConvolvedPk      = (double *)  malloc((evalConvPkNumb-1)*sizeof(double));
    ConvolvedPk2d    = (double *)  malloc((evalConvPkNumb-1)*sizeof(double));
    
    for(j=0; j<evalConvPkNumb-1; j++)    ConvolvedPk[j] = 0.0;
    
    return 0;
}*/

/*
int prepAnisoConvolution(){
    inputPk                          = (double *)   malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*inputPk));
    
    convolvedPk3d                    = (double *)   malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*convolvedPk3d)); 

    flattenedConvolvedPk3D           = (double **)  malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*flattenedConvolvedPk3D));

    for(j=0; j<(n0+1)*(n1+1)*(n2+1); j++){
        flattenedConvolvedPk3D[j]    = (double *)   malloc(2*sizeof(double));

        flattenedConvolvedPk3D[j][0] = 0.0;
        flattenedConvolvedPk3D[j][1] = 0.0;
    }

    return 0;
}*/

/*
int prep_wfKernelminAmp(){
    wfKernel_minAmpIndices            = (int **)   malloc(largeAmpIndices*sizeof(*wfKernel_minAmpIndices));
    
    windowFunc3D                      = (double *) malloc(largeAmpIndices*sizeof(*windowFunc3D));       
    
    for(j=0; j<largeAmpIndices; j++){  
        wfKernel_minAmpIndices[j]     = (int *)  malloc(4*sizeof(int));                     // columns 
        
        wfKernel_minAmpIndices[j][0]  = 0;
        wfKernel_minAmpIndices[j][1]  = 0;
        wfKernel_minAmpIndices[j][2]  = 0;
        wfKernel_minAmpIndices[j][3]  = 0;
    }

    return 0;
}*/



/*
int ClippingModelling(){
    PkCube            = (double *) realloc(PkCube,          n0*n1*n2*sizeof(*PkCube));
    
    Corrfn            = (double *)  realloc(Corrfn,             n0*n1*n2*sizeof(*Corrfn));
    suppressedCorrfn  = (double *)  realloc(suppressedCorrfn,   n0*n1*n2*sizeof(*suppressedCorrfn));
    distortedCorrfn   = (double *)  realloc(distortedCorrfn,    n0*n1*n2*sizeof(*distortedCorrfn));
    clippedPk         = (double *)  realloc(clippedPk,          n0*n1*n2*sizeof(*clippedPk));
    
    rmodulus_vec      =             realloc(rmodulus_vec,       (1 + n0*n1*n2)*sizeof(*rmodulus_vec));
    
    for(j=0; j<=n0*n1*n2; j++){ 
        rmodulus_vec[j] = malloc(2*sizeof(double));
    }
    
    monoCorr          = (double *)  realloc(monoCorr,           4096*sizeof(*monoCorr));
    rmonoCorr         = (double *)  realloc(rmonoCorr,          4096*sizeof(*rmonoCorr));
    monoCorr2d        = (double *)  realloc(monoCorr2d,         4096*sizeof(*monoCorr2d));

    return 0;
}*/

/*
int assign_mixingmatrix(int kBinNumb, double kMonopole[], double kQuadrupole[]){
    theta                     = malloc(2*(kBinNumb-1)*sizeof(double));
        
    for(j=0; j<kBinNumb-1; j++){  
        theta[j]              =   kMonopole[j];   
        theta[kBinNumb-1 + j] = kQuadrupole[j];
    }
    
    mixingmatrix              = malloc(2*(kBinNumb-1)*sizeof(*mixingmatrix));
    
    for(j=0; j<2*(kBinNumb-1); j++){  
        mixingmatrix[j]       = malloc(2*(kBinNumb-1)*sizeof(double));

        for(k=0; k<2*(kBinNumb-1); k++){
            mixingmatrix[j][k] = 0.0;
        }
    }
    
    return 0;
}*/


/*
int prepFFTw1D(int n0){
    in1D               = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0);
    out1D              = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0);
    
    p1D                = fftw_plan_dft_1d(n0, in1D, out1D, FFTW_FORWARD, FFTW_ESTIMATE);

    return 0;
}*/
