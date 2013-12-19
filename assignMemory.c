int prepNGP(){
    densityArray            =  (double *)  realloc(densityArray, n0*n1*n2*sizeof(*densityArray));
    FKPweights              =  (double *)  realloc(FKPweights,   n0*n1*n2*sizeof(*FKPweights));
    booldensity             =  (double *)  realloc(booldensity,  n0*n1*n2*sizeof(*booldensity));

    Cell_rotatedXvals       =  (double *)  realloc(Cell_rotatedXvals,       n0*n1*n2*sizeof(*Cell_rotatedXvals));
    Cell_rotatedYvals       =  (double *)  realloc(Cell_rotatedYvals,       n0*n1*n2*sizeof(*Cell_rotatedYvals));
    Cell_rotatedZvals       =  (double *)  realloc(Cell_rotatedZvals,       n0*n1*n2*sizeof(*Cell_rotatedZvals));

    Cell_raVIPERSsystem     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_raVIPERSsystem));
    Cell_decVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_decVIPERSsystem));
    Cell_chiVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_chiVIPERSsystem));

    Cell_VIPERSweights      =  (double  *) malloc(n0*n1*n2*sizeof(*Cell_VIPERSweights));
    Cell_VIPERSbools        =  (double  *) malloc(n0*n1*n2*sizeof(*Cell_VIPERSbools));
    Cell_SurveyLimitsMask   =  (double  *) malloc(n0*n1*n2*sizeof(*Cell_SurveyLimitsMask));

    // Allocate memory for NGP density arrays of both ZADE galaxies and randoms. 
    for(j=0; j<n0*n1*n2; j++)            FKPweights[j] = 1.0;
    for(j=0; j<n0*n1*n2; j++)          densityArray[j] = 0.0;
    for(j=0; j<n0*n1*n2; j++)      Cell_VIPERSbools[j] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  Cell_SurveyLimitsMask[j] =0.0;
    
    return 0;
}


int prepFFTw(int n0, int n1, int n2){    
    in                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
    out                = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2); 

    // PkArray[row][column]. row < n0*n1*n2. column < 2 {Real, Imaginary}.
    printf("\n\nAssigning Pk array.");
    PkArray            = (double **)     malloc(n0*n1*n2*sizeof(double*));             // rows
    
    for(j=0; j<n0*n1*n2; j++){  
        PkArray[j]     = (double *)      malloc(2*sizeof(double));                     // columns 
        
        PkArray[j][0]  = 0.0;
        PkArray[j][1]  = 0.0;
    }

    printf("\nCreating FFTw plan.");
    p                  = fftw_plan_dft_3d(n0, n1, n2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    return 0;
}


int prepFFTbinning(){    
    // Memory assigned for Pk binning.  
    kBinLimits        = (double *) malloc(kBinNumb*sizeof(*kBinLimits));

    del2              = (double *) malloc((kBinNumb-1)*sizeof(*del2));
    midKBin           = (double *) malloc((kBinNumb-1)*sizeof(*midKBin));
    meanKBin          = (double *) malloc((kBinNumb-1)*sizeof(*meanKBin));
    modesPerBin       = (int *)    malloc((kBinNumb-1)*sizeof(*modesPerBin));
    binnedPk          = (double *) malloc((kBinNumb-1)*sizeof(*binnedPk));
    linearErrors      = (double *) malloc((kBinNumb-1)*sizeof(*linearErrors));
    
    return 0;
}


int assign2DPkMemory(){
    TwoDpkArray           = (double **) malloc(n0*n1*n2*sizeof(double*));                         

    for(j=0; j<n0*n1*n2; j++){
        TwoDpkArray[j]    = (double *)  malloc(3*sizeof(double));
        
        TwoDpkArray[j][0] = 0.0; 
        TwoDpkArray[j][1] = 0.0;
        TwoDpkArray[j][2] = 0.0;
    }

    legendre2weights      = (double *)  malloc(n0*n1*n2*sizeof(*legendre2weights));

    zSpaceBinnedPk                                   = (double **) malloc((kBinNumb-1)*sizeof(double*));
    for(j=0; j<kBinNumb-1; j++) zSpaceBinnedPk[j]    = (double  *) malloc((kBinNumb-1)*sizeof(double));
    
    zSpacemodesPerBin                                = (int **)    malloc((kBinNumb-1)*sizeof(int*));
    for(j=0; j<kBinNumb-1; j++) zSpacemodesPerBin[j] = (int  *)    malloc((kBinNumb-1)*sizeof(int));
    
    mean_perpk                                       = (double **) malloc((kBinNumb-1)*sizeof(double*));
    for(j=0; j<kBinNumb-1; j++) mean_perpk[j]        = (double *)  malloc((kBinNumb-1)*sizeof(double));
    
    mean_losk                                        = (double **) malloc((kBinNumb-1)*sizeof(double*));
    for(j=0; j<kBinNumb-1; j++) mean_losk[j]         = (double  *) malloc((kBinNumb-1)*sizeof(double));


    for(j=0; j<kBinNumb-1; j++){
        for(i=0; i<kBinNumb-1; i++){
            zSpaceBinnedPk[j][i]    = 0.0;
            zSpacemodesPerBin[j][i] =   0;
            mean_perpk[j][i]        = 0.0;
            mean_losk[j][i]         = 0.0;
        }
    }

    return 0;
}


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
}


int prepAnisoConvolution(){
    inputPk                          = (double *)   malloc(n0*n1*n2*sizeof(*inputPk));
    
    windowFunc3D                     = (double *)   malloc(wfKernelsize*wfKernelsize*wfKernelsize*sizeof(*windowFunc3D));       
    
    convolvedPk3d                    = (double *)   malloc((n0-2*wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize)*sizeof(*convolvedPk3d)); 

    flattenedConvolvedPk3D           = (double **)  malloc((n0-2*wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize)*sizeof(*flattenedConvolvedPk3D));

    for(j=0; j<(n0-2*wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize); j++){
        flattenedConvolvedPk3D[j]    = (double *)   malloc(2*sizeof(double));

        flattenedConvolvedPk3D[j][0] = 0.0;
        flattenedConvolvedPk3D[j][1] = 0.0;
    }

    return 0;
}


int assignAnisoWfKernel(){
    AnisoWfKernel                    = (double *) malloc(wfKernelsize*wfKernelsize*wfKernelsize*sizeof(*AnisoWfKernel));
    AnisoWfKernel_ModeNumb           = (int *)    malloc(wfKernelsize*wfKernelsize*wfKernelsize*sizeof(*AnisoWfKernel_ModeNumb));

    for(j=0; j<wfKernelsize*wfKernelsize*wfKernelsize; j++){
        AnisoWfKernel[j]             = 0.0;
        AnisoWfKernel_ModeNumb[j]    = 0;
    }

    return 0;
}
