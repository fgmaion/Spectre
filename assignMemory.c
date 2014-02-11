int prepNGP(){
    densityArray            =  (double *)  realloc(densityArray,     n0*n1*n2*sizeof(*densityArray));
    meanCellRedshift        =  (double *)  realloc(meanCellRedshift, n0*n1*n2*sizeof(*meanCellRedshift));  

    FKPweights              =  (double *)  realloc(FKPweights,       n0*n1*n2*sizeof(*FKPweights));
    booldensity             =  (double *)  realloc(booldensity,      n0*n1*n2*sizeof(*booldensity));

    Cell_rotatedXvals       =  (double *)  realloc(Cell_rotatedXvals,       n0*n1*n2*sizeof(*Cell_rotatedXvals));
    Cell_rotatedYvals       =  (double *)  realloc(Cell_rotatedYvals,       n0*n1*n2*sizeof(*Cell_rotatedYvals));
    Cell_rotatedZvals       =  (double *)  realloc(Cell_rotatedZvals,       n0*n1*n2*sizeof(*Cell_rotatedZvals));

    Cell_raVIPERSsystem     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_raVIPERSsystem));
    Cell_decVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_decVIPERSsystem));
    Cell_chiVIPERSsystem    =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_chiVIPERSsystem));

    Cell_VIPERSweights      =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_VIPERSweights));
    Cell_VIPERSbools        =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_VIPERSbools));
    Cell_SurveyLimitsMask   =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_SurveyLimitsMask));
    Cell_SurveyEdge         =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_SurveyEdge));
    Cell_ApodiseWeights     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_ApodiseWeights));
    Cell_ShortDist2edge     =  (double *)  malloc(n0*n1*n2*sizeof(*Cell_ShortDist2edge));

    // Allocate memory for NGP density arrays.
    for(j=0; j<n0*n1*n2; j++)             FKPweights[j] =                  1.0;
    for(j=0; j<n0*n1*n2; j++)           densityArray[j] =                  0.0;
    for(j=0; j<n0*n1*n2; j++)        Cell_SurveyEdge[j] =                  0.0;
    for(j=0; j<n0*n1*n2; j++)       Cell_VIPERSbools[j] =                  0.0;
    for(j=0; j<n0*n1*n2; j++)       meanCellRedshift[j] =                  0.0;
    for(j=0; j<n0*n1*n2; j++)    Cell_ShortDist2edge[j] =  10.0*GibbsSkinDepth;
    for(j=0; j<n0*n1*n2; j++)    Cell_ApodiseWeights[j] =                  1.0;
    for(j=0; j<n0*n1*n2; j++)  Cell_SurveyLimitsMask[j] =                  0.0;
    
    return 0;
}


int CleanNGP(){
    for(j=0; j<n0*n1*n2; j++)             FKPweights[j] = 1.0;
    for(j=0; j<n0*n1*n2; j++)           densityArray[j] = 0.0;
    for(j=0; j<n0*n1*n2; j++)       meanCellRedshift[j] = 0.0;

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


int cleanFFTbinning(){
    for(j=0; j<n0*n1*n2; j++){
        PkArray[j][0] = 0.0;
        PkArray[j][1] = 0.0;
    }
    
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
    polar2Dpk             = (double **) malloc(n0*n1*n2*sizeof(double*));
    TwoDpkArray           = (double **) malloc(n0*n1*n2*sizeof(double*));                         
 
    for(j=0; j<n0*n1*n2; j++){
          polar2Dpk[j]    = (double *)  malloc(3*sizeof(double)); 
        TwoDpkArray[j]    = (double *)  malloc(3*sizeof(double));
    
          polar2Dpk[j][0] = 0.0;
          polar2Dpk[j][1] = 0.0;
          polar2Dpk[j][2] = 0.0;
          
        TwoDpkArray[j][0] = 0.0; 
        TwoDpkArray[j][1] = 0.0;
        TwoDpkArray[j][2] = 0.0;
    }


    legendre2weights                                     = (double *)  malloc(n0*n1*n2*sizeof(*legendre2weights));

    zSpaceBinnedPk                                       = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    for(j=0; j<loskBinNumb-1; j++) zSpaceBinnedPk[j]     = (double  *) malloc((perpkBinNumb-1)*sizeof(double));
    
    zSpacemodesPerBin                                    = (int **)    malloc((loskBinNumb-1)*sizeof(int*));
    for(j=0; j<loskBinNumb-1; j++) zSpacemodesPerBin[j]  = (int  *)    malloc((perpkBinNumb-1)*sizeof(int));
    
    mean_perpk                                           = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    for(j=0; j<loskBinNumb-1; j++) mean_perpk[j]         = (double *)  malloc((perpkBinNumb-1)*sizeof(double));
    
    mean_losk                                            = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    for(j=0; j<loskBinNumb-1; j++) mean_losk[j]          = (double  *) malloc((perpkBinNumb-1)*sizeof(double));

    clean2Dpk();
    
    loskBinLimits                                        = (double *) malloc(loskBinNumb*sizeof(*loskBinLimits));
    perpkBinLimits                                       = (double *) malloc(perpkBinNumb*sizeof(*perpkBinLimits));

    muBinLimits                                          = (double *) malloc(muBinNumb*sizeof(*muBinLimits));

    polar_modesPerBin                                    = (int   **) malloc((muBinNumb-1)*sizeof(int*));
    for(j=0; j<muBinNumb-1; j++)   polar_modesPerBin[j]  = (int    *) malloc((kBinNumb-1)*sizeof(int));

    mean_mu                                              = (double **) malloc((muBinNumb-1)*sizeof(double*));
    for(j=0; j<muBinNumb-1; j++)   mean_mu[j]            = (double  *) malloc((kBinNumb-1)*sizeof(double));

    mean_modk                                            = (double **) malloc((muBinNumb-1)*sizeof(double*));
    for(j=0; j<muBinNumb-1; j++)   mean_modk[j]          = (double  *) malloc((kBinNumb-1)*sizeof(double));

    polar2DBinnedPk                                      = (double **) malloc((muBinNumb-1)*sizeof(double*));
    for(j=0; j<muBinNumb-1; j++)   polar2DBinnedPk[j]    = (double  *) malloc((kBinNumb-1)*sizeof(double));
    
    kMonopole                                            = (double  *) malloc((kBinNumb-1)*sizeof(*kMonopole));
    kQuadrupole                                          = (double  *) malloc((kBinNumb-1)*sizeof(*kQuadrupole));

    return 0;
}


int clean2Dpk(){
    for(j=0; j<loskBinNumb-1; j++){
        for(i=0; i<perpkBinNumb-1; i++){
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
    
    windowFunc3D                     = (double *)   malloc(xwfKernelsize*ywfKernelsize*zwfKernelsize*sizeof(*windowFunc3D));       
    
    convolvedPk3d                    = (double *)   malloc((n0-2*zwfKernelsize)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize)*sizeof(*convolvedPk3d)); 

    flattenedConvolvedPk3D           = (double **)  malloc((n0-2*zwfKernelsize)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize)*sizeof(*flattenedConvolvedPk3D));

    for(j=0; j<(n0-2*zwfKernelsize)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize); j++){
        flattenedConvolvedPk3D[j]    = (double *)   malloc(2*sizeof(double));

        flattenedConvolvedPk3D[j][0] = 0.0;
        flattenedConvolvedPk3D[j][1] = 0.0;
    }

    return 0;
}


int assignAnisoWfKernel(){
    AnisoWfKernel                    = (double *) malloc(xwfKernelsize*ywfKernelsize*zwfKernelsize*sizeof(*AnisoWfKernel));
    AnisoWfKernel_ModeNumb           = (int *)    malloc(xwfKernelsize*ywfKernelsize*zwfKernelsize*sizeof(*AnisoWfKernel_ModeNumb));

    for(j=0; j<xwfKernelsize*ywfKernelsize*zwfKernelsize; j++){
        AnisoWfKernel[j]             = 0.0;
        AnisoWfKernel_ModeNumb[j]    = 0;
    }

    return 0;
}
