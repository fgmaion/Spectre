int prepNGP(){
    densityArray            =  (double *)  realloc(densityArray,     n0*n1*n2*sizeof(*densityArray));
    meanCellRedshift        =  (double *)  realloc(meanCellRedshift, n0*n1*n2*sizeof(*meanCellRedshift));  

    FKPweights              =  (double *)  realloc(FKPweights,       n0*n1*n2*sizeof(*FKPweights));

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
    
    muIntervalPk       = (double **)     malloc(n0*n1*n2*sizeof(*muIntervalPk));
    
    for(j=0; j<n0*n1*n2; j++){  
        PkArray[j]     = (double *)      malloc(2*sizeof(double));                     // columns 
        
        PkArray[j][0]  = 0.0;
        PkArray[j][1]  = 0.0;
    
        muIntervalPk[j]    = (double *)     malloc(2*sizeof(double));
        
        muIntervalPk[j][0] = 0.0;
        muIntervalPk[j][1] = 0.0; 
    }

    printf("\nCreating FFTw plan.");
    p                  = fftw_plan_dft_3d(n0, n1, n2, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    iplan              = fftw_plan_dft_3d(n0, n1, n2, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    return 0;
}


int prepFFTw1D(int n0){
    in1D               = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0);
    out1D              = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0);
    
    p1D                = fftw_plan_dft_1d(n0, in1D, out1D, FFTW_FORWARD, FFTW_ESTIMATE);

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


int assign2DPkMemory(int muBinNumb, int kBinNumb){
    
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
    
    // legendre2weights                                     = (double *)  malloc(n0*n1*n2*sizeof(*legendre2weights));

    //zSpaceBinnedPk                                       = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    //for(j=0; j<loskBinNumb-1; j++) zSpaceBinnedPk[j]     = (double  *) malloc((perpkBinNumb-1)*sizeof(double));
    
    //zSpacemodesPerBin                                    = (int **)    malloc((loskBinNumb-1)*sizeof(int*));
    //for(j=0; j<loskBinNumb-1; j++) zSpacemodesPerBin[j]  = (int  *)    malloc((perpkBinNumb-1)*sizeof(int));
    
    mean_perpk                                           = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    for(j=0; j<loskBinNumb-1; j++) mean_perpk[j]         = (double *)  malloc((perpkBinNumb-1)*sizeof(double));
    
    mean_losk                                            = (double **) malloc((loskBinNumb-1)*sizeof(double*));
    for(j=0; j<loskBinNumb-1; j++) mean_losk[j]          = (double  *) malloc((perpkBinNumb-1)*sizeof(double));

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
    kHexadecapole                                        = (double  *) malloc((kBinNumb-1)*sizeof(*kHexadecapole));

    f_meanKBin                                           = (float  *) malloc((kBinNumb-1)*sizeof(*f_meanKBin));
    f_kMonopole                                          = (float  *) malloc((kBinNumb-1)*sizeof(*kMonopole));
    f_kQuadrupole                                        = (float  *) malloc((kBinNumb-1)*sizeof(*kQuadrupole));

    f_kMonopole2d                                        = (float  *) malloc((kBinNumb-1)*sizeof(*f_kMonopole2d));
    f_kQuadrupole2d                                      = (float  *) malloc((kBinNumb-1)*sizeof(*f_kQuadrupole2d));
    
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
    inputPk                          = (double *)   malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*inputPk));
    
    convolvedPk3d                    = (double *)   malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*convolvedPk3d)); 

    flattenedConvolvedPk3D           = (double **)  malloc((n0+1)*(n1+1)*(n2+1)*sizeof(*flattenedConvolvedPk3D));

    for(j=0; j<(n0+1)*(n1+1)*(n2+1); j++){
        flattenedConvolvedPk3D[j]    = (double *)   malloc(2*sizeof(double));

        flattenedConvolvedPk3D[j][0] = 0.0;
        flattenedConvolvedPk3D[j][1] = 0.0;
    }

    return 0;
}


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
}


int LikelihoodMemory(){
    // Parameter grid of evaluated Chi sq values. 
    ChiSqGrid  = (double ***) malloc(Res*sizeof(*ChiSqGrid));
  
    for(j=0; j<Res; j++)  ChiSqGrid[j] = (double **) malloc(Res*sizeof(**ChiSqGrid));
  
    for(j=0; j<Res; j++){
        for(k=0; k<Res; k++){
            ChiSqGrid[j][k]    = (double *) malloc(Res*sizeof(***ChiSqGrid));
        }
    }

    // Likelihood values. 
    lnLikelihoodGrid  = (double ***) malloc(Res*sizeof(*lnLikelihoodGrid));
  
    for(j=0; j<Res; j++)  lnLikelihoodGrid[j] = (double **) malloc(Res*sizeof(**lnLikelihoodGrid));
  
    for(j=0; j<Res; j++){
        for(k=0; k<Res; k++){
            lnLikelihoodGrid[j][k]    = (double *) malloc(Res*sizeof(***lnLikelihoodGrid));
        }
    }


    betaPosterior = (double *)  malloc(Res*sizeof(*betaPosterior));    
    
    for(i=0; i<Res; i++)  betaPosterior[i] = 0.0;
    
    sigmaPosterior = (double *)  malloc(Res*sizeof(*sigmaPosterior));    

    for(i=0; i<Res; i++)  sigmaPosterior[i] = 0.0;
    
    
    betaSigmaPosterior = (double **) malloc(Res*sizeof(*betaSigmaPosterior));
    
    for(j=0; j<Res; j++)  betaSigmaPosterior[j] = malloc(Res*sizeof(**betaSigmaPosterior));
    
    
    for(i=0; i<Res; i++){
        for(j=0; j<Res; j++){
            betaSigmaPosterior[i][j] = 0.0;
        }
    }
    
    return 0;
}


int assignCovMat(int mockNumber, int kBinNumb, int hiMultipoleOrder){
    // Multipoles is a [CatalogNumber][hiMultipoleOrder][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
    Multipoles                                             = (double ***) malloc(mockNumber*sizeof(*Multipoles));
        
    for(j=0; j<mockNumber; j++) Multipoles[j]              = (double  **) malloc(hiMultipoleOrder*sizeof(**Multipoles));
    
    for(j=0; j<mockNumber; j++){
        for(k=0; k<hiMultipoleOrder; k++){      
            Multipoles[j][k]                               = (double   *) malloc(kBinNumb*sizeof(***Multipoles));
        }
    }

    for(i=0; i<mockNumber; i++){
        for(j=0; j<hiMultipoleOrder; j++){
            for(k=0; k<kBinNumb; k++){
                Multipoles[i][j][k] = 0.0;
            }
        }
    }
    
    kMultipoles                                            = (double  *) malloc(kBinNumb*sizeof(*kMultipoles));

    ModeNumber                                             = (int     *) malloc(kBinNumb*sizeof(*ModeNumber));

    MeanMultipoles                                         = (double **) malloc(hiMultipoleOrder*sizeof(*MeanMultipoles));    
    
    for(j=0; j<hiMultipoleOrder; j++)  MeanMultipoles[j]   = (double  *) malloc(kBinNumb*sizeof(**MeanMultipoles));
    
    for(j=0; j<hiMultipoleOrder; j++){
        for(k=0; k<kBinNumb; k++){
            MeanMultipoles[j][k] = 0.0;
        }
    }
    
    
    // Covariance is an N x N matrix, where N corresponds to hiMultipoleOrder*(kBinNumb-1), here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
         
    Covariance                                                      = (double **) malloc(2*kBinNumb*sizeof(*Covariance));
    
    // Covariance_CofactorMatrix                                    = (double **) malloc(2*(kBinNumb-1)*sizeof(*Covariance_CofactorMatrix));
    
    for(j=0; j<2*kBinNumb; j++) Covariance[j]                       = (double  *) malloc(2*kBinNumb*sizeof(**Covariance));
    // for(j=0; j<2*(kBinNumb-1); j++) Covariance_CofactorMatrix[j] = (double  *) malloc(2*(kBinNumb-1)*sizeof(**Covariance_CofactorMatrix));
    
    // flatCov     = malloc(2*(kBinNumb-1)*2*(kBinNumb-1)*sizeof(*flatCov));
    
    for(k=0; k<2*kBinNumb; k++){
        for(j=0; j<2*kBinNumb; j++){
            Covariance[j][k]                = 0.0;
            // Covariance_CofactorMatrix[j][k] = 0.0;
        }
    }
    
    /*
    invCov      = (double **) malloc(2*(kBinNumb-1)*sizeof(*invCov));
    
    for(j=0; j<2*(kBinNumb-1); j++){
        invCov[j] = malloc(2*(kBinNumb-1)*sizeof(**invCov));
    }
    */
    
    return 0;
}


int ClippingModelling(){
    PkCube             = (double *) realloc(PkCube,          n0*n1*n2*sizeof(*PkCube));
    
    Corrfn            = (double *)  realloc(Corrfn,             n0*n1*n2*sizeof(*Corrfn));
    suppressedCorrfn  = (double *)  realloc(suppressedCorrfn,   n0*n1*n2*sizeof(*suppressedCorrfn));
    distortedCorrfn   = (double *)  realloc(distortedCorrfn,    n0*n1*n2*sizeof(*distortedCorrfn));
    clippedPk         = (double *)  realloc(clippedPk,          n0*n1*n2*sizeof(*clippedPk));

    return 0;
}
