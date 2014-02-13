int ComputeMixingMatrix(){
    // Convolves input theory Monopole with window fn.
    double MixingMatrix[kBinNumb][kBinNumb];
    
    double unitTheory_lowk;
    double unitTheory_hik;

    MeasureAnisoWfKernel();

    prepAnisoConvolution();

    pt2Pk = &unitTheoryVector;
    
    printf("\n\nCalculating mixing matrix.");

    setMeasuredWfKernel();

    // Calculate matrix coefficient, M_ij.
    for(ii=0; ii<kBinNumb; ii++){
        unitTheory_lowk = kBinLimits[i];
        unitTheory_hik  = kBinLimits[i+1];
        
        setInputPk();

        convolve3DInputPk(convolvedPk3d, inputPk);

        }
    }
    inputPK();

    return 0;
}


double unitTheoryVector(double k){
    // Return a unit theory vector into which the matter P(k) will be decomposed. 
    // i.e return 1. if k lies between unitTheory_hik and unitTheory_lowk.
    
    if((unitTheory_lowk<k) && (k<unitTheory_hik)){
        return 1.;
    }
    
    else{
        return 0.0;
    }
}
