int analyticConvTest(){
    printf("\n\nImplementing anisotropic window fn. analytic test case.");

    prepAnisoConvolution();

    pt2Pk                       = &Analytic2powerlaw;
    pt2AnisoWf                  = &anisokGauss;
    
    assignAnisoWfKernel();
    
    xminKernelshift              = -(xwfKernelsize-1)/2;
    xmaxKernelshift              =  (xwfKernelsize-1)/2;
    
    yminKernelshift              = -(ywfKernelsize-1)/2;
    ymaxKernelshift              =  (ywfKernelsize-1)/2;
    
    zminKernelshift              = -(zwfKernelsize-1)/2;
    zmaxKernelshift              =  (zwfKernelsize-1)/2;
    
    prep_wfKernelminAmp();
    
    // Largest size of the kernel, in k.
    KernelMaxk = fmax(fmax(xmaxKernelshift*kIntervalx, ymaxKernelshift*kIntervaly), zmaxKernelshift*kIntervalz);

    setInputPk();

    // SetWfKernel();
  
    SetWfKernel_minAmp();
    
    ConvNorm = minAmp_ConvolutionNorm(windowFunc3D);
  
    printf("\nConvolution norm. is %e ", ConvNorm);
  
    convolve3DInputPk(convolvedPk3d, inputPk);


    sprintf(filepath, "%s/Data/Del2k/midK_Pk_ConvolvedAnisoGaussTestCase2.dat", root_dir);

    output = fopen(filepath, "w");
  
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %d \n", midKBin[j], binnedPk[j], modesPerBin[j]);
  
    fclose(output);  
      
    return 0;
}

/*
int analyticConvTest_MeasuredWf(){
  printf("\n\nImplementing anisotropic window fn. analytic test case, measured Window fn.");

  prepAnisoConvolution();

  pt2Pk                       = &Analytic2powerlaw;

  minKernelshift              = -0.5*(wfKernelsize-1);
  maxKernelshift              =  0.5*(wfKernelsize-1);

  setInputPk();

  setMeasuredWfKernel();

  convolve3DInputPk(convolvedPk3d, inputPk);

  flattenConvolvedPk();

  sprintf(filepath, "%s/Data/Del2k/midK_Pk_ConvolvedAnisoGaussTestCase_MeasuredWf.dat", root_dir);

  output = fopen(filepath, "w");
  
  for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
  
  fclose(output);

  printWfKernel();
    
  return 0;
}
*/
