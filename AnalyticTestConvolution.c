int analyticConvTest(){
  printf("\n\nImplementing anisotropic window fn. analytic test case.");

  prepAnisoConvolution();

  pt2Pk                       = &Analytic2powerlaw;
  pt2AnisoWf                  = &anisokGauss;

  minKernelshift              = -0.5*(wfKernelsize-1);
  maxKernelshift              =  0.5*(wfKernelsize-1);

  setInputPk();

  SetWfKernel();
  
  convolve3DInputPk(convolvedPk3d, inputPk);

  flattenConvolvedPk();

  sprintf(filepath, "%s/Data/Del2k/midK_Pk_ConvolvedAnisoGaussTestCase.dat", root_dir);

  output = fopen(filepath, "w");
  
  for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
  
  fclose(output);
  
  printWfKernel();
  
  return 0;
}


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
