int freeHOD(){
// Free dynamic memory allocated for ZADE catalogue + derived parameters.                                                                                                     
  free(id);
  free(ra);
  free(dec);
  free(zcos);
  free(zpec);
  free(zobs);
  free(zphot);
  free(M_B);
  free(type);
  free(csr);
  free(sampling);
  free(sampling35);
  free(flag_Nagoya);
  free(flag_SSPOC);
  free(flag_SSPOC35);
  free(rand_sel);
  
//  for(j=0; j<Vipers_Num; j++){
//    free(pointing[j]);
//    free(quadrant[j]);  
//  }

  //free(pointing);
  //free(quadrant);

  // derived parameters.                                                                                                 
  free(polarAngle);
  free(rDist);
  free(xCoor);
  free(yCoor);
  free(zCoor);
  return 0;
}


int freeNGP(){
 // NGP density arrays for both ZADE galaxies and randoms.                                                              
  free(densityArray);
  free(FKPweights);
  
  free(Cell_rotatedXvals);
  free(Cell_rotatedYvals);
  free(Cell_rotatedZvals);
  
  free(Cell_raVIPERSsystem);
  free(Cell_decVIPERSsystem);
  free(Cell_chiVIPERSsystem);
  
  free(Cell_VIPERSweights);
  free(Cell_VIPERSbools);
  
  free(booldensity);
  return 0;
}


int freeFFTw(){
    fftw_free(in);
    fftw_free(out);
    
    for(j=0; j<n0*n1*n2; j++)  free(PkArray[j]);
    free(PkArray);
    
    fftw_destroy_plan(p);
    return 0;
}


int freeBinning(){

  // binned Pk arrays.                                                                                                   
  free(kBinLimits); 
  free(del2);   
  free(midKBin);
  free(meanKBin);
  free(binnedPk);
  free(modesPerBin);
  free(linearErrors);
  
  return 0;
}


int free_sdltInterp(){
  free(sdltk);
  free(sdltPk);
  free(sdlt2d);
  return 0;
}


int free_linear(){
    free(lineark);
    free(linearPk);
    free(linear2d);
    return 0;
}


int freeClipped(){  
  free(PkCube);
  free(Corrfn);
  free(suppressedCorrfn);
  free(distortedCorrfn);
  free(clippedPk);
  return 0;
}


int free2dPk(){
    for(j=0; j<n0*n1*n2; j++) free(TwoDpkArray[j]);
    free(TwoDpkArray);
    
    for(j=0; j<(kBinNumb-1); j++){
        free(zSpaceBinnedPk[j]);
        free(zSpacemodesPerBin[j]);
        free(mean_perpk[j]);
        free(mean_losk[j]);
    }
    
    free(zSpaceBinnedPk);
    free(zSpacemodesPerBin);
    free(mean_losk);
    free(mean_perpk);
    
    return 0;
}


int freeConvolutionMemory(){
    free(kVals);
    free(interpolatedPk);
    
    free(midKmodesperbin);
    free(midKBinNR);
    free(WindowFuncNR);
    free(WindowFunc2d);
    
    free(ConvolvedPk);
    free(ConvolvedPk2d);
    
    return 0;
}
