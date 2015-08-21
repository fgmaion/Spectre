int free_HOD(){
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
  free(Acceptanceflag);                                                                                           
  free(polarAngle);
  free(rDist);
  free(xCoor);
  free(yCoor);
  free(zCoor);
  
  return 0;
}


int free_grid(){
  free(overdensity);
  
  return 0;
}


int free_pkRegression(){
  fftw_free(H_k);
  
  fftw_destroy_plan(p);
    
  for(j=0; j<n0*n1*n2; j++)  free(polar_pk[j]);
    
  free(polar_pk);
    
  // binned Pk arrays.                                                                                                   
  free(logk_limits); 
  free(mean_modk);
  free(binnedPk);
  free(modes_perbin);
  
  free(kMonopole);
  free(kQuadrupole);
  free(kHexadecapole);
  
  return 0;
}


int free2dPk(){
    for(j=0; j<n0*n1*n2; j++){
        free(polar_pk[j]);
        // free(TwoDpkArray[j]);
    }
    
    free(polar_pk);
    // free(TwoDpkArray);
    
    /*
    free(muBinLimits);
    
    for(j=0; j<muBinNumb-1; j++)  free(polar_modesPerBin[j]);
    free(polar_modesPerBin);

    for(j=0; j<muBinNumb-1; j++)   free(mean_mu[j]);
    free(mean_mu);
    
    for(j=0; j<muBinNumb-1; j++)  free(mean_modk[j]);
    free(mean_modk);
    
    for(j=0; j<muBinNumb-1; j++)  free(polar2DBinnedPk[j]);
    free(polar2DBinnedPk);
    */   
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


int free_rand(){
    free(rand_ra); 
    free(rand_dec); 
    free(rand_chi);
    
    free(rand_x); 
    free(rand_y); 
    free(rand_z);

    return 0;
}
