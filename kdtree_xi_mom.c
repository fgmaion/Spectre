int xiMonopole(double* xi, double* mu, double* Monopole){
  for(i=0; i<nlogbins; i++){
    Monopole[i] = 0.0;
    logrbins[i] = pow(10., zerolog + (i + 0.5)*logbinsz);
    
    for(j=0; j<nlinbins; j++){
      // second equality tests for NaNs in mu[i][j] (due to zero DD[i][j] pairs). NaN is the only number not equal to itself, not respected by all compilers. 
      // if((gr[i][j] <10) || (mu[i][j] != mu[i][j])){
      //    quality_flag = 0;
      // }
      
      // mu[i][j]     = zerolin + (j+0.5)*linbinsz;

      Index           = j + nlinbins*i;
      
      Monopole[i]    += xi[Index]; //*LegendrePolynomials(mu[i][j], 0)*linbinsz;
    }
  }    
    
  return 0;
}

int assignMemory_xi(){
  point_rands = (Particle *) realloc(point_rands, rand_number*sizeof(Particle));

  nlogbins    = (int) ceil((maxlog - zerolog)/logbinsz);
  nlinbins    = (int) ceil((maxlin - zerolin)/linbinsz);

  rr          = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_0        = calloc(nlogbins*nlinbins, sizeof(double*)); // window. 
  rr_2        = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_4        = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_6        = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_8        = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_10       = calloc(nlogbins*nlinbins, sizeof(double*));
    
  rr_meanr    = calloc(nlogbins*nlinbins, sizeof(double*));
  rr_meanmu   = calloc(nlogbins*nlinbins, sizeof(double*));
    
  xi0         = calloc(nlogbins, sizeof(*xi0));
  xi2         = calloc(nlogbins, sizeof(*xi2));
  xi4         = calloc(nlogbins, sizeof(*xi4));
  xi6         = calloc(nlogbins, sizeof(*xi6));
  xi8         = calloc(nlogbins, sizeof(*xi8));
  xi10        = calloc(nlogbins, sizeof(*xi10));

  logrbins    = calloc(nlogbins, sizeof(*logrbins));
  pairsperbin = calloc(nlogbins, sizeof(*pairsperbin));
              
  return 0;
}

int randWindow_pairCount(){    
    grow_randTree();
    
    printf("\n\nCounting RR pairs.");
    
    // bruteforce_nonodes(rr_0, rr_2, rr_4, rr_6, rr_8, rr_10, rr_meanr, rr_meanmu, point_rands, point_rands, rand_number, rand_number, 1);
    findSuitableNodePairs_bruteforcePairCount(rr_0, rr_2, rr_4, rr_6, rr_8, rr_10, rr_meanr, rr_meanmu, randTree, randTree, 1);

    postprocess_pairs(rr_0, rr_2, rr_4, rr_6, rr_8, rr_10, rr_meanr, rr_meanmu, randTree, randTree, 1);
    
    xiMonopole(rr_0,  rr_meanmu,  xi0);  // monopole calc. of *weighted* pairs.
    xiMonopole(rr_2,  rr_meanmu,  xi2);  // PREVIOUSLY gg_meanmu.
    xiMonopole(rr_4,  rr_meanmu,  xi4);
    xiMonopole(rr_6,  rr_meanmu,  xi6);
    xiMonopole(rr_8,  rr_meanmu,  xi8);  // masked RSD work.
    
    print_xiMultipoles();
    
    return 0;
}

int print_xiMultipoles(){
  sprintf(filepath, "%s/Qmultipoles/%s.dat", outputdir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++)  fprintf(output, "%.4le \t %.4le \t %.4le \t %.4le \t %.4le \t %.4le \n", logrbins[j], xi0[j], xi2[j], xi4[j], xi6[j], xi8[j]);
  
  fclose(output);

  return 0;
}
