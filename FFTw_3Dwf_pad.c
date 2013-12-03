int padwfPkCalc(int sidepad){
    printf("\n\nBeginning padded window function calculation.");
    // Padfactor must be be odd.
    padfactor = 2*sidepad + 1;
    
    prepFFTw(padfactor*n0, padfactor*n1, padfactor*n2);
    
    assignbinninginterval(kbinInterval/(float) padfactor);
    prepFFTbinning();
    
    padcellsNumber = padfactor*n0*padfactor*n1*padfactor*n2;
     
    for(j=0; j<padcellsNumber; j++) in[j][0] = 0.0;
    for(j=0; j<padcellsNumber; j++) in[j][1] = 0.0; 
    
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index    =   k*n1*n2 + j*n2 + i;
                padIndex =  (sidepad*n0 + k)*padfactor*n1*padfactor*n2 + (sidepad*n1 + j)*padfactor*n2 + (sidepad*n2 + i);
            
                in[padIndex][0] = (double) pow(padfactor, 3.)*(TotalVolume/TotalSurveyedVolume)*booldensity[Index]*FKPweights[Index];
            }
        }
    } 
    
    printf("\nPerforming FFT for padded window fn.");
    fftw_execute(p);
    
    printf("\nFFT for padd window fn. complete.");


    for(k=0; k<n0*padfactor; k++){
      for(j=0; j<n1*padfactor; j++){
        for(i=0; i<n2*padfactor; i++){
          k_x = (kIntervalx/padfactor)*i;
          k_y = (kIntervaly/padfactor)*j;
          k_z = (kIntervalz/padfactor)*k;

          if(k_x>NyquistWaveNumber)  k_x    -= (padfactor*n2)*(kIntervalx/padfactor);
          if(k_y>NyquistWaveNumber)  k_y    -= (padfactor*n1)*(kIntervaly/padfactor);
          if(k_z>NyquistWaveNumber)  k_z    -= (padfactor*n0)*(kIntervalz/padfactor);

          kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

          Index                              = k*padfactor*n1*padfactor*n2 + j*padfactor*n2 + i;

          H_kReal                            = pow(padcellsNumber, -1.0)*out[Index][0];
          H_kImag                            = pow(padcellsNumber, -1.0)*out[Index][1];

          PkArray[Index][0]                  = pow(kSq, 0.5);
          PkArray[Index][1]                  = pow(H_kReal, 2.) + pow(H_kImag, 2.);

        }
      }
    }

    printf("\nBinning padded W2.");
    PkBinningCalc(padcellsNumber, PkArray);
    
    sprintf(filepath, "%s/Data/WindowfuncSpherical/midK_pad%dW2k_%s.dat", root_dir, padfactor, surveyType);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++)  fprintf(output, "%f \t %f\n", meanKBin[j], binnedPk[j]);

    printf("\npadded Window function P(k) calculation complete.");
    
    fclose(output);
    
    return 0;
}
