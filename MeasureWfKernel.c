int MeasureAnisoWfKernel(){
    // Delivers a kernel of W^2(q) for convolution. 
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");
    
    printf("\nAssigning window function kernel.");
    
    double W2k = 0.0;
    
    largeAmpIndices = 0;
    
    int W2x, W2y, W2z;
    
    // sprintf(filepath, "%s/Data/windowfuncHist.data", root_dir);
    
    // output = fopen(filepath, "w");
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index                              = k*n1*n2 + j*n2 + i;

                H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
                H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];
                
                W2k                                = pow(H_kReal, 2.) + pow(H_kImag, 2.);
                
                if(W2k > wfKernel_minAmp){
                    largeAmpIndices               += 1;  
                }             
            }
        }                
    }

    printf("\n\nRetained window fn. elements: %d", largeAmpIndices);
    
    prep_wfKernelminAmp();
    
    int modeCount = 0;
        
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus                           = pow(kSq, 0.5);

                H_kReal                            = pow(n0*n1*n2, -1.0)*out[Index][0];
                H_kImag                            = pow(n0*n1*n2, -1.0)*out[Index][1];
                
                W2k                                = pow(H_kReal, 2.) + pow(H_kImag, 2.);
                
                if((i==0) && (j==0) && (k==0)){
                    Wfzeropoint = W2k;
                }
                                        
                W2x                                = i;
                W2y                                = j;
                W2z                                = k;
                
                if(i>n2/2)  W2x                   -= n2;
                if(j>n1/2)  W2y                   -= n1;
                if(k>n0/2)  W2z                   -= n0;
                
                if(W2k > wfKernel_minAmp){
                    wfKernel_minAmpIndices[modeCount][0] = W2x;
                    wfKernel_minAmpIndices[modeCount][1] = W2y;
                    wfKernel_minAmpIndices[modeCount][2] = W2z;
                                                                
                    windowFunc3D[modeCount]              = W2k;
                                                                  
                    modeCount                           += 1; 
                    
                    printf("\n%d \t %d \t %d \t %e", W2x, W2y, W2z, W2k);
                    
                    if(W2x>max_xShift) max_xShift = W2x;
                    if(W2x<min_xShift) min_xShift = W2x;
    
                    if(W2y>max_yShift) max_yShift = W2y;
                    if(W2y<min_yShift) min_yShift = W2y;
                    
                    if(W2z>max_zShift) max_zShift = W2z;
                    if(W2z<min_zShift) min_zShift = W2z;
                }             
            }
        }                
    }
    
    printf("\nx shift range: %d, %d", min_xShift, max_xShift);
    printf("\ny shift range: %d, %d", min_yShift, max_yShift);
    printf("\nz shift range: %d, %d", min_zShift, max_zShift);
        
    return 0;
}
