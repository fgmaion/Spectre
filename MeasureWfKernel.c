int MeasureAnisoWfKernel(){
    // Delivers a kernel of W^2(q) for convolution. 
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");

    // Allocates memory for (the averaging of) the measured window fn. kernel.
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
    
    printf("\nAssigning window function kernel.");
    
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
                
                // The kernel only extends to kx, ky, kz limited by e.g 0.037699, can ignore modes which have a modulus much greater than this, (safety factor of 1.1) 
                if((fabs(k_x) < xmaxKernelshift*kIntervalx) && (fabs(k_y) < ymaxKernelshift*kIntervaly) && (fabs(k_z) < zmaxKernelshift*kIntervalz)){
                    
                    kk = (int) ceil(k_z/kIntervalz) + zmaxKernelshift;
                    
                    jj = (int) ceil(k_y/kIntervaly) + ymaxKernelshift;
                    
                    ii = (int) ceil(k_x/kIntervalx) + xmaxKernelshift;
                    
                    Index = kk*ywfKernelsize*xwfKernelsize + jj*xwfKernelsize + ii; 
                 
                   /*
                    for(kkshift= zminKernelshift; kkshift<zmaxKernelshift + 1; kkshift++){
                        for(jjshift= yminKernelshift; jjshift<ymaxKernelshift + 1; jjshift++){ 
                            for(iishift= xminKernelshift; iishift<xmaxKernelshift + 1; iishift++){
                   */         
                                
                    // Index    = (kkshift + zmaxKernelshift)*ywfKernelsize*xwfKernelsize + (jjshift + ymaxKernelshift)*xwfKernelsize + (iishift + xmaxKernelshift); 
                                
                                /* if(    ((iishift - 0.1)*kIntervalx < k_x) && (k_x<(iishift + 0.1)*kIntervalx) 
                                    && ((jjshift - 0.1)*kIntervaly < k_y) && (k_y<(jjshift + 0.1)*kIntervaly) 
                                    && ((kkshift - 0.1)*kIntervalz < k_z) && (k_z<(kkshift + 0.1)*kIntervalz)){
                                */     
                                            AnisoWfKernel[Index]                      += pow(H_kReal, 2.) + pow(H_kImag, 2.);
                                            AnisoWfKernel_ModeNumb[Index]             += 1;
                                            
                                            if(AnisoWfKernel[Index] > pow(10., -9.)){
                                                wfKernel_minAmpIndices[largeAmpIndices][0] = ii - xmaxKernelshift;
                                                wfKernel_minAmpIndices[largeAmpIndices][1] = jj - ymaxKernelshift;
                                                wfKernel_minAmpIndices[largeAmpIndices][2] = kk - zmaxKernelshift;
                                                
                                                wfKernel_minAmpIndices[largeAmpIndices][3] = Index;
                 
                                                largeAmpIndices                           += 1;  
                                            } 
                                            
                                            
                                // }
                        /*
                            }                        
                        }
                    }
                    */
                }        
            }
        }                
    }
    
    /*
    for(kk=0; kk<zwfKernelsize; kk++){
        for(jj=0; jj<ywfKernelsize; jj++){
            for(ii=0; ii<xwfKernelsize; ii++){
                Index = kk*ywfKernelsize*xwfKernelsize + jj*xwfKernelsize + ii;   
            
                if(AnisoWfKernel_ModeNumb[Index] != 0){
                    AnisoWfKernel[Index] /= (double) AnisoWfKernel_ModeNumb[Index];
                }
            }
        }
    }
    */

    printf("\n\nUpper Wf kernel x limit: %g", xmaxKernelshift*kIntervalx);
    printf("\nUpper Wf kernel y limit: %g",   ymaxKernelshift*kIntervaly);
    printf("\nUpper Wf kernel z limit: %g\n\n",   zmaxKernelshift*kIntervalz);
    
    printf("\nPercentage of required cells in wk Kernel: %f", (float) 100.*largeAmpIndices/(xwfKernelsize*ywfKernelsize*zwfKernelsize));

    return 0;
}
