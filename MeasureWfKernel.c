int MeasureAnisoWfKernel(){
    // Delivers a kernel of W^2(q) for convolution. 
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = Cell_AppliedWindowFn[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nFFT complete.");

    // Allocates memory for (the averaging of) the measured window fn. kernel.
    assignAnisoWfKernel();
    
    // Largest fundamental mode.
    MaxkInterval = fmax(fmax(kIntervalx, kIntervaly), kIntervalz);
    
    minKernelshift              = -0.5*(wfKernelsize-1);
    maxKernelshift              =  0.5*(wfKernelsize-1);
    
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
                
                // The kernel only extends to kx, ky, kz limited by e.g 0.037699, can ignore modes which have a modulus much greater than this. 
                if(fabs(kmodulus) < pow(3., 0.5)*(wfKernelsize-1.)*MaxkInterval){
                
                    for(kkshift= minKernelshift; kkshift<maxKernelshift + 1; kkshift++){
                        for(jjshift= minKernelshift; jjshift<maxKernelshift + 1; jjshift++){ 
                            for(iishift= minKernelshift; iishift<maxKernelshift + 1; iishift++){
                            
                                Index    = (kkshift + maxKernelshift)*wfKernelsize*wfKernelsize + (jjshift + maxKernelshift)*wfKernelsize + (iishift + maxKernelshift); 
                                
                                if(    ((iishift-0.1)*kIntervalx < k_x) && (k_x<(iishift+0.1)*kIntervalx) 
                                    && ((jjshift-0.1)*kIntervaly < k_y) && (k_y<(jjshift+0.1)*kIntervaly) 
                                    && ((kkshift-0.1)*kIntervalz < k_z) && (k_z<(kkshift+0.1)*kIntervalz)){
                                            AnisoWfKernel[Index]           += pow(H_kReal, 2.) + pow(H_kImag, 2.);
                                            AnisoWfKernel_ModeNumb[Index]  += 1;
                                }
                            }                        
                        }
                    }
                    
                }
                    
                    
            }
        }                
    }
    
    for(kk=0; kk<wfKernelsize; kk++){
        for(jj=0; jj<wfKernelsize; jj++){
            for(ii=0; ii<wfKernelsize; ii++){
                Index = kk*wfKernelsize*wfKernelsize + jj*wfKernelsize + ii;   
            
                if(AnisoWfKernel_ModeNumb[Index] != 0){
                    AnisoWfKernel[Index] /= (double) AnisoWfKernel_ModeNumb[Index];
                }
            }
        }
    }

    printf("\n\nUpper Wf kenel x limit: %g", 0.5*(wfKernelsize-1)*kIntervalx);
    printf("\nUpper Wf kenel y limit: %g", 0.5*(wfKernelsize-1)*kIntervaly);
    printf("\nUpper Wf kenel z limit: %g", 0.5*(wfKernelsize-1)*kIntervalz);

    return 0;
}
