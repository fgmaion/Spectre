int EstimateAnisoWfKernel(){
    printf("\n\nBeginning measurement of window function.");
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = (double) Cell_AppliedWindowFn[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = (double) 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);

    printWindowfuncSlices();
    
    printf("\nFFT complete.");
    
    assignAnisoWfKernel();

    float iishift, kkshift, jjshift;
    int   ii, jj, kk;
    
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
                
                
                if(fabsf(kmodulus) < 0.6*(wfKernelsize-1.)*fmaxf(fmaxf(kIntervalx, kIntervaly), kIntervalz)){
                
                    kkshift = -(wfKernelsize-1)/2;
                    
                    for(kk=0; kk<wfKernelsize; kk++){
                        
                        jjshift = -(wfKernelsize-1)/2;
                        
                        for(jj=0; jj<wfKernelsize; jj++){
                        
                            iishift = -(wfKernelsize-1)/2;
                            
                            for(ii=0; ii< wfKernelsize; ii++){
                                Index    = kk*wfKernelsize*wfKernelsize + jj*wfKernelsize + ii; 
                                
                                if(((iishift-0.05)*kIntervalx < k_x) && (k_x<(iishift+0.05)*kIntervalx)){
                                    if(((jjshift-0.05)*kIntervaly < k_y) && (k_y<(jjshift+0.05)*kIntervaly)){
                                        if(((kkshift-0.05)*kIntervalz < k_z) && (k_z<(kkshift+0.05)*kIntervalz)){
                                    
                                            AnisoWfKernel[Index]           += pow(H_kReal, 2.) + pow(H_kImag, 2.);
                                            AnisoWfKernel_ModeNumb[Index]  += 1;
                                        }
                                    }                        
                                }
                        
                                iishift += 1;
                            }
                    
                            jjshift += 1;
                        }
                    
                        kkshift += 1;
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
                    AnisoWfKernel[Index] /= (float) AnisoWfKernel_ModeNumb[Index];
                }
            }
        }
    }
    
    wf3Dnorm = filter3Dnorm(AnisoWfKernel);

    return 0;
}


float anisoGauss(float x, float y, float z){
    float xsig =  80.0;
    float ysig =  100.0;
    float zsig =  120.0;
    
    return expf(-1.0*(pow(xsig*x, 2.) + pow(ysig*y, 2.) + pow(zsig*z, 2.)));
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. calculation.");

	prepAnisoConvolution();

	pt2Pk                       = &splintMatterPk;
    // pt2AnisoWf               = &anisoGauss;

	setMeasuredWfKernel();

	inputPK();

	// inputLinearPk();

	setInputPk();

	convolve3DInputPk(convolvedPk3d, inputPk);

	// AnisoICC();

	ConvolvedPkZeroPoint     = ConvolveCell(n2/2, n1/2, n0/2);
    printf("\nConvolved P(vec k) zero point calculated to be: %e", ConvolvedPkZeroPoint);

    Index = 0.5*(wfKernelsize - 1)*wfKernelsize*wfKernelsize + 0.5*(wfKernelsize - 1)*wfKernelsize + 0.5*(wfKernelsize - 1);
    printf("\nWindow fn. zero point calculated to be:  %f", windowFunc3D[Index]);

	flatten3dConvolvedPk();

	sprintf(filepath, "%s/Data/Del2k/midK_Pk_Convolved%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
    fclose(output);

    return 0;	
}


int AnisoICC(){    
	// Integral constraint correction in the 3D anisotropic case. 
    int   blahIndex;
	int   bbIndex;

	int   ishift = -(wfKernelsize-1)/2;
	int   jshift = -(wfKernelsize-1)/2;
	int   kshift = -(wfKernelsize-1)/2;

	// -q_z to q_z
	for(k=0; k<wfKernelsize; k++){
	  // -q_y to q_y

	  jshift = -(wfKernelsize-1)/2;

	  for(j=0; j<wfKernelsize; j++){
	  	// -q_x to q_x 
	    
	    ishift = -(wfKernelsize-1)/2;
	    
	    for(i=0; i<wfKernelsize; i++){
			blahIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;
			bbIndex     = (n0/2 + kshift)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + (n1/2 + jshift)*(n2-2*wfKernelsize) + (n2/2 + ishift);

			convolvedPk3d[bbIndex] -= windowFunc3D[blahIndex];

			// k indexing.                                                                                      
            ishift  += 1;
        }
	  
	    jshift += 1;
	  }
	  
	  kshift += 1;
	}

    return 0;
}


int flatten3dConvolvedPk(){
	int qqIndex = 0;
	int kkIndex = 0;

	int totalModes = 0;

	for(k=n0/2; k<n0-wfKernelsize; k++){
		for(j=n1/2; j<n1-wfKernelsize; j++){
			for(i=n2/2; i<n2-wfKernelsize; i++){
				qqIndex 		= (k - wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + (j - wfKernelsize)*(n2-2*wfKernelsize) + (i - wfKernelsize);

				kkIndex 		= k*n1*n2 + j*n2 + i;

				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

				flattenedConvolvedPk3D[totalModes][0] =  (double) kmodulus;
				flattenedConvolvedPk3D[totalModes][1] =  (double) convolvedPk3d[qqIndex];

				totalModes     += 1;
			}
		}
	}

	PkBinningCalc(totalModes, flattenedConvolvedPk3D);

	return 0;
}


int convolve3DInputPk(float convolvedPk[], float inputPk[]){
	int ii, jj, kk;

	for(kk=0; kk<n0-2*wfKernelsize; kk++){
	    for(jj=0; jj<n1-2*wfKernelsize; jj++){
		    for(ii=0; ii<n2-2*wfKernelsize; ii++){
				Index = kk*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + jj*(n2-2*wfKernelsize) + ii;
				convolvedPk3d[Index]  = (1./wf3Dnorm)*ConvolveCell(ii + wfKernelsize, jj + wfKernelsize, kk + wfKernelsize);
			}
		}
	}

	return 0;
}


int setInputPk(){
	for(k=0; k<n0; k++){
		for(j=0; j<n1; j++){
			for(i=0; i<n2; i++){
				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                Index 	 		= k*n1*n2 + j*n2 + i;

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

                inputPk[Index]  = (*pt2Pk)(kmodulus);
			}
		}
	}

	return 0;
}


float ConvolveCell(int x, int y, int z){
	float Interim = 0.0;

	int   aaIndex;
	int   bbIndex;

	int   ishift = -(wfKernelsize-1)/2;
	int   jshift = -(wfKernelsize-1)/2;
	int   kshift = -(wfKernelsize-1)/2;

	// -q_z to q_z
	for(k=0; k<wfKernelsize; k++){
	  // -q_y to q_y

	  jshift = -(wfKernelsize-1)/2;

	  for(j=0; j<wfKernelsize; j++){
	  	// -q_x to q_x 
	    
	    ishift = -(wfKernelsize-1)/2;
	    
	    for(i=0; i<wfKernelsize; i++){
			aaIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;
			bbIndex  = (z + kshift)*n1*n2 + (y + jshift)*n2 + (x + ishift);

			Interim += windowFunc3D[aaIndex]*inputPk[bbIndex];

			// k indexing.                                                                                      
                        ishift  += 1;
                 }
	  
	    jshift += 1;
	  }
	  
	  kshift += 1;
	}

	return Interim;
}


int setMeasuredWfKernel(){
    for(j=0; j<wfKernelsize*wfKernelsize*wfKernelsize; j++) windowFunc3D[j] = AnisoWfKernel[j];	

return 0;
}


int SetWfKernel(){
	for(k=0; k<wfKernelsize; k++){
		for(j=0; j<wfKernelsize; j++){
			for(i=0; i<wfKernelsize; i++){
			  k_x			   		= kIntervalx*(i-(wfKernelsize-1)/2.);
			  k_y 		   		    = kIntervaly*(j-(wfKernelsize-1)/2.);
			  k_z			   		= kIntervalz*(k-(wfKernelsize-1.)/2.);

			  Index 		   		= k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;

			  windowFunc3D[Index]   = (*pt2AnisoWf)(k_x, k_y, k_z);
			}
		}
	}
	
	wf3Dnorm  					= filter3Dnorm(windowFunc3D);
	
	for(j=0; j<wfKernelsize*wfKernelsize*wfKernelsize; j++) windowFunc3D[j] /= wf3Dnorm;

	return 0;
}


float filter3Dnorm(float array[]){
	int   qIndex;
	float Interim = 0.0;

	for(k=0; k<wfKernelsize; k++){
		for(j=0; j<wfKernelsize; j++){
			for(i=0; i<wfKernelsize; i++){
				qIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;
				Interim += array[qIndex];
			}
		}
	}

	printf("\n3D filter normalisation:  %f", Interim);

	return Interim;
}