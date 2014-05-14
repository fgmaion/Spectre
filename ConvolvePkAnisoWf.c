double minAmp_ConvolveCell(int x, int y, int z){
	double Interim            = 0.0;
    
    for(k=0; k<largeAmpIndices; k++){
    	// PkIndex  = (z + wfKernel_minAmpIndices[k][2])*(n1+1)*(n2+1) + (y + wfKernel_minAmpIndices[k][1])*(n2+1) + (x + wfKernel_minAmpIndices[k][0]);
    
        k_x   	 		= kIntervalx*(x + wfKernel_minAmpIndices[k][0] - n2/2.);
        k_y   	 		= kIntervaly*(y + wfKernel_minAmpIndices[k][1] - n1/2.);
        k_z   	 		= kIntervalz*(z + wfKernel_minAmpIndices[k][2] - n0/2.);
    
        kSq             = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

        kmodulus        = pow(kSq, 0.5);
        
        mu                                 = k_x/kmodulus;
        if(kmodulus < 0.000001)       mu   = 0.0;  
	
	    Interim        += windowFunc3D[k]*(*pt2Pk)(kmodulus)*pow(1. + beta*mu*mu, 2.)/TotalVolume;
	}

	return Interim/ConvNorm;
}


int convolve3DInputPk(){
    // Now convolves solely the positive modes. 
    
    int    modeCount =   0;
    
    double kx, ky, kz, kSq;

    for(kk=0;     kk<(n0+1); kk++){
        printf("\n%d", kk);
    
	    for(jj=0; jj<(n1+1); jj++){
		    
		    // Only half the modes of P*(k) are independent, P*(-k) = P*(k). Choose kx >= 0.
		    for(ii=n2/2;      ii<(n2+1); ii++){ 	
		       	k_x   	 		= kIntervalx*(ii - n2/2.);
                k_y   	 		= kIntervaly*(jj - n1/2.);
                k_z   	 		= kIntervalz*(kk - n0/2.);

                Index 	 		= kk*(n1+1)*(n2+1) + jj*(n2+1) + ii;

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);
		       	
                if(kmodulus <= convolution_modkmax){
    			    convolvedPk3d[Index]                 = minAmp_ConvolveCell(ii, jj, kk);
	            }
			}
		}
	}
	
    	
	// Integral constraint correction. 
	for(k=0; k<largeAmpIndices; k++){
	    Index 	 		        = (n0/2 + wfKernel_minAmpIndices[k][2])*(n1+1)*(n2+1) + (n1/2 + wfKernel_minAmpIndices[k][1])*(n2+1) + (n2/2 + wfKernel_minAmpIndices[k][0]);
	
	    convolvedPk3d[Index]    = convolvedPk3d[Index] - (windowFunc3D[k]/Wfzeropoint)*ConvPkZeroPoint; 
	}
    	
	
	modeCount         = 0;
	polarPk_modeCount = 0;
	
	for(kk=0;     kk<(n0+1); kk++){    
	    for(jj=0; jj<(n1+1); jj++){
		    // Only half the modes of P*(k) are independent, P*(-k) = P*(k). Choose kx >= 0.
		    for(ii=n2/2;      ii<(n2+1); ii++){ 	
		       	k_x   	 		= kIntervalx*(ii - n2/2.);
                k_y   	 		= kIntervaly*(jj - n1/2.);
                k_z   	 		= kIntervalz*(kk - n0/2.);

                Index 	 		= kk*(n1+1)*(n2+1) + jj*(n2+1) + ii;

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);
                
                mu              = k_x/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0; 
		       	
                if(kmodulus <= convolution_modkmax){
    			    flattenedConvolvedPk3D[modeCount][0] = kmodulus;
	    			flattenedConvolvedPk3D[modeCount][1] = convolvedPk3d[Index];
	    			
	    			modeCount                           += 1;
	            }
	            
	            
	            if((kmodulus > 0.000001) && (kmodulus <= convolution_modkmax)){	            
		            // Issue with mu for a zeroth length vector being ill defined. 
		            polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		            polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		            polar2Dpk[polarPk_modeCount][2]    = convolvedPk3d[Index];
                    
                    polarPk_modeCount                 += 1;
	            }
			}
		}
	}
	
	observedQuadrupole(polarPk_modeCount);

	return 0;
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. convolution.");
    
    prepAnisoConvolution();
    
    MeasureAnisoWfKernel();
  
    ConvNorm                       = minAmp_ConvolutionNorm(windowFunc3D);
    
    ConvPkZeroPoint                = minAmp_ConvolveCell(n2/2, n1/2, n0/2);
    
    printf("\nConvolved P(vec k) zero point calculated to be: %e", ConvPkZeroPoint);
    printf("\nWindow fn. zero point calculated to be:  %e",        Wfzeropoint);
    
    convolve3DInputPk(convolvedPk3d, inputPk);
    
    // PkBinningCalc(modeCount, flattenedConvolvedPk3D);
        
    // sprintf(filepath, "%s/Data/Del2k/midK_Pk_ConvolvedHOD_pk_10e-8.dat", root_dir);

    // output = fopen(filepath, "w");
  
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %d \t %e\n", meanKBin[j], binnedPk[j], modesPerBin[j], splintHODpk(midKBin[j]));
  
    // fclose(output);  
    
    return 0;
}
