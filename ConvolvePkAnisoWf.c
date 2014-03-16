double ConvolveCell(int x, int y, int z){
	double Interim            = 0.0;
	
	ConvNorm = ConvolutionNorm(windowFunc3D);

	for(kkshift=zminKernelshift; kkshift<zmaxKernelshift+1; kkshift++){
	    for(jjshift=yminKernelshift; jjshift<ymaxKernelshift + 1; jjshift++){
	        for(iishift=xminKernelshift; iishift<xmaxKernelshift + 1; iishift++){
	        
			    filterIndex  = (kkshift + zmaxKernelshift)*ywfKernelsize*xwfKernelsize + (jjshift + ymaxKernelshift)*xwfKernelsize + (iishift + xmaxKernelshift);
			    PkIndex      = (z + kkshift)*n1*n2 + (y + jjshift)*n2 + (x + iishift);

			    Interim     += pow(ConvNorm, -1.)*windowFunc3D[filterIndex]*inputPk[PkIndex];
            }
	    }
	}

    Interim /= xwfKernelsize*ywfKernelsize*zwfKernelsize;

    // Looks perfect for full cube. 
	return Interim;
}


double minAmp_ConvolveCell(int x, int y, int z){
	double Interim            = 0.0;
    
    for(k=0; k<largeAmpIndices; k++){
    	PkIndex  = (z + wfKernel_minAmpIndices[k][2])*n1*n2 + (y + wfKernel_minAmpIndices[k][1])*n2 + (x + wfKernel_minAmpIndices[k][0]);
    
		Interim += pow(ConvNorm, -1.)*windowFunc3D[wfKernel_minAmpIndices[k][3]]*inputPk[PkIndex];
	}
    
    Interim /= largeAmpIndices;

    // Looks perfect for full cube. 
	return Interim;
}


int convolve3DInputPk(){
    // Now convolves solely the positive modes. 
    
    double kx, ky, kz, kSq;

    for(kk=n0/2 - zwfKernelsize; kk<n0-2*zwfKernelsize; kk++){
	    for(jj=n1/2 - ywfKernelsize; jj<n1-2*ywfKernelsize; jj++){
		    for(ii=n2/2 - xwfKernelsize; ii<n2-2*xwfKernelsize; ii++){
		        kx  = (ii -n2/2 + xwfKernelsize)*kIntervalx; 
		        
		        ky  = (jj -n1/2 + ywfKernelsize)*kIntervaly; 
		        
		        kz  = (kk -n0/2 + zwfKernelsize)*kIntervalz; 
		        
		        kSq = kx*kx + ky*ky + kz*kz; 
		    
                if((0.01*0.01 < kSq) && (kSq < 0.50*0.50)){
    				Index = kk*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize) + jj*(n2-2*xwfKernelsize) + ii;
				
	    			// convolvedPk3d[Index]  = ConvolveCell(ii + xwfKernelsize, jj + ywfKernelsize, kk + zwfKernelsize);
	    			
	    			convolvedPk3d[Index]     = minAmp_ConvolveCell(ii + xwfKernelsize, jj + ywfKernelsize, kk + zwfKernelsize);
	            }
			}
		}
	}

	return 0;
}


int Theory2DPk(){
	int qqIndex = 0;
	int kkIndex = 0;

	int totalModes = 0;

	for(k=n0/2; k<n0-zwfKernelsize; k++){
		for(j=n1/2; j<n1-ywfKernelsize; j++){
			for(i=n2/2; i<n2-xwfKernelsize; i++){
				qqIndex 		= (k-zwfKernelsize)*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize) + (j-ywfKernelsize)*(n2-2*xwfKernelsize) + (i-xwfKernelsize);
				
				Index 	 		= k*n1*n2 + j*n2 + i;

				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

                if((0.01*0.01 < kSq) && (kSq < 0.50*0.50)){
    				flattenedConvolvedPk3D[totalModes][0] = (double) kmodulus;
	    			flattenedConvolvedPk3D[totalModes][1] = (double) convolvedPk3d[qqIndex];

                    /*
                    TwoDpkArray[totalModes][0]            = fabs(k_x);
                    TwoDpkArray[totalModes][1]            = pow(k_y*k_y + k_z*k_z, 0.5);
                    TwoDpkArray[totalModes][2]            = (double) convolvedPk3d[qqIndex];
                    */
                
	    			totalModes     += 1;
	    		}
			}
		}
	}
	
	PkBinningCalc(totalModes, flattenedConvolvedPk3D);
	
	sprintf(filepath, "%s/Data/ConvolvedPk/IntegralConstraint_Corrected/%s_kbinInterval_%.2f.dat", root_dir, surveyType, kbinInterval);
	
    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %e \t %d \t %e \t %e\n", meanKBin[j], del2[j], binnedPk[j], modesPerBin[j], linearErrors[j], (*pt2Pk)(meanKBin[j]));

    fclose(output);
    
	return 0;
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. convolution.");

    MeasureAnisoWfKernel();

    prepAnisoConvolution();

    setInputPk();

    setMeasuredWfKernel();
    
    // printWfKernel();
    
    ConvNorm = minAmp_ConvolutionNorm(windowFunc3D);

    convolve3DInputPk(convolvedPk3d, inputPk);
    
    printf("\nConvolution norm.:  %e", minAmp_ConvolutionNorm(windowFunc3D));
    
    // minAmp_AnisoICC();
    
    Theory2DPk();
    
    return 0;	
}
