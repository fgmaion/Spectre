int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. convolution.");

    prepAnisoConvolution();

    pt2Pk                       = &splintMatterPk;

    minKernelshift              = -0.5*(wfKernelsize-1);
    maxKernelshift              =  0.5*(wfKernelsize-1);

    inputPK();

    setInputPk();

    setMeasuredWfKernel();

    convolve3DInputPk(convolvedPk3d, inputPk);

    // AnisoICC();

    flattenConvolvedPk();

    sprintf(filepath, "%s/Data/Del2k/midK_Pk_Convolved_ICCUncorrected%s.dat", root_dir, surveyType);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %e \t %d \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
    fclose(output);
    
    printWfKernel();

    return 0;	
}


double ConvolveCell(int x, int y, int z){
	double Interim            = 0.0;
	
	ConvNorm = ConvolutionNorm(windowFunc3D);

	for(kkshift=minKernelshift; kkshift<maxKernelshift+1; kkshift++){
	    for(jjshift=minKernelshift; jjshift<maxKernelshift + 1; jjshift++){
	        for(iishift=minKernelshift; iishift<maxKernelshift + 1; iishift++){
	        
			    filterIndex  = (kkshift + maxKernelshift)*wfKernelsize*wfKernelsize + (jjshift + maxKernelshift)*wfKernelsize + (iishift + maxKernelshift);
			    PkIndex      = (z + kkshift)*n1*n2 + (y + jjshift)*n2 + (x + iishift);

			    Interim     += pow(ConvNorm, -1.)*windowFunc3D[filterIndex]*inputPk[PkIndex];
            }
	    }
	}

    Interim /= wfKernelsize*wfKernelsize*wfKernelsize;

	return Interim;
}


int convolve3DInputPk(){
    for(kk=0; kk<n0-2*wfKernelsize; kk++){
	    for(jj=0; jj<n1-2*wfKernelsize; jj++){
		    for(ii=0; ii<n2-2*wfKernelsize; ii++){
				Index = kk*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + jj*(n2-2*wfKernelsize) + ii;
				
				convolvedPk3d[Index]  = ConvolveCell(ii + wfKernelsize, jj + wfKernelsize, kk + wfKernelsize);
			}
		}
	}

	return 0;
}


int flattenConvolvedPk(){
	int qqIndex = 0;
	int kkIndex = 0;

	int totalModes = 0;

	for(k=n0/2; k<n0-wfKernelsize; k++){
		for(j=n1/2; j<n1-wfKernelsize; j++){
			for(i=n2/2; i<n2-wfKernelsize; i++){
				qqIndex 		= (k-wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + (j-wfKernelsize)*(n2-2*wfKernelsize) + (i-wfKernelsize);

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
