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

	return Interim;
}


int convolve3DInputPk(){
    for(kk=0; kk<n0-2*zwfKernelsize; kk++){
	    for(jj=0; jj<n1-2*ywfKernelsize; jj++){
		    for(ii=0; ii<n2-2*xwfKernelsize; ii++){
				Index = kk*(n1-2*ywfKernelsize)*(n2-2*xwfKernelsize) + jj*(n2-2*xwfKernelsize) + ii;
				
				convolvedPk3d[Index]  = ConvolveCell(ii + xwfKernelsize, jj + ywfKernelsize, kk + zwfKernelsize);
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

				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

				flattenedConvolvedPk3D[totalModes][0] = (double) kmodulus;
				flattenedConvolvedPk3D[totalModes][1] = (double) convolvedPk3d[qqIndex];

                TwoDpkArray[totalModes][0]            = fabs(k_x);
                
                TwoDpkArray[totalModes][1]            = pow(k_y*k_y + k_z*k_z, 0.5);
                
                TwoDpkArray[totalModes][2]            = (double) convolvedPk3d[qqIndex];

				totalModes     += 1;
			}
		}
	}

	return 0;
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. convolution.");

    MeasureAnisoWfKernel();

    prepAnisoConvolution();

    pt2Pk = &splintMatterPk;

    inputPK();

    setInputPk();

    setMeasuredWfKernel();
    
    printWfKernel();

    convolve3DInputPk(convolvedPk3d, inputPk);
    
    // AnisoICC();

    Theory2DPk();
    
    return 0;	
}
