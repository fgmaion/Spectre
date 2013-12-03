float anisoGauss(float x, float y, float z){
    float xsig = 1.0;
    float ysig = 0.5;
    float zsig = 0.2;
    
    return expf(-1.0*(pow(xsig*x, 2.) + pow(ysig*y, 2.) + pow(zsig*z, 2.)));
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. calculation.");

	prepAnisoConvolution();

	pt2Pk                       = &Analytic2powerlaw;
    pt2AnisoWf                  = &anisoGauss;

	SetWfKernel();

	inputLinearPk();

	setInputPk();

	wf3Dnorm  					= filter3Dnorm();

	printf("\nConvolved cell:  %f", ConvolveCell(n2/2, n1/2, n0/2));

	// convolve3DInputPk(convolvedPk3d, inputPk);

	// PkBinningCalc((n0-2*wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize), flattenedConvolvedPk3D);

	// sprintf(filepath, "%s/Data/Pk/midK_Pk_ConvolvedAnisoGauss.dat", root_dir);
    
    // output = fopen(filepath, "w");
    // for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
    // fclose(output);

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

	int   qIndex;
	int   kkIndex;

	// -q_z to q_z
	for(k=0; k<wfKernelsize; k++){
	  // -q_y to q_y
	  for(j=0; j<wfKernelsize; j++){
	  	// -q_x to q_x 
	    for(i=0; i<wfKernelsize; i++){
			qIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;

			// k indexing. 
			i       -= (wfKernelsize-1)/2;
			j       -= (wfKernelsize-1)/2;
			k       -= (wfKernelsize-1)/2;

			// kkIndex  = (z + 0.*k)*n1*n2 + (y + j*0.)*n2 + (x + i*0.);

			Interim += windowFunc3D[qIndex]; // *inputPk[kIndex];
		}
	  }
	}

	return Interim;
}


int flatten3dConvolvedPk(){
	int qIndex;

	for(k=0; k<n0-2*wfKernelsize; k++){
		for(j=0; j<n1-2*wfKernelsize; j++){
			for(i=0; i<n2-2*wfKernelsize; i++){
				qIndex = k*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + j*(n2-2*wfKernelsize) + (i - 2*wfKernelsize);

				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

				flattenedConvolvedPk3D[qIndex][0] =  kmodulus;
				flattenedConvolvedPk3D[qIndex][1] =  convolvedPk3d[Index];
			}
		}
	}

	return 0;
}


int convolve3DInputPk(float convolvedPk[], float inputPk[]){
	for(k=0; k<n0-2*wfKernelsize; k++){
		for(j=0; j<n1-2*wfKernelsize; j++){
			for(i=0; i<n2-2*wfKernelsize; i++){
				Index = k*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + j*(n2-2*wfKernelsize) + i;

				convolvedPk[Index]  = ConvolveCell(i + wfKernelsize, j + wfKernelsize, k + wfKernelsize);

				// convolvedPk[Index] /= wf3Dnorm;
			}
		}
	}

	return 0;
}


int SetWfKernel(){
	for(k=0; k<wfKernelsize; k++){
		for(j=0; j<wfKernelsize; j++){
			for(i=0; i<wfKernelsize; i++){
				k_x			   		= kIntervalx*(i-wfKernelsize/2.);
				k_y 		   		= kIntervaly*(j-wfKernelsize/2.);
				k_z			   		= kIntervalz*(k-wfKernelsize/2.);

				Index 		   		= k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;

				windowFunc3D[Index] = (*pt2AnisoWf)(k_x, k_y, k_z);
			}
		}
	}

	return 0;
}


float filter3Dnorm(){
	int   qIndex;
	float Interim = 0.0;

	for(k=0; k<wfKernelsize; k++){
		for(j=0; j<wfKernelsize; j++){
			for(i=0; i<wfKernelsize; i++){
				qIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;
				Interim += windowFunc3D[qIndex];
			}
		}
	}

	return Interim;
}