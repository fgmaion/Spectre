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

	convolve3DInputPk(convolvedPk3d, inputPk);

	flatten3dConvolvedPk();

	sprintf(filepath, "%s/Data/Pk/midK_Pk_ConvolvedAnisoGauss.dat", root_dir);
    
    output = fopen(filepath, "w");
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%g \t %g \t %g \t %d \t %g \n", midKBin[j], del2[j], binnedPk[j], modesPerBin[j]);
    fclose(output);

	return 0;	
}


int flatten3dConvolvedPk(){
	int qqIndex;
	int kkIndex;

	printf("\nkIntervals, x: %f, y: %f, z: %f", kIntervalx, kIntervaly, kIntervalz);

	int totalModes = 0;

	// k indexing! ie not k=wfKernelsize to blah, not 0 to blah. 
	for(k=n0/2; k<n0-wfKernelsize; k++){
		for(j=n1/2; j<n1-wfKernelsize; j++){
			for(i=n2/2; i<n2-wfKernelsize; i++){
				// qqIndex 		= (k - wfKernelsize)*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + (j - wfKernelsize)*(n2-2*wfKernelsize) + (i - wfKernelsize);

				kkIndex 		= k*n1*n2 + j*n2 + i;

				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

				flattenedConvolvedPk3D[totalModes][0] =  kmodulus;
				flattenedConvolvedPk3D[totalModes][1] =  convolvedPk3d[qqIndex];

				totalModes     += 1;
			}
		}
	}

	for(j=0; j<100; j++) printf("\n %f \t %f", flattenedConvolvedPk3D[j][0], flattenedConvolvedPk3D[j][1]);

	PkBinningCalc(totalModes, flattenedConvolvedPk3D);

	return 0;
}


int convolve3DInputPk(float convolvedPk[], float inputPk[]){
	int ii, jj, kk;

	for(kk=0; kk<n0-2*wfKernelsize; kk++){
	    for(jj=0; jj<n1-2*wfKernelsize; jj++){
		    for(ii=0; ii<n2-2*wfKernelsize; ii++){
				Index = kk*(n1-2*wfKernelsize)*(n2-2*wfKernelsize) + jj*(n2-2*wfKernelsize) + ii;

				convolvedPk3d[Index]  = ConvolveCell(ii + wfKernelsize, jj + wfKernelsize, kk + wfKernelsize);

				convolvedPk3d[Index] /= wf3Dnorm;
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
	  for(j=0; j<wfKernelsize; j++){
	  	// -q_x to q_x 
	    for(i=0; i<wfKernelsize; i++){
			aaIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;

			// k indexing. 
			ishift  += 1;
			jshift  += 1;
			kshift  += 1; 

			bbIndex  = (z + kshift)*n1*n2 + (y + jshift)*n2 + (x + ishift);

			Interim += windowFunc3D[aaIndex]*inputPk[bbIndex];
		}
	  }
	}

	return Interim;
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