float anisoGauss(float x, float y, float z){
    float xsig = 1.0;
    float ysig = 0.5;
    float zsig = 0.2;
    
    return expf(-1.0*(pow(xsig*x, 2.) + pow(ysig*y, 2.) + pow(zsig*z, 2.)));
}


int AnisoConvolution(){
    printf("\n\nImplementing anisotropic window fn. calculation.");

	prepAnisoConvolution();

	inputLinearPk();

	pt2Pk                       = &Analytic2powerlaw;
    pt2AnisoWf                  = &anisoGauss;

	SetWfKernel();

	setInputPk();
/*
	for(k=0; k<n0-wfKernelsize; k++){
		for(j=0; j<n1-wfKernelsize; j++){
			for(i=0; i<n2-wfKernelsize; i++){
				Index = k*n1*n2 + j*n2 + i;

				convolvedPk3d[Index] = ConvolveCell(inputPk, i, j, k);
			}
		}
	}
*/
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

                inputPk[Index]  = Analytic2powerlaw(kmodulus);
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


float ConvolveCell(float array[], int x, int y, int z){
	float Interim = 0.0;

	int   qIndex;

	Index = z*n2*n1 + y*n2 + x; 

	for(k=0; k<10; k++){
		for(j=0; j<10; j++){
			for(i=0; i<10; i++){
				qIndex   = (z - 5 + k)*n2*n1 + (y - 5 + j)*n2 + (x - 5 + i);
				Interim += array[qIndex];
			}
		}
	}

	return Interim;
}
