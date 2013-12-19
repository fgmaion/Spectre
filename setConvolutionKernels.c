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

                inputPk[Index]  = (double) (*pt2Pk)(kmodulus);
			}
		}
	}

	return 0;
}


int printWfKernel(){
    printf("\nApplied window fn. normalisation(convolution):  %e", ConvNorm);
    printf("\nApplied window fn. normalisation(zero point) :  %e", ZeroPointNorm());
    
    printf("\n\nFilter measured as: \n\n");
 
    for(i=0; i<wfKernelsize*wfKernelsize*wfKernelsize; i=i+wfKernelsize){
        for(j=0; j<wfKernelsize; j++){
            printf("%e \t", windowFunc3D[i + j]);
        }
        
        printf("\n");
    }
 
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

	return 0;
}


double ConvolutionNorm(double array[]){
	int   qIndex;
	
	double Interim = 0.0;

	for(k=0; k<wfKernelsize; k++){
		for(j=0; j<wfKernelsize; j++){
			for(i=0; i<wfKernelsize; i++){
				qIndex   = k*wfKernelsize*wfKernelsize + j*wfKernelsize + i;
				Interim += array[qIndex];
			}
		}
	}
	
	Interim /= wfKernelsize*wfKernelsize*wfKernelsize;

	return Interim;
}


double ZeroPointNorm(){
    return windowFunc3D[(wfKernelsize-1)/2*wfKernelsize*wfKernelsize + wfKernelsize*(wfKernelsize-1)/2 + (wfKernelsize-1)/2];
}


int setMeasuredWfKernel(){
    for(j=0; j<wfKernelsize*wfKernelsize*wfKernelsize; j++) windowFunc3D[j] = AnisoWfKernel[j];	
    
    return 0;
}


double AnalyticSpherical(double k){
    float y = k*250.;

    return 3.*pow(y, -3.)*(sin(y) - y*cos(y));
}


double AnalyticGaussian(double q){
    float r0 = 10.;

    return exp(-q*q*r0*r0);
}


double Analytic2powerlaw(double k){
    float p0 = 150.;
    float  n =   2.;

    return p0*pow(k, n);
}


double anisokGauss(double x, double y, double z){
    double xsig =   80.0;
    double ysig =  100.0;
    double zsig =  120.0;
    
    return exp(-1.0*(pow(xsig*x, 2.) + pow(ysig*y, 2.) + pow(zsig*z, 2.)));
}
