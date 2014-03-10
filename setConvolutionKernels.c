int setInputPk(){
    int totalModes = 0;

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
                
                /*
                TwoDpkArray[totalModes][0] = fabs(k_x);
                
                TwoDpkArray[totalModes][1] = pow(k_y*k_y + k_z*k_z, 0.5);
                
                TwoDpkArray[totalModes][2] = (double) inputPk[Index];
                */
                
                totalModes     += 1;
			}
		}
	}
	
	// 2D Power spectrum.
	// sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/InputTheory2Dpk_%s.dat", root_dir, surveyType);
    // Cartesian2Dpk(filepath);

	return 0;
}


int printWfKernel(){
    ConvNorm = ConvolutionNorm(windowFunc3D);

    printf("\nApplied window fn. normalisation(convolution):  %e", ConvNorm);
    printf("\nApplied window fn. normalisation(zero point) :  %e", ZeroPointNorm());
    
    printf("\n\nFilter measured as: \n\n");
 
    for(i=0; i<xwfKernelsize*ywfKernelsize*zwfKernelsize; i=i+zwfKernelsize){
        for(j=0; j<zwfKernelsize; j++){
            printf("%e \t", windowFunc3D[i + j]);
        }
        
        printf("\n");
    }
 
    return 0;
}


int SetWfKernel(){
	for(k=0; k<zwfKernelsize; k++){
		for(j=0; j<ywfKernelsize; j++){
			for(i=0; i<xwfKernelsize; i++){
			  k_x			   		= kIntervalx*(i-(xwfKernelsize-1)/2.);
			  k_y 		   		    = kIntervaly*(j-(ywfKernelsize-1)/2.);
			  k_z			   		= kIntervalz*(k-(zwfKernelsize-1.)/2.);

			  Index 		   		= k*ywfKernelsize*xwfKernelsize + j*xwfKernelsize + i;

			  windowFunc3D[Index]   = (*pt2AnisoWf)(k_x, k_y, k_z);
			}
		}
	}

	return 0;
}


double splintMatterPk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 0.96)*pow(linearBias/1.495903, 2.); 

    else{
        float Interim;
    
        splint(sdltk, sdltPk, sdlt2d, 469, (float) k, &Interim);
    
        return (double) Interim*pow(linearBias/1.495903, 2.);
    }
}


double ConvolutionNorm(double array[]){
	int   qIndex;
	
	double Interim = 0.0;

	for(k=0; k<zwfKernelsize; k++){
		for(j=0; j<ywfKernelsize; j++){
			for(i=0; i<xwfKernelsize; i++){
				qIndex   = k*ywfKernelsize*xwfKernelsize + j*xwfKernelsize + i;
				
				Interim += array[qIndex];
			}
		}
	}
	
	Interim /= xwfKernelsize*ywfKernelsize*zwfKernelsize;

	return Interim;
}


double ZeroPointNorm(){
    return windowFunc3D[(zwfKernelsize-1)/2*ywfKernelsize*xwfKernelsize + xwfKernelsize*(ywfKernelsize-1)/2 + (xwfKernelsize-1)/2];
}


int setMeasuredWfKernel(){
    for(j=0; j<xwfKernelsize*ywfKernelsize*zwfKernelsize; j++) windowFunc3D[j] = AnisoWfKernel[j];	
    
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
