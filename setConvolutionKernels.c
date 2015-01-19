int setInputPk(){
    double checkMax = 0.0;

	for(k=0; k<(n0+1); k++){
		for(j=0; j<(n1+1); j++){
			for(i=0; i<(n2+1); i++){
				k_x   	 		= kIntervalx*(i - n2/2.);
                k_y   	 		= kIntervaly*(j - n1/2.);
                k_z   	 		= kIntervalz*(k - n0/2.);

                Index 	 		= k*(n1+1)*(n2+1) + j*(n2+1) + i;

                kSq      		= pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                kmodulus 		= pow(kSq, 0.5);

                inputPk[Index]  = (double) (*pt2Pk)(kmodulus);
                
                if(inputPk[Index] >= checkMax){
                    checkMax    = inputPk[Index];
                }
                
                /*
                TwoDpkArray[totalModes][0] = fabs(k_x);
                
                TwoDpkArray[totalModes][1] = pow(k_y*k_y + k_z*k_z, 0.5);
                
                TwoDpkArray[totalModes][2] = (double) inputPk[Index];
                */
                
			}
		}
	}
	
	printf("\nMaximum input P(k) is %e", checkMax);
	
	return 0;
}

/*
int SetWfKernel_minAmp(){
    double amplitude = 0.0;
    
    largeAmpIndices  =   0;

    for(k=0; k<zwfKernelsize; k++){
		for(j=0; j<ywfKernelsize; j++){
			for(i=0; i<xwfKernelsize; i++){
			    k_x			   		= kIntervalx*(xminKernelshift + i);
			    k_y 		   		= kIntervaly*(yminKernelshift + j);
			    k_z			   		= kIntervalz*(zminKernelshift + k);

			    Index 		   		= k*ywfKernelsize*xwfKernelsize + j*xwfKernelsize + i;
			    
			    windowFunc3D[Index] = (*pt2AnisoWf)(k_x, k_y, k_z);

			    amplitude           = (*pt2AnisoWf)(k_x, k_y, k_z);
			    
			    printf("\n %e \t %e \t %e \t %e", k_x, k_y, k_z, amplitude);
			
			    if(amplitude > pow(10., -9.)){
                    wfKernel_minAmpIndices[largeAmpIndices][0] = xminKernelshift + i;
                    wfKernel_minAmpIndices[largeAmpIndices][1] = yminKernelshift + j;
                    wfKernel_minAmpIndices[largeAmpIndices][2] = zminKernelshift + k;
                                                   
                    wfKernel_minAmpIndices[largeAmpIndices][3] = Index;
                 
                    largeAmpIndices                           += 1;  
                } 
			}
		}
	}      
                                            
    return 0;
}
*/

double splintHODpk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    if(k>70.)          return 2.51*pow((k/70.), -1.8)*pow(linearBias/1.495903, 2.);
    
    else if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 0.96)*pow(linearBias/1.495903, 2.); 

    else{
        double Interim;
    
        splint(sdltk, sdltPk, sdlt2d, 469, k, &Interim);
    
        return Interim*pow(linearBias/1.495903, 2.);
    }
}

/*
int  printHODMatterPk(){
    double waveNumber    = 0.0;
    
    sprintf(filepath, "%s/Data/Del2k/spline_HODmatterPk.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<100000; j++){
        waveNumber = pow(10., -5.)*j;

        fprintf(output, "%le \t %le \n", waveNumber, splintHODpk(waveNumber));
    }
    
    fclose(output);

    return 0;
}
*/

double splintLinearPk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    
    // if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 1.00)*pow(linearBias/1.495903, 2.); 
    if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 0.96)*pow(linearBias/1.495903, 2.); 

    else{
        splint(lineark, linearPk, linear2d, 469, k, &Interim);
    
        return Interim*pow(linearBias/1.495903, 2.);
    }
}

/*
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
*/

double minAmp_ConvolutionNorm(double array[]){
    double Interim = 0.0;

	for(k=0; k<largeAmpIndices; k++){
		Interim += array[k];
	}

    printf("\nConvolution norm. is %e ", Interim);

    return Interim;
}
/*


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

*/

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
