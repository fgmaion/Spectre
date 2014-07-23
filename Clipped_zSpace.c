int clipDensity(double threshold){
    int    CellsClipped =  0;

    printf("\n\nClipping threshold: %e",    appliedClippingThreshold);

    printf("\nMax overdensity: %f", arrayMax(densityArray, n0*n1*n2));

    for(j=0; j<n0*n1*n2; j++){
        if(densityArray[j]  > threshold){ 
	        densityArray[j] = threshold;
	        
		    CellsClipped   += 1; 
        }
    }

    clippedVolume = CellVolume*CellsClipped;

    printf("\nUnclipped volume: %e", 1. - clippedVolume/TotalVolume);
    
    return 0;
}


int inputHODPk(){
    // Interpolated theoretical P(k) on a regular grid in k. 
    sdltk  = (float *) realloc(sdltk,          470*sizeof(*sdltk));
    sdltPk = (float *) realloc(sdltPk,         470*sizeof(*sdltPk));
              
    // Second derivates of HOD P(k) for cubic spline.           
    sdlt2d = (float *) realloc(sdlt2d,         470*sizeof(*sdlt2d));

    sprintf(filepath, "%s/Data/HODTheoryPk/cambExtendedPk_hod_20.00.dat", root_dir);
    
    inputfile = fopen(filepath, "r");
    
    for(j=0; j<470; j++) fscanf(inputfile, "%f \t %f\n", &sdltk[j], &sdltPk[j]);
    
    fclose(inputfile);

    spline(sdltk, sdltPk, 469, 1.0e31, 1.0e31, sdlt2d);
    
    pt2Pk   = &splintHODpk;
    
    sprintf(theoryPk_flag, "HOD_-20.0");
   
    return 0;
}


int inputLinearPk(){
    // Interpolated theoretical P(k) on a regular grid in k. 
    lineark  = (float *) realloc(lineark,          470*sizeof(*lineark));
    linearPk = (float *) realloc(linearPk,         470*sizeof(*linearPk));
              
    // Second derivates of HOD P(k) for cubic spline.           
    linear2d = (float *) realloc(linear2d,         470*sizeof(*linear2d));

    sprintf(filepath, "%s/Data/HODTheoryPk/camb_matterPk.dat", root_dir);
    
    inputfile = fopen(filepath, "r");
    
    for(j=0; j<470; j++) fscanf(inputfile, "%f \t %f \n", &lineark[j], &linearPk[j]);
    
    for(j=0; j<470; j++) linearPk[j] *= linearBias*linearBias;
    
    fclose(inputfile); 
    
    spline(lineark, linearPk, 469, 1.0e31, 1.0e31, linear2d);
    

    pt2Pk = &splintLinearPk;
    
    sprintf(theoryPk_flag, "linear_-20.0");
    
    return 0;
}


int formPkCube(){
    // hz                 = pow(Om_v + Om_m*pow(1. + 0.8, 3.), 0.5);  // Units of h, [h]
 
    // Om_mz              = Om_m*pow(1. + 0.8, 3.)/pow(hz, 2.);
    
    // f                  = pow(Om_mz, gamma_GR);

    // beta               = f/linearBias;

    // Take line of sight wavevector to be (1, 0, 0), i.e unit vector along kx. 

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kmodulus                           = pow(pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.), 0.5);
                
                mu                                 = k_z/(double) kmodulus;
                
                if(kmodulus == 0.0)          mu    = 0.0;  

		        PkCube[Index]                      = Pk_powerlaw(kmodulus, 5., 1.8);         // Convert from CAMB units for P(k), [P_CAMB] = Volume, to [P(k)] dimensionless.
                PkCube[Index]                     *= 1.0/TotalVolume;
                
                // Impose spherical filter to calculate sigma_8.
                /*
                y                                  = kmodulus*8.;  

                if(kmodulus != 0.0){                
                    PkCube[Index]                 *= 3.*pow(y, -3.)*(sin(y) - y*cos(y));
                    PkCube[Index]                 *= 3.*pow(y, -3.)*(sin(y) - y*cos(y));
                }
                */    
                
                PkCube[Index]                     *= pow(1. + beta*mu*mu, 2.);
                
                // Lorentzian factor for non-linear redshift space distortions. 
                // PkCube[Index]                     /= 1. + 0.5*pow(kmodulus*mu*velDispersion, 2.);
                
                WindowFunc                         = 1.;

                if(k_x != 0.){
		            WindowFunc                    *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);}
                
                if(k_y != 0.){
		            WindowFunc                    *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);}
                
                if(k_z != 0.){
		            WindowFunc                    *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);}		      
	        
	            PkCube[Index]                     *= pow(WindowFunc, 2.);
	        	PkCube[Index]                     *= pow(WindowFunc, 2.);
                
                // Gaussian factor for non-linear redshift space distortion.
                // PkCube[Index]                  *= exp(-kmodulus*kmodulus*mu*mu*velDispersion*velDispersion);
             }
        }
    }
    
    // zero mean density field. 
    PkCube[0] = 0.0;

    return 0;
}


double HermitePolynomial(double x, int n){
    switch(n){
        case 0:
            return 1.0;
        case 1:
            return 2.0*x;
        case 2:
            return 4.0*pow(x, 2.0)   - 2.0;
        case 3:
            return 8.0*pow(x, 3.0)   - 12.0*x;
        case 4:
            return 16.0*pow(x, 4.0)  - 48.*pow(x, 2.0)   + 12.;
        case 5:
            return 32.0*pow(x, 5.0)  - 160.*pow(x, 3.0)  + 120.*x;
        case 6:
            return 64.0*pow(x, 6.0)  - 480.*pow(x, 4.0)  + 720.*pow(x, 2.0)   - 120.;
        case 7:
            return 128.0*pow(x, 7.0) - 1344.*pow(x, 5.0) + 3360.*pow(x, 3.0)  - 1680.*x;
        case 8:
            return 256.*pow(x, 8.)   - 3584.*pow(x, 6.)  + 13440.*pow(x, 4.)  - 13440.*pow(x, 2.)  + 1680.;
        case 9:
            return 512.*pow(x, 9.)   - 9216.*pow(x, 7.)  + 48384.*pow(x, 5.)  - 80640.*pow(x, 3.0) + 30240.*x;
        case 10:
            return 1024.*pow(x, 10.) - 23040.*pow(x, 8.) + 161280.*pow(x, 6.) - 403200.*pow(x, 4.) + 302400.*pow(x, 2.0) - 30240.;   
    }
}


double C_n(double x, int n){                                                // (n+1)! = Gamma (n+2)
    return pow(HermitePolynomial(x, n-1), 2.)*exp(-2.*x*x)/(pi*pow(2., n)*gsl_sf_gamma(n + 2));
}


int clipCorrfn(){
    for(j=0; j<n0*n1*n2; j++) in[j][0] = (double) PkCube[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = (double) 0.0;
   
    printf("\nPerforming FFT.");
    
    fftw_execute(p);

    for(j=0; j<n0*n1*n2; j++) Corrfn[j] = out[j][0];

    double variance;
    double u0;
    
    variance =  Corrfn[0];
    
    // inverse error fn. defined between -1 and 1. | A11Sq | < 1, i.e suppression factor. 
    u0       = appliedClippingThreshold/(sqrt(2.*variance));
    
    printf("\n\nu0:  %e, variance: %e", u0, variance);
    
    for(j=0; j<n0*n1*n2; j++) suppressedCorrfn[j]      = 0.25*pow(1.0 + gsl_sf_erf(u0), 2.)*Corrfn[j]; 
    for(j=0; j<n0*n1*n2; j++)  distortedCorrfn[j]      = suppressedCorrfn[j]; 
    
    for(i=1; i<10; i++){      
        for(j=0; j<n0*n1*n2; j++){  
            distortedCorrfn[j]   += variance*pow(Corrfn[j]/variance, i + 1)*C_n(u0, i);
        }
    }
    
    iplan   = fftw_plan_dft_3d(n0, n1, n2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = distortedCorrfn[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    fftw_execute(iplan);

    for(j=0; j<n0*n1*n2; j++) clippedPk[j] = pow(n0*n1*n2, -1.)*out[j][0];

    printf("\n\nClipped P(k).");
    
    ClippedMultipole();
   
    return 0;    
}


double splint_kMonopole(double k){
    float Interim;
    
    splint(f_meanKBin, f_kMonopole, f_kMonopole2d, kBinNumb-2, (float) k, &Interim);

    return (double) Interim*TotalVolume;
}


double splint_kQuadrupole(double k){
    float Interim;
    
    splint(f_meanKBin, f_kQuadrupole, f_kQuadrupole2d, kBinNumb-2, (float) k, &Interim);
    
    return (double) Interim*TotalVolume;
}


int ClippedMultipole(){
    polarPk_modeCount = 0;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;  
                
                WindowFunc                         = 1.;

                if(k_x != 0.){
		            WindowFunc                    *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);}
                
                if(k_y != 0.){
		            WindowFunc                    *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);}
                
                if(k_z != 0.){
		            WindowFunc                    *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);}

                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2.
                
                clippedPk[Index]                  /= pow(WindowFunc, 2.)*pow(WindowFunc, 2.);

                if(kmodulus > 0.000001){	            
		            // Issue with mu for a zeroth length vector being ill defined. 
		            polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		            polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		            polar2Dpk[polarPk_modeCount][2]    = clippedPk[Index];

                    // printf("\n%e \t %e \t %e", kmodulus, fabs(mu), clippedPk[Index]);

                    polarPk_modeCount                 += 1;
	            }
            }
        }
    }
    
    polar2DpkBinning(polarPk_modeCount);
    
    sprintf(filepath, "%s/Data/Multipoles/Multipoles_clippedPrediction_%s_%s_beta_%.2f_sigma_%.2f_kbin_%.3f.dat", root_dir, surveyType, theoryPk_flag, beta, velDispersion, kbinInterval);

    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, modesPerBin, polar2Dpk, polarPk_modeCount, filepath);
    /*
    for(j=0; j<kBinNumb; j++){
        f_meanKBin[j]    = (float)  meanKBin[j];
    
        f_kMonopole[j]   = (float)  kMonopole[j];
        
        f_kQuadrupole[j] = (float) kQuadrupole[j];
    }  
    
    spline(f_meanKBin, f_kMonopole,   kBinNumb-2, 1.0e31, 1.0e31, f_kMonopole2d);
    
    spline(f_meanKBin, f_kQuadrupole, kBinNumb-2, 1.0e31, 1.0e31, f_kQuadrupole2d);
    */
    return 0;
}
