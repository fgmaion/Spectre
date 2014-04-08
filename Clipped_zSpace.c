int clipDensity(double threshold){
    int CellsClipped = 0;

    for(j=0; j<n0*n1*n2; j++){
        if(densityArray[j] > threshold){ 
	        densityArray[j] = threshold;
		    CellsClipped   += 1; 
        }
    }
    
    double clippedVolume;

    clippedVolume = CellVolume*CellsClipped;

    printf("\nClipped volume: %e", clippedVolume/TotalVolume);

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
    
    return 0;
}


int formPkCube(){
    PkCube             = (double *) realloc(PkCube,          n0*n1*n2*sizeof(*PkCube));

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

		        PkCube[Index]                      = (*pt2Pk)(kmodulus);
  
                // Convert from CAMB units for P(k), [P_CAMB] = Volume, to [P(k)] dimensionless.
                PkCube[Index]                     *= 1.0/TotalVolume;

                // Impose spherical filter to calculate sigma_8.
                /*
                y                                  = kmodulus*8.;  

                if(kmodulus != 0.0){                
                    PkCube[Index]                 *= 3.*pow(y, -3.)*(sin(y) - y*cos(y));
                    PkCube[Index]                 *= 3.*pow(y, -3.)*(sin(y) - y*cos(y));
                }
                */
                
                mu                                 = k_x/(double) kmodulus;
                
                if(kmodulus == 0.0)          mu   = 0.0;      
                
                kaiserFactor                       = pow(1. + beta*mu*mu, 2.);
                
                PkCube[Index]                     *= kaiserFactor;
                
                // Lorentzian factor for non-linear redshift space distortions. 
                // PkCube[Index]                     /= 1. + 0.5*pow(kmodulus*mu*velDispersion, 2.);
                
                // Gaussian factor for non-linear redshift space distortion.
                PkCube[Index]                     *= exp(-kmodulus*kmodulus*mu*mu*velDispersion*velDispersion);
                
                polar2Dpk[Index][0]                = kmodulus;
		        polar2Dpk[Index][1]                = fabs(mu);
		        polar2Dpk[Index][2]                = PkCube[Index];
             }
        }
    }

    printf("\n\nMinimum Kaiser factor:  %e", pow(1. + 0., 2.));
    printf("\nMaximum Kaiser factor:  %e", pow(1. + beta, 2.));

    return 0;
}


int theoryHexadecapole(){
    
    HexadecapoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, kHexadecapole, modesPerBin, polar2Dpk);

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
    Corrfn            = (double *) realloc(Corrfn,             n0*n1*n2*sizeof(*Corrfn));
    suppressedCorrfn  = (double *) realloc(suppressedCorrfn,   n0*n1*n2*sizeof(*suppressedCorrfn));
    distortedCorrfn   = (double *) realloc(distortedCorrfn,    n0*n1*n2*sizeof(*distortedCorrfn));
    clippedPk         = (double *) realloc(clippedPk,          n0*n1*n2*sizeof(*clippedPk));

    for(j=0; j<n0*n1*n2; j++) in[j][0] = (double) PkCube[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = (double) 0.0;
   
    printf("\nPerforming FFT.");
    
    fftw_execute(p);

    for(j=0; j<n0*n1*n2; j++) Corrfn[j] = out[j][0];

    double variance    =  Corrfn[0];
    double u0;

    u0      = inverse_erf(2.*sqrt(A11) -1.);
    
    printf("\nu0 estimated  from suppression factor:  %f", u0);
    
    for(j=0; j<n0*n1*n2; j++) suppressedCorrfn[j]      = 0.25*pow(1.0 + gsl_sf_erf(u0), 2.)*Corrfn[j]; 
    for(j=0; j<n0*n1*n2; j++) distortedCorrfn[j]       = suppressedCorrfn[j]; 
    
    for(i=1; i<10; i++){
                                 distortedCorrfn[0]   += variance*C_n(u0, i);
      for(j=1; j<n0*n1*n2; j++)  distortedCorrfn[j]   += variance*pow(Corrfn[j]/variance, i + 1)*C_n(u0, i);
    }
    
    iplan   = fftw_plan_dft_3d(n0, n1, n2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(j=0; j<n0*n1*n2; j++) in[j][0] = distortedCorrfn[j]/A11;
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    fftw_execute(iplan);

    for(j=0; j<n0*n1*n2; j++) clippedPk[j] = pow(n0*n1*n2, -1.)*out[j][0];

    printf("\n\nClipped P(k).");
     
    return 0;    
}
