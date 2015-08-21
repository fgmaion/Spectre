int gal_slowDFT(double kx, double ky, double kz, double* Real, double* Imag){
    *Real  = 0.0;
    *Imag  = 0.0;
    
    double theta = 0.0;
    double weight;

    for(jj=0; jj<Vipers_Num; jj++){
       if(Acceptanceflag[jj] == true){ 
           theta  = kx*xCoor[jj] + ky*yCoor[jj] + kz*zCoor[jj];         
                    
           weight = sampling[jj];
                    
           *Real += weight*cos(theta);
           *Imag += weight*sin(theta);
       }
    }
    
    // weight normalisation. 
    *Real /= sqrt(alpha);       
    *Real /= fkp_weight;     
    
    *Imag /= sqrt(alpha);       
    *Imag /= fkp_weight;         
    
    // printf("\n%e \t %e", *Real, *Imag);
    
    return 0;
}


int rand_slowDFT(double kx, double ky, double kz, double* Real, double* Imag){
    *Real  = 0.0;
    *Imag  = 0.0;
    
    double theta = 0.0;
    double weight;

    for(jj=0; jj<rand_number; jj++){
       if(rand_accept[jj] == true){ 
           theta   = kx*rand_x[jj] + ky*rand_y[jj] + kz*rand_z[jj];         
                
            weight = rand_weight[jj];  
                    
           *Real  += weight*cos(theta);
           *Imag  += weight*sin(theta);
       }
    }
    
    // alpha = 1 normalisation. 
    *Real /= fkp_weight;     
    
    // alpha = 1 normalisation. 
    *Imag /= fkp_weight;        
    
    return 0;
}


int randoms_slowDFTcalc(){
    // change in gal_slowDFT aswell. 
    int mode_number = 20;
    
    double logk_interval, real, imag;
    
    Rr = malloc(mode_number*mode_number*mode_number*sizeof(double));
    Ir = malloc(mode_number*mode_number*mode_number*sizeof(double));
                
    logk_interval = 2./mode_number;
    
    printf("\n\nCalculating randoms slow DFT calc:");    
    
    for(k=0; k<mode_number; k++){    
        printf("\n%d", k);
        
        for(j=0; j<mode_number; j++){
            for(i=0; i<mode_number; i++){
                k_x    = pow(10., -2. + logk_interval*i);
                k_y    = pow(10., -2. + logk_interval*j);
                k_z    = pow(10., -2. + logk_interval*k);
                
                Index  = k*mode_number*mode_number + j*mode_number + i;
                
                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                if(kmodulus > 0.000001){
	                // Only half the modes are independent. 
	            	if(k_z>0.){       
                        rand_slowDFT(k_x, k_y, k_z, &real, &imag);
                        
                        Rr[Index] = real;
                        Ir[Index] = imag;
		            }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        rand_slowDFT(k_x, k_y, k_z, &real, &imag);
                       
                        Rr[Index] = real;
                        Ir[Index] = imag;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){ 
                        rand_slowDFT(k_x, k_y, k_z, &Rr, &Ir);

                        Rr[Index] = real;
                        Ir[Index] = imag;
		            }
		            		            
		            // else no dice.    
	            }
            }
        }
    }
    
    return 0;
}   


int slowDFTcalc(){
    polar_pkcount   =  0;
    
    // change in rand_slowDFT aswell. 
    int mode_number = 20;
    
    double pk, GaussianFilter, WindowFunc, logk_interval, Rg, Ig, Shot;
    
    prep_pkRegression(-2., log10(modkMax), kbin_no);
    
    logk_interval = 2./mode_number;
    
    
    fkp_weight   = fkp_weightsnorm();
    
    Shot         = fkp_shotnoise();
    
    printf("\n\nCalculating slow DFT calc:");    
    
    for(k=0; k<mode_number; k++){    
        printf("\n%d", k);
        
        for(j=0; j<mode_number; j++){
            for(i=0; i<mode_number; i++){
                k_x    = pow(10., -2. + logk_interval*i);
                k_y    = pow(10., -2. + logk_interval*j);
                k_z    = pow(10., -2. + logk_interval*k);
                
                Index  = k*mode_number*mode_number + j*mode_number + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                if(kmodulus > 0.000001){
	                // Only half the modes are independent. 
	            	if(k_z>0.){
	            	    // One hemi-sphere is independent, e.g. k_z >= 0.
	            	    mu                         = k_z/kmodulus;
                        if(kmodulus < 0.000001) mu = 0.0;      
                
                         gal_slowDFT(k_x, k_y, k_z, &Rg, &Ig);

                        pk                         = pow(Rg - alpha*fkp_weight*Rr[Index], 2.) + pow(Ig - alpha*fkp_weight*Ir[Index], 2.);
                    
                        // Account for the affect of the window on the amplitude of the measured P(k).
                        // pk                        /= fkpSqWeightsVolume*pow(TotalVolume, -1.);
	            	    
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk - Shot;
		            
		                polar_pkcount               += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
                        mu                         = k_z/kmodulus;
                        if(kmodulus < 0.000001) mu = 0.0;      
                
                
                        gal_slowDFT(k_x, k_y, k_z, &Rg, &Ig);
                        
                        rand_slowDFT(k_x, k_y, k_z, &Rr, &Ir);

                        pk                         = pow(Rg - alpha*Rr[Index], 2.) + pow(Ig - alpha*Ir[Index], 2.);
                    
                        // Account for the affect of the window on the amplitude of the measured P(k).
                        // pk                        /= fkpSqWeightsVolume*pow(TotalVolume, -1.);
                        
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk - Shot;
		            
		                polar_pkcount                += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
		                // on the line k_z=k_y=0, one half is independent, k_x>=0.
		                // in the k_z=0 plane one semi-circle is independent, k_y>0.

                    	mu                         = k_z/kmodulus;
                        if(kmodulus < 0.000001) mu = 0.0;      
                
                         gal_slowDFT(k_x, k_y, k_z, &Rg, &Ig);
                        
                        rand_slowDFT(k_x, k_y, k_z, &Rr, &Ir);

                        pk                         = pow(Rg - alpha*Rr[Index], 2.) + pow(Ig - alpha*Ir[Index], 2.);
                    
                        // Account for the affect of the window on the amplitude of the measured P(k).
                        // pk                        /= fkpSqWeightsVolume*pow(TotalVolume, -1.);		                
		                
		                polar_pk[polar_pkcount][0]    = kmodulus;
		                polar_pk[polar_pkcount][1]    = fabs(mu);
		                polar_pk[polar_pkcount][2]    = pk - Shot;
		            
		                polar_pkcount                 += 1;
		            }
		            		            
		            // else no dice.    
	            }
            }
        }
    }
                
    observedQuadrupole(polar_pkcount);
    
    return 0;
}   


int modeSampling(double length){
    // Monte Carlo sample the modes supplied by an FFT of the box with 
    // fundamental mode 2.pi/length.
    double dk = 2.*pi/length;

    double draw, dk2, real, imag, Rg, Ig, Shot, pk;

    int    maxint, jjj, numModes;

    polar_pkcount = 0;

    numModes = 2000;

    // fundamental period to k of 1. 
    maxint = (int) ceil(0.6/dk);

    dk2 = dk*dk;
    
    prep_pkRegression(-2., log10(modkMax), kbin_no);

    double   kx[numModes];
    double   ky[numModes];
    double   kz[numModes];
    double  kmu[numModes];
    double kmod[numModes];
    
    Rr = malloc(numModes*sizeof(double));
    Ir = malloc(numModes*sizeof(double));

    for(j=0; j<rand_number; j++) rand_weight[j] = 1./interp_nz(rand_chi[j]);

    // normalisation of the weights. assuming alpha = 1
    fkp_weight = fkp_weightsnorm();
    
    Shot       =   fkp_shotnoise();
    
    printf("\n\nfkp weight norm: %e, shot: %e  (alpha = 1 for both)", fkp_weight, Shot);
    
    while(polar_pkcount<numModes){
        i      = gsl_rng_uniform_int(gsl_ran_r, maxint);
        j      = gsl_rng_uniform_int(gsl_ran_r, maxint);
        k      = gsl_rng_uniform_int(gsl_ran_r, maxint);

        k_x    = dk*i;
        k_y    = dk*j;
        k_z    = dk*k;
    
        kSq    = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
    
        draw = gsl_rng_uniform(gsl_ran_r);
    
        if(draw < dk2/kSq){
            kmodulus = pow(kSq, 0.5);
        
            mu = k_z/kmodulus;
        
            if(kmodulus < 0.000001) mu = 0.0;      
                
            kx[polar_pkcount]        = k_x;
            ky[polar_pkcount]        = k_y;
            kz[polar_pkcount]        = k_z;
            
            rand_slowDFT(k_x, k_y, k_z, &real, &imag);
            
            Rr[polar_pkcount] = real;
            Ir[polar_pkcount] = imag;
	            	    
		    kmod[polar_pkcount]   = kmodulus;
		     kmu[polar_pkcount]   = fabs(mu);
		            
		    // printf("\n%e \t %e \t %e", polar_pk[polar_pkcount][0], polar_pk[polar_pkcount][1], polar_pk[polar_pkcount][2]);
		            
		    polar_pkcount               += 1;
		    
		    if(polar_pkcount%100 == 0)  printf("\n%d", polar_pkcount);
        }   
    }

    for(loopCount=1; loopCount<307; loopCount++){
        printf("\n\n%d", loopCount);

        // new 500s. limit to 0.7<z<0.8, limit to linear bias of 1.495903 corresponding to ~ -20.0 mag galaxies. at this mag. volume limited to z=0.85
        if(loopCount<10)        sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_00%d_specmask.dat", vipersHOD_dir, loopCount);
        else if(loopCount<100)  sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_0%d_specmask.dat",  vipersHOD_dir, loopCount);
        else                    sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_%d_specmask.dat",   vipersHOD_dir, loopCount);
    
        CatalogueInput_500s(filepath);

        // Choice of redshift from zcos, zpec, zphot, zobs.
        gal_z = &zobs[0];
    
        // Set redshift and absolute mag. cuts. 
        assignAcceptance();

        // Convert from (ra, dec, redshift) to (x, y, z) in Stefano's basis. basis choice must be consistent with that used for the mask defined by the randoms. 
        StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
        
        for(j=0; j<Vipers_Num; j++)  sampling[j] = 1./interp_nz(rDist[j]);
        
        printf("\n\nalpha: %e, fkp weight norm: %e, shot: %e", alpha, fkp_weight, (1. + alpha)*Shot);
        
        printf("\n\nCalculating monte carlo mode sampled slow DFT calc:");    
        
        for(jjj=0; jjj<numModes; jjj++){
             gal_slowDFT(kx[jjj], ky[jjj], kz[jjj], &Rg, &Ig);   
            
            pk                 = pow(Rg - alpha*(Rr[jjj]/sqrt(alpha)), 2.) + pow(Ig - alpha*(Ir[jjj]/sqrt(alpha)), 2.);
    
            // polar pk is sorted in MultipoleCalc, must be reinitialised. 
            polar_pk[jjj][0]   =              kmod[jjj];
            polar_pk[jjj][1]   =               kmu[jjj];
            polar_pk[jjj][2]   = pk - (1. + alpha)*Shot;
        
            // printf("\n%e", pk);
        }
        
        observedQuadrupole(numModes);
    }
    
    return 0;
}


int vipers_fkpCalc(double length){
    // Monte Carlo sample the modes supplied by an FFT of the box with 
    // fundamental mode 2.pi/length.
    double dk = 2.*pi/length;

    double draw, dk2, real, imag, Rg, Ig, Shot, pk;

    int    maxint, jjj, numModes;

    polar_pkcount = 0;

    numModes = 20000;

    // fundamental period to k of 1. 
    maxint = (int) ceil(1.0/dk);

    dk2 = dk*dk;
    
    prep_pkRegression(-2., log10(modkMax), kbin_no);

    double   kx[numModes];
    double   ky[numModes];
    double   kz[numModes];
    double  kmu[numModes];
    double kmod[numModes];
    
    Rr = malloc(numModes*sizeof(double));
    Ir = malloc(numModes*sizeof(double));

    for(j=0; j<rand_number; j++) rand_weight[j] = 1.; // /interp_nz(rand_chi[j]);

    // normalisation of the weights. assuming alpha = 1, alpha scaling dealth with later. 
    fkp_weight = fkp_weightsnorm();
    
    Shot       =   fkp_shotnoise();
    
    printf("\n\nfkp weight norm: %e, shot: %e  (alpha = 1 for both)", fkp_weight, Shot);
        
    double modk, mu, phi;
         
    while(polar_pkcount<numModes){   
        i      = gsl_rng_uniform_int(gsl_ran_r, maxint);
        j      = gsl_rng_uniform_int(gsl_ran_r, maxint);
        k      = gsl_rng_uniform_int(gsl_ran_r, maxint);
    
        k_x    = dk*i;
        k_y    = dk*j;
        k_z    = dk*k;

        kSq    = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
    
        if( (kSq>pow(10., -4.)) && (kSq<1.)){
          draw = gsl_rng_uniform(gsl_ran_r);
    
          if(draw < dk2/kSq){
            kmodulus = pow(kSq, 0.5);
        
            mu = k_z/kmodulus;
        
            if(kmodulus < 0.000001) mu = 0.0;      
                
            kx[polar_pkcount]        = k_x;
            ky[polar_pkcount]        = k_y;
            kz[polar_pkcount]        = k_z;
            
            rand_slowDFT(k_x, k_y, k_z, &real, &imag);
            
            Rr[polar_pkcount] = real;
            Ir[polar_pkcount] = imag;
	            	    
	    	kmod[polar_pkcount]   = kmodulus;
	    	kmu[polar_pkcount]   = fabs(mu);
		            
	    	polar_pkcount               += 1;
		    
	    	if(polar_pkcount%100 == 0) printf("\n%d", polar_pkcount);
          }
        }   
    }
    
    
    sprintf(filepath, "%s/W1_Spectro_V5_0/W1_SPECTRO_V5_0.txt", root_dir);
  
    DataInput(filepath);

    // Choice of redshift from zcos, zpec, zphot, zobs.
    gal_z = &zobs[0];
    
    // Set redshift and absolute mag. cuts. 
    assignAcceptance();

    // Convert from (ra, dec, redshift) to (x, y, z) in Stefano's basis. basis choice must be consistent with that used for the mask defined by the randoms. 
    StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
        
    // weight w used in fkp calc. 
    for(j=0; j<Vipers_Num; j++)  sampling[j] = 1.;
        
    
    printf("\n\nalpha: %e, fkp weight norm: %e, shot: %e", alpha, fkp_weight, (1. + alpha)*Shot);
        
    printf("\n\nCalculating monte carlo mode sampled slow DFT calc:");    
        
    for(jjj=0; jjj<numModes; jjj++){
        gal_slowDFT(kx[jjj], ky[jjj], kz[jjj], &Rg, &Ig);   
            
        pk                 = pow(Rg - alpha*(Rr[jjj]/sqrt(alpha)), 2.) + pow(Ig - alpha*(Ir[jjj]/sqrt(alpha)), 2.);
    
        // polar pk is sorted in MultipoleCalc, must be reinitialised. 
        polar_pk[jjj][0]   =              kmod[jjj];
        polar_pk[jjj][1]   =               kmu[jjj];
        polar_pk[jjj][2]   = pk - (1. + alpha)*Shot;
        
	    if(jjj%100 == 0) printf("\n%d", jjj);
    }
        
    observedQuadrupole(numModes);
    
    return 0;
}
