double slowDFT(double kx, double ky, double kz){
    double Real  = 0.0;
    double Imag  = 0.0;
    double theta = 0.0;
    double Pk    = 0.0;

    for(jj=0; jj<Vipers_Num; jj++){
       if(Acceptanceflag[jj] == true){ 
           theta = kx*xCoor[jj] + ky*yCoor[jj] + kz*zCoor[jj];         
                    
           Real += cos(theta);
           Imag += sin(theta);
       }
    }

    // spec mock 001. 
    Real /= 20768.;
    Imag /= 20768.;
    
    return pow(Real, 2.) + pow(Imag, 2.);
}


int slowDFTcalc(){
    polar_pkcount = 0;
    
    double pk, GaussianFilter, WindowFunc;
    
    prep_pkRegression(-2., log10(modkMax), kbin_no);
    
    for(k=0; k<n0; k+=10){
        printf("\n%.2lf percentage complete.", 100.*k/n0);
        
        for(j=0; j<n1; j+=10){
            for(i=0; i<n2; i+=10){
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
                
                pk                                 = slowDFT(k_x, k_y, k_z);
                    
                // Account for the affect of the window on the amplitude of the measured P(k).
                pk                                /= fkpSqWeightsVolume*pow(TotalVolume, -1.);
                
                if(kmodulus > 0.000001){
	                // Only half the modes are independent. 
	            	if(k_z>0.){
	            	    // One hemi-sphere is independent, e.g. k_z >= 0.
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk;
		            
		                polar_pkcount               += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk;
		            
		                polar_pkcount                += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
		                // on the line k_z=k_y=0, one half is independent, k_x>=0.
		                                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]    = kmodulus;
		                polar_pk[polar_pkcount][1]    = fabs(mu);
		                polar_pk[polar_pkcount][2]    = pk;
		            
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
