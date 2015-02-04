int CatalogueInput_Cube(char filepath[]){
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    
    if(inputfile == NULL){  
        printf("\nError opening %s\n", filepath); 
        
        return 1;
    }
    
    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    xCoor          =  (double *)  malloc(Vipers_Num*sizeof(*xCoor));
    yCoor          =  (double *)  malloc(Vipers_Num*sizeof(*yCoor));
    zCoor          =  (double *)  malloc(Vipers_Num*sizeof(*zCoor));
    
    for(j=0; j<Vipers_Num; j++)  fscanf(inputfile, "%le \t %le \t %le \t %*le \t %*le \t %*le \t %*le \t %*d \t %*d \t %*le \n", &xCoor[j], &yCoor[j], &zCoor[j]);
    
    // for(j=0; j<10; j++)  printf("\n%.4lf \t %.4lf \t %.4lf", xCoor[j], yCoor[j], zCoor[j]);
    
    fclose(inputfile);
    
    printf("\nHOD Cube catalogue input successful.");
    
    return 0;
}


int calc_overdensityCube(){
    prep_grid();
    
    for(j=0; j<n0*n1*n2; j++)           overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)           overdensity[j][1] = 0.0;

    for(j=0; j<Vipers_Num; j++){
        // no selection necessary. 
        boxlabel = boxCoordinates(xCoor, yCoor, zCoor, j);
            
        overdensity[boxlabel][0]   += 1.;
    }
    
    double cube_nbar = Vipers_Num/TotalVolume;

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] /= CellVolume*cube_nbar;
    
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] -=                  1.0;

    return 0;
}


int cube_PkCorrections(){
    polar_pkcount = 0;

    double pk, GaussianFilter, WindowFunc;

    double cube_nbar = Vipers_Num/TotalVolume;

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
                
                GaussianFilter                     = exp(-1.*kSq*0.5*(pow(1., 2.)));

                WindowFunc                         = 1.;

                if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);
                
                if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);
                
                if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);
                
                H_k[Index][0]                     *= pow(n0*n1*n2, -1.0);
                H_k[Index][1]                     *= pow(n0*n1*n2, -1.0);
                
                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                H_k[Index][0]                     /= pow(WindowFunc, 1.);
                H_k[Index][1]                     /= pow(WindowFunc, 1.);
                
                pk                                 = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);
        
                // Rescale measured P(k) in the cube to the V=1 camb convention.
                pk                                *= TotalVolume; 
                    
                pk                                -= 1./cube_nbar;

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
    
    return 0;
}


int cube_PkCalc(){
    // assign memory for H_k and plan. 
    prep_fftw();
    
    fftw_execute(p);
    
    // free overdensity.
    free_grid();
    
    prep_pkRegression(-2., log10(modkMax), kbin_no);
    
    cube_PkCorrections();
    
    
    sprintf(filepath,"%s/Data/500s/HOD_mocks/cube_-20.0_512.dat", root_dir);

    MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 1);

    return 0;
}
