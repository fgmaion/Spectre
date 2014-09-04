double Si_minus_Si_QuadTerm_usingTaylorAboutkr(double r, double k, double d){
    // Quadrupole corr. fns correponding to a unit step in P_2(k) contains a term differing Sine integrals, Si(r(k+d)) - Si(r(k-d)).  
    // Treat this by Taylor expanding both functions about k*r. 

    Interim  = 0.0;
    
    Interim += 2.*d*sin(k*r)/k;
    
    Interim -= 2.*pow(d, 3.)*(2.*k*r*cos(k*r) + (k*r*r - 2.)*sin(k*r))/(6.*k);
    
    Interim += 2.*pow(d, 5.)*(4.*k*r*(k*r*r - 6.)*cos(k*r) + (24. - 12.*k*r*r + k*pow(r, 4.))*sin(k*r))/(120.*k);
    
    Interim -= 2.*pow(d, 7.)*(6.*k*r*(120. - 20.*k*r*r + k*pow(r, 4.)*cos(k*r) + (-720. + 360.*k*r*r -30.*k*pow(r, 4.) + k*pow(r, 6.))*sin(k*r)))/(5040.*k);
    
    return Interim;
}


double unitTheory(double r, double k, double d, int order){
  // Evaluates the monopole/quadrupole correlation function (order={0,2}) corresponding to a unit theory vector, it P_l(k) = 1. for k_l<k<k_u, =0 otherwise. 

  switch(order){
    case 0:
      return  pow(pi, -2.)*(pow(r, -2.)*(k*sin(k*r)*sin(d*r) - d*cos(k*r)*cos(d*r)) + pow(r, -3.)*cos(k*r)*sin(d*r));

    case 2:
      return (-1./(2.*pi*pi))*(3.*pow(r, -3.)*Si_minus_Si_QuadTerm_usingTaylorAboutkr(r, k, d) -2.*pow(r, -2.)*(k*sin(k*r)*sin(d*r) - d*cos(k*r)*cos(d*r) + (1./r)*cos(k*r)*sin(d*r)) -6.*pow(r, -3.)*cos(k*r)*sin(d*r));
    }
}


int calc_mixingmatrix(double kbinLimits[], int kBinNumb){
    // theta is the column vector consisting of monopole then quadrupole.
    // kbinLimits contain the limit of each k bin in which both monopole and quadrupole have been binned.  

    // theta*_i represents convolved (monopole, quadrupole)^T
    // theta*_i = Sum M_ij theta_j

    int mm, nn, transformOrder;
    
    free2DBinning();

    assign2DPkMemory(muBinNumb, kBinNumb);

    assign_mixingmatrix(kBinNumb, kMonopole, kQuadrupole);
    
    int m0, m1, m2;
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = (j+1)*kbinInterval;

    // unsort rmodulus_vec
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               m0 = k;
               m1 = j;
               m2 = i;

               if(m2>n2/2)        m2 -= n2;
               if(m1>n1/2)        m1 -= n1;
               if(m0>n0/2)        m0 -= n0;

               Index                  = k*n1*n2 + j*n2 + i;
               
               rmodulus_vec[Index][0] =           CellSize*sqrt(m2*m2 + m1*m1 + m0*m0);
               rmodulus_vec[Index][1] = (double)                                 Index;
            }
        }
    }
    
    //append an 'end of array' element to exit loop successfully. 
    rmodulus_vec[n0*n1*n2][0] =           9999.;
    rmodulus_vec[n0*n1*n2][1] = (double)   9999;
    
    printf("\n\nCalculating mixing matrix.\n");
    
    for(transformOrder=0; transformOrder<4; transformOrder += 2){
        for(nn=0; nn<(kBinNumb-1); nn++){         
            printf("%d \t %d \n", transformOrder, nn);
            
            // calculation of correlation function for given unit theory vector, 1 for the k interval, zero otherwise, either monopole or quadrupole. 
            for(k=0; k<n0*n1*n2; k++)       Corrfn[k]   = unitTheory(rmodulus_vec[k][0], kBinLimits[nn], 0.5*kbinInterval, transformOrder)*TotalVolume;
            
                                            Corrfn[0]   = 0.0;
            
            // Carry out the convolution by inverse FFT of the correlation fn. x FFT of W^2(k).  Renormalise to prevent underflow. 
            for(ii=0; ii<n0*n1*n2; ii++) in[ii][0]     = Corrfn[ii]*FFTW2_vecr_re[ii]*LegendrePolynomials(mu_vec[ii], transformOrder);
            for(ii=0; ii<n0*n1*n2; ii++) in[ii][1]     = Corrfn[ii]*FFTW2_vecr_im[ii]*LegendrePolynomials(mu_vec[ii], transformOrder);
    
            // Calculate the convolved P(vec k).
            fftw_execute(p);
                
            // Set up for multipole calculation. 
            for(kk=0; kk<n0*n1*n2; kk++){
                polar2Dpk[kk][0]  =                  kmodulus_vec[kk];
	    	    polar2Dpk[kk][1]  =                        mu_vec[kk];
	    	    polar2Dpk[kk][2]  =     pow(n0*n1*n2, -1.0)*out[kk][0];
            }
    
            // Calculate monopole and quadrupole of convolved P(k) resulting from unit theory. 
            MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polarPk_modeCount, filepath, 0.0, 1.0, kbinInterval, 0);
        
            for(mm=0; mm<(kBinNumb-1); mm++){
                mixingmatrix[mm][(transformOrder/2)*(kBinNumb-1) + nn]                =   kMonopole[mm];
           
                mixingmatrix[(kBinNumb-1) + mm][(transformOrder/2)*(kBinNumb-1) + nn] = kQuadrupole[mm];  
                
                printf("\n%e \t %e", kMonopole[mm], kQuadrupole[mm]);
            }  
        }        
    }     
    
    sprintf(filepath, "%s/Data/SpectralDistortion/mixingmatrix_VipersMask.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<2*(kBinNumb-1); j++){
      for(k=0; k<2*(kBinNumb-1); k++){
        fprintf(output, "%e\n", mixingmatrix[j][k]);       
      }
    }
    
    fclose(output);
    
    return 0;
}


int prepConvolution(double kbinLimits[], int kBinNumb){
    free2DBinning();

    assign2DPkMemory(muBinNumb, kBinNumb);

    assign_mixingmatrix(kBinNumb, kMonopole, kQuadrupole);

    sprintf(filepath, "%s/Data/SpectralDistortion/mixingmatrix_VipersMask.dat", root_dir);
    
    inputfile = fopen(filepath, "r");
    
    for(j=0; j<2*(kBinNumb-1); j++){
        for(k=0; k<2*(kBinNumb-1); k++){
            fscanf(inputfile, "%le \n", &mixingmatrix[j][k]);       
        }
    }
    
    fclose(inputfile);

    return 0;
}


int convolvePk(){
    sprintf(filepath, "%s/Data/SpectralDistortion/mixingmatrix_convolvedPk_VipersMask_IntConCorr.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    int    transformOrder;
    
    double MonoQuadleak;
    
    for(j=0; j<kBinNumb; j++)  kBinLimits[j]  = (j+1)*kbinInterval;
    
    for(j=0; j<(kBinNumb-1); j++){
        convMono[j]         = 0.0;
        convQuad[j]         = 0.0;
        
        MonoQuadleak = 0.0;
        
        for(transformOrder=0; transformOrder<4; transformOrder += 2){
            for(k=0; k<(kBinNumb-1); k++)  convMono[j] += (float) mixingmatrix[j][(transformOrder/2)*(kBinNumb-1) + k]*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, transformOrder);
        
            for(k=0; k<(kBinNumb-1); k++)  convQuad[j] += (float) mixingmatrix[(kBinNumb-1) + j][(transformOrder/2)*(kBinNumb-1) + k]*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, transformOrder);
        }
        
        // for(k=0; k<(kBinNumb-1); k++)  MonoQuadleak += mixingmatrix[(kBinNumb-1) + j][k]*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, 0);
        
        // fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e\n", kBinLimits[j], (*pt2Pk)(kBinLimits[j])*(*pt2RSD_k)(kBinLimits[j]*velDispersion, beta, 0), (*pt2Pk)(kBinLimits[j])*(*pt2RSD_k)(kBinLimits[j]*velDispersion, beta, 2), Mono, Quad, MonoQuadleak);
    }
    
    // for(j=0; j<kBinNumb-1; j++)  fkBinLimits[j] = (float) kBinLimits[j];

    // spline(fkBinLimits, convMono, kBinNumb-2, 1.0e31, 1.0e31, convMono2d);

    // spline(fkBinLimits, convQuad, kBinNumb-2, 1.0e31, 1.0e31, convQuad2d);
    
    for(j=0; j<(kBinNumb-1); j++)  fprintf(output, "%e \t %e \t %e \n", kBinLimits[j], convMono[j], convQuad[j]);

    fclose(output);

    return 0;
}


double splintConvMono(double k){
    float Interim;

    splint(fkBinLimits, convMono, convMono2d, kBinNumb-2, (float) k, &Interim);

    return (double)  Interim;
}


double splintConvQuad(double k){
    float Interim;

    splint(fkBinLimits, convQuad, convQuad2d, kBinNumb-2, (float) k, &Interim);

    return (double)  Interim;
}


int printf_unitTheoryxi(){
    sprintf(filepath, "%s/Data/SpectralDistortion/unitTheoryxi_5kvals_quad.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double dj, kj, Interim;
    
    for(j=0; j<1000; j++){ 
        dj  =   (double) j;
        dj /=         200.;
        
        dj -=           4.;
        
        dj  = pow(10., dj);
        
        fprintf(output, "%e \t ", dj);
    
        for(k=0; k<1; k++){  
            kj = (double) k;
                        
            kj -=        2.;
            
            fprintf(output, "%e \t ", unitTheory(dj, pow(10., kj), 0.5*kbinInterval, 2));
        }
    }
    
    fclose(output);

    return 0;
}


int printf_xi(){
    sprintf(filepath, "%s/Data/SpectralDistortion/unitTheory_xi_mono_quad.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double dj, kj, mono, quad;
    
    for(j=0; j<kBinNumb; j++)      kBinLimits[j]  =     (j+1)*kbinInterval;
    
    for(j=0; j<1000; j++){ 
        dj      =   (double) j;
        dj     /=         200.;
        
        dj     -=           2.;
        
        dj      = pow(10., dj);
        
        mono = 0.0;
        quad = 0.0;
    
        // Cannot use kMonopole values due to aliased measurement. 
        for(k=0; k<(kBinNumb-1); k++) mono += unitTheory(dj, kBinLimits[k], 0.5*kbinInterval, 0)*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, 0);
        
        for(k=0; k<(kBinNumb-1); k++) quad += unitTheory(dj, kBinLimits[k], 0.5*kbinInterval, 2)*(*pt2Pk)(kBinLimits[k])*(*pt2RSD_k)(kBinLimits[k]*velDispersion, beta, 2);

        fprintf(output, "%e \t %e \t %e \n", dj, mono, quad);
    }
    
    fclose(output);
    
    return 0;
}
