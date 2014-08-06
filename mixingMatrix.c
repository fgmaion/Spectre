double unitTheory(double r, double k, double d, int order){
    // Evaluates the monopole/quadrupole correlation function (order={0,2}) corresponding to a unit theory vector, it P_l(k) = 1. for k_l<k<k_u, =0 otherwise. 

    switch(order){
        case 0:
            // include Taylor expansion if statement. 
            return  pow(pi, -2.)*(pow(r, -2.)*(k*sin(k*r)*sin(d*r) - d*cos(k*r)*cos(d*r)) + pow(r, -3.)*cos(k*r)*sin(d*r));

        /*
        case 2:
            // include Taylor expansion if statement.
            return -(1./(2.*pow(pi, 2.)))*( 3.*pow(r, -3.)*(gsl_sf_Si(r*ku) -gsl_sf_Si(r*kl)) - kl*pow(r, -2.)*cos(kl*r) + ku*pow(r, -2.)*cos(ku*r) - pow(r, -3.)*sin(ku*r) + pow(r, -3.)*sin(kl*r) - 3.*pow(r, -3.)*sin(ku*r) + 3.*pow(r, -3.)*sin(kl*r));
        */
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
    
    // for(ii=0; ii<n0*n1*n2; ii++)  printf("\n%e \t %e \t %e \t %e", rmodulus_vec[ii], Corrfn[ii], FFTW2_vecr_re[ii], FFTW2_vecr_im[ii]);
    
    int m0, m1, m2;
    
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
    
    //append a 'end of array' element to exit loop successfully. 
    rmodulus_vec[n0*n1*n2][0] =           9999.;
    rmodulus_vec[n0*n1*n2][1] = (double)   9999;
    
    printf("\n\nCalculating mixing matrix.\n");
    
    // ** Change to calculate quadrupole contribution.
    // for(transformOrder=0; transformOrder<2; transformOrder += 2){
        transformOrder = 0;
    
        for(nn=0; nn<(kBinNumb-1); nn++){         
            printf("%d \t %d \n", transformOrder, nn);
            
            // calculation of correlation function for given unit theory vector, 1 for the k interval, zero otherwise, either monopole or quadrupole. 
            for(k=0; k<n0*n1*n2; k++)       Corrfn[k]   = unitTheory(rmodulus_vec[k][0], kBinLimits[nn], 0.5*kbinInterval, transformOrder);
            
                                            Corrfn[0]   = 0.0; // **Include Taylor expansion in unit theory and get rid of this! ** nan from 1/0 will cover convolved Pk array after FFT.
            
            // Carry out the convolution by inverse FFT of the correlation fn. x FFT of W^2(k).  Renormalise to prevent underflow. 
            for(ii=0; ii<n0*n1*n2; ii++) out[ii][0]     = pow(10., 15.)*Corrfn[ii]*FFTW2_vecr_re[ii]; //*LegendrePolynomials(mu_vec[ii], transformOrder);
            for(ii=0; ii<n0*n1*n2; ii++) out[ii][1]     =               Corrfn[ii]*FFTW2_vecr_im[ii]; //*LegendrePolynomials(mu_vec[ii], transformOrder);
    
            // Calculate the convolved P(vec k).
            fftw_execute(iplan);
                
            // Set up for multipole calculation. 
            for(kk=0; kk<n0*n1*n2; kk++){
                polar2Dpk[kk][0]  =                  kmodulus_vec[kk];
	    	    polar2Dpk[kk][1]  =                        mu_vec[kk];
	    	    polar2Dpk[kk][2]  =     pow(n0*n1*n2, -1.0)*in[kk][0];
            }
    
            // Calculate monopole and quadrupole of convolved P(k) resulting from unit theory. 
            MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, modesPerBin, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0);
        
            for(mm=0; mm<(kBinNumb-1); mm++){
                mixingmatrix[mm][(transformOrder/2)*(kBinNumb-1) + nn]                =   kMonopole[mm];
           
                mixingmatrix[(kBinNumb-1) + mm][(transformOrder/2)*(kBinNumb-1) + nn] = kQuadrupole[mm];  
            }  
        }        
    // }     
    
    sprintf(filepath, "%s/Data/SpectralDistortion/mixingmatrix_VipersMask.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    // for(j=0; j<2*(kBinNumb-1); j++){
    //    for(k=0; k<2*(kBinNumb-1); k++){
    for(j=0; j<(kBinNumb-1); j++){
      for(k=0; k<(kBinNumb-1); k++){
        fprintf(output, "%e\n", mixingmatrix[j][k]);       
      }
    }
    
    fclose(output);
    
    sprintf(filepath, "%s/Data/SpectralDistortion/mixingmatrix_convolvedPk_VipersMask.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double Interim;
    
    for(j=0; j<(kBinNumb-1); j++){
        Interim = 0.0;
    
        for(k=0; k<(kBinNumb-1); k++)  Interim += mixingmatrix[j][k]*(*pt2Pk)(kBinLimits[k]);

        fprintf(output, "%e \t %e \n", kBinLimits[j], Interim);       
    }
    
    fclose(output);
    
    return 0;
}


int printf_unitTheoryxi(){
    sprintf(filepath, "%s/Data/SpectralDistortion/unitTheoryPk_monoxi_quadxi.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double dj, kj, Interim;
    
    for(j=0; j<1000; j++){ 
        dj  =   (double) j;
        dj /=         200.;
        
        dj -=           2.;
        
        dj  = pow(10., dj);
        
        fprintf(output, "%e \t ", dj);
    
        for(k=0; k<5; k++){  
            kj = (double) k;
                        
            kj -=        2.;
            
            fprintf(output, "%e \t ", TotalVolume*unitTheory(dj, pow(10., kj), 0.5*kbinInterval, 0));
        }
    
        fprintf(output, "\n");
    }
    
    fclose(output);

    return 0;
}


int printf_xi(){
    sprintf(filepath, "%s/Data/SpectralDistortion/unitTheory_xi_selectModes.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double dj, kj, Interim;
    
    // for(k=0; k<(kBinNumb-1); k++)  printf("%e \n", kBinLimits[k]);
    
    for(j=0; j<1000; j++){ 
        dj      =   (double) j;
        dj     /=         200.;
        
        dj     -=           2.;
        
        dj      = pow(10., dj);
        
        Interim = 0.0;
    
        // Cannot use kMonopole values due to aliased measurement. 
        for(k=0; k<(kBinNumb-1); k++) Interim += unitTheory(dj, kBinLimits[k], 0.5*kbinInterval, 0)*(*pt2Pk)(kBinLimits[k]);

        fprintf(output, "%e \t %e \n", dj, Interim);
    }
    
    fclose(output);
    
    return 0;
}
