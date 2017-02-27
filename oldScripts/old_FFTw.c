// Created 10/02/2017

#include "qSortCompare.c"

int PkCalc(){  
    fftw_execute(p);
    
    // H2_k assigned to hold fft of gals.
    for(j=0; j<n0*n1*n2; j++){
      H2_k[j][0] = H_k[j][0]; 
      H2_k[j][1] = H_k[j][1]; 
    }

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = surveyMask[j];

    fftw_execute(p);
    
    // free overdensity.
    // free_grid();

    assign2DPkMemory();
    
    PkCorrections();
    
    observedQuadrupole(polar_pkcount);
    
    return 0;
}


int PkCorrections(){
    double pk, WindowFunc;

    double rand_shot = 0.0, gal_shot = 0.0;

    // shot noise from randoms cat. calculation.    
    for(j=0; j<rand_number; j++){
        if(rand_accept[j]    == true)  rand_shot += pow(rand_weight[j], 2.);
    } 
                    
    rand_shot *= alpha*alpha;

    // shot noise from galaxies cat. calculation, including angular sampling. 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j]/sampling[j], 2.);
    }

    printf("\n\nShot noise contributions: Randoms %.4lf, Galaxies %.4lf", rand_shot, gal_shot);

    //** NOTE:  Clipping is assumed to not alter the galaxy shot noise 
    // as (a)   The volume affected is very small (~ 1 %) 
    //    (b)   Galaxies removed account for linear fluctuations to be 
    //          smaller, as opposed to a change in number density. 
    
    double newZero = H_k[0][0];

    polar_pkcount = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;
                
                if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
                if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
                if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;     

                WindowFunc                         = 1.;

		// 0/0 - > Taylor expand. 
                if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/xNyquistWaveNumber)/(pi*k_x*0.5/xNyquistWaveNumber);

                if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/yNyquistWaveNumber)/(pi*k_y*0.5/yNyquistWaveNumber);
    
                if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/zNyquistWaveNumber)/(pi*k_z*0.5/zNyquistWaveNumber);
                
                // Correct for mass assignment of randoms.
                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                H_k[Index][0]                     /= pow(WindowFunc, 2.);
                H_k[Index][1]                     /= pow(WindowFunc, 2.);

                // If smoothing the field, do not correct for mass assignment.  Obvious in hindsight.    
                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                H2_k[Index][0]                    /= pow(WindowFunc, 2.);
                H2_k[Index][1]                    /= pow(WindowFunc, 2.);
                
                                                   // gals          // normalised rands. 
                H_k[Index][0]                      = H2_k[Index][0] - alpha*H_k[Index][0]; 
		H_k[Index][1]                      = H2_k[Index][1] - alpha*H_k[Index][1]; 

		// enforce k=0 mode has P(0) = 0.
		// H_k[Index][0]                      = H2_k[Index][0] - H_k[Index][0]*(H2_k[0][0]/newZero);
		// H_k[Index][1]                      = H2_k[Index][1] - H_k[Index][1]*(H2_k[0][0]/newZero);
		
                pk                                 = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

                pk                                -= rand_shot;

		// fit for constant p(k) shotnoise of galaxies when clipping
                // pk                                -=  gal_shot;

		// printf("\n%.6lf \t %.6lf \t %.6lf %le", k_x, k_y, k_z, pk);
		
	            if(kmodulus > 0.000001){
	                //// Only half the modes are independent. ////
	            	if(k_z>0.){
			  // Restrict radial modes allowed into regresison calc. (hi-z clipping problem -> radial modes).
			  // if(k_z/kmodulus <= 0.8){
	            	    // One hemi-sphere is independent, e.g. k_z >= 0.
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk;		            	   
		            	            		            		            
		                twodim_pk[polar_pkcount][0]  = fabs(k_z);                       
	                        twodim_pk[polar_pkcount][1]  = pow(k_y*k_y + k_x*k_x, 0.5);    
                                twodim_pk[polar_pkcount][2]  = pk;
		            
		                polar_pkcount               += 1;
			    }
			  // }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]   = kmodulus;
		                polar_pk[polar_pkcount][1]   = fabs(mu);
		                polar_pk[polar_pkcount][2]   = pk;

		                twodim_pk[polar_pkcount][0]  = fabs(k_z);                       
	                        twodim_pk[polar_pkcount][1]  = pow(k_y*k_y + k_x*k_x, 0.5);    
                                twodim_pk[polar_pkcount][2]  = pk;
		            
		                polar_pkcount                += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
		                // on the line k_z=k_y=0, one half is independent, k_x>=0.
		                                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar_pk[polar_pkcount][0]    = kmodulus;
		                polar_pk[polar_pkcount][1]    = fabs(mu);
		                polar_pk[polar_pkcount][2]    = pk;
		                
		                twodim_pk[polar_pkcount][0]   = fabs(k_z);                       
	                        twodim_pk[polar_pkcount][1]   = pow(k_y*k_y + k_x*k_x, 0.5);    
                                twodim_pk[polar_pkcount][2]   = pk;
		            
		                polar_pkcount                += 1;
		            }
		            
		            // if((0.01923825<kmodulus) && (kmodulus<0.02407176))  printf("\n%e \t %e \t %e", kmodulus, mu, pk);
	            }
	        }
        }
    }
    
    // free_HOD();
    
    return 0;
}


int observedQuadrupole(int modeCount){
    // polar2DpkBinning(modeCount);

    // sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_spec_%d.dat", root_dir, loopCount);
    // sprintf(filepath,"%s/Data/500s/spoc_zobs_allgals/HOD_mock_512_specmask_%d.dat", root_dir, loopCount); 

    // sprintf(filepath,"%s/Data/500s/spec/HOD_mock_256_specmask_nbarweight_%d.dat", root_dir, loopCount);
    // sprintf(filepath,"%s/Data/500s/slowDFT/HOD_mocks_slowDFT_invnbar_weight_monteCarlo_%d.dat", root_dir, loopCount);
    
    // sprintf(filepath,"%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/Pk_ell/HOD_mock_256_subSample_localweighted_%d.dat", root_dir, loopCount);

    // sprintf(filepath,"%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/Pk_ell/HOD_mock_256_subSample_nonpoisson_Granettweighted_%d.dat", root_dir, loopCount);

    // sprintf(filepath,"%s/Data/500s/Nagoyav4_GaussianFields/nagoya_v4_%d.dat", root_dir, loopCount);
    
    // sprintf(filepath,"%s/Data/500s/spec/HOD_mock_256_specmask_clipped_%d.dat", root_dir, loopCount);
    
    // sprintf(filepath,"%s/Data/Jenkins_fold/HOD_cube_256_nofold_exprands_counts_sphere.dat", root_dir);
    
    // sprintf(filepath,"%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/HOD_cube_256_pk_galweights_depletionfactor_%.2lf_foldfactor_%.1lf_d0_%.1lf_R_%.1lf.dat", 
    //                                                                                                                             root_dir, depletion_factor, 
    //                                                                                                                             Jenkins_foldfactor, appliedClippingThreshold, clipping_smoothing_radius);
    
    // if(loopCount<10)        sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_00%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat", root_dir, loopCount);
    // else if(loopCount<100)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_0%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",  root_dir, loopCount);
    // else                    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/pk/clipped_d0_6/mocks_W1_v8.0_500_%d_Nagoya_v6_Samhain_specweight_fkpweight_multipoles.dat",   root_dir, loopCount);

    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/clipped_lnnormal_0.2_mock_001.dat", root_dir, loopCount);

    // sprintf(filepath,"%s/W1_Nagoya_v6_mocks_work/new_pk_15_06_15/clipped_d0_5/mock_%d_256_pk_d0_%.2lf.dat", root_dir, loopCount, appliedClippingThreshold);

    // sprintf(filepath,"%s/Data/500s/hod_cube/hod_cube_pk_filament_envrandoms_folded.dat", root_dir, fieldFlag, nz_smoothRadius);
  
  if(data_mock_flag == 0){
    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/true_nbar/declim_5.97/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
  
    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/stefano_comp/mock_%d_W%d_zlim_%.1lf_%.1lf_Jf_%d_permocknbar.dat", root_dir, loopCount, fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
  
    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/thesis_nbar_test/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/true_nbar/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    // sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/pk/parent_specz_err/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    // Clipping. 10 June 2016
    sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/pk/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    sprintf(filepath, "%s/John/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
  }

  else if(data_mock_flag == 1){
    // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/d0_%d/W%d/data_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) floor(appliedClippingThreshold), fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    // investigating rand_number dependence and rigorous angular limits of dodgy W1, hi z slice. 
    // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/d0_%d/W%d/data_morerand_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

    // commented 25th Jan. 
    // sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/pk/d0_%d/W%d/data_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) floor(appliedClippingThreshold), fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
    
    sprintf(filepath, "%s/John/data_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
  }

    /* COMMENTED 11 JAN FOR CUBE WORK. 
    if(data_mock_flag == 1){
        sprintf(filepath,"%s/W1_Spectro_V7_0/W%d_Spectro_v7_smoothed_%.1lf_nbar_pk.dat", root_dir, fieldFlag, nz_smoothRadius);
    }
    
    else if(data_mock_flag ==0){    
        // sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/pk_smoothed_%.1lf_reflected_2field_nbar/mock_%d_256_pk.dat", root_dir, fieldFlag, nz_smoothRadius, loopCount); 
        sprintf(filepath, "%s/W1_Spectro_V7_1/W%d_pk/mock_%d_d0_%.2lf_256_Jf_1_pk.dat", root_dir, fieldFlag, loopCount, appliedClippingThreshold); 
    }
    */
    // MonopoleCalc(kbin_no, mean_modk, kMonopole, polar_pk, modeCount, filepath, 0.0, 1.0, 1);

    MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, modeCount, filepath, 0.0, 1.0, 1);
    
    // HexadecapoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, kHexadecapole, polar2Dpk, modeCount, filepath, kbinInterval, 0.0, 1.0, 1);

    // if(data_mock_flag == 1) Cartesian2Dpk(modeCount);

    return 0;
}


int Cartesian2Dpk(int modeCount){    
    // for(j=0; j<modeCount; j++)  printf("\n%.3lf \t %.3lf \t %.3lf", twodim_pk[j][0], twodim_pk[j][1], twodim_pk[j][2]);    
    DualBinning(modeCount, twodim_pk, d2_binnedpk);

    printf("\n\nWriting 2D pk.");

    sprintf(filepath,"%s/W1_Spectro_V7_2/data_v1.7/pk_2d/pk_2d_W%d_%.1lf_%.1lf_Jf_%d.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
    
    output = fopen(filepath, "w");
    
    // corresponds to k=0.5 for **current** binning.
    for(k=0; k<kbin_no; k++){
        for(j=0; j<kbin_no; j++){
            fprintf(output, "%e \t", d2_binnedpk[k][j]);
        }
        
        fprintf(output, "\n");
    }
    
    fclose(output);
    
    return 0;
}


/*
int polar2DpkBinning(int modeCount){
    for(j=0; j<kBinNumb;  j++)        kBinLimits[j]  =                     kbinInterval*(j+1);
    for(j=0; j<muBinNumb; j++)       muBinLimits[j]  =   (1.0/(double) (muBinNumb - 1))*(j+1);

    DualBinning(modeCount, muBinNumb, muBinLimits, kBinNumb, kBinLimits, polar2Dpk, polar2DBinnedPk, mean_mu, mean_modk, polar_modesPerBin);

    sprintf(filepath, "%s/Data/SpectralDistortion/VipersMask_2Dpk_W2k.dat", root_dir);
    
    output = fopen(filepath, "w");

    printf("\nNumber of mu bins: %d", muBinNumb);

    for(j=0; j<muBinNumb-1; j++){
        for(k=0; k<kBinNumb-1; k++){
                    fprintf(output, "%e \t %e \t %e \t %e \n", mean_mu[j][k], mean_modk[j][k], TotalVolume*polar2DBinnedPk[j][k], (*pt2Pk)(mean_modk[j][k])*pow(1. + beta*pow(mean_mu[j][k], 2.), 2.)/(1. + 0.5*pow(mean_modk[j][k]*mean_mu[j][k]*velDispersion, 2.)));
        }
    }

    fclose(output);
    
    return 0;
}*/


int MonopoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    printf("\n\nPerforming multipole calculation. (to Monopole order)");
    
    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    loIndex = 0;
    hiIndex = 0;
    
    for(j=0; j<modBinNumb-1; j++){
        mean_modBin[j]    = 0.0;
        
        Monopole[j]       = 0.0;
    }
    
    for(j=0; j<modeCount; j++){
        if(Array[j][0]   >= logk_limits[0]){
            loIndex = j; 
            
            break;
        }
    }
    
    if(fileOutput==1)  output = fopen(filepath, "w");
        
    for(j=0; j<modBinNumb-1; j++){
        modes_perbin = 0;
    
        for(i=loIndex; i<modeCount; i++){
            if(Array[i][0] > logk_limits[j+1]){
                hiIndex = i;
                
                break;
            } 
        }
        
        for(i=loIndex; i<hiIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){      
                mean_modBin[j] += Array[i][0];
                   
                   Monopole[j] += Array[i][2];
                
                modes_perbin   += 1;
            }
        }
        	    
        if(modes_perbin != 0)  mean_modBin[j]  /= modes_perbin;
        if(modes_perbin != 0)     Monopole[j]  /= modes_perbin;

        // Peacock and Nicholson 1991, pg 313. above eqn (20).
        // Result of summing over a shell in k space containing m modes, should be a Gaussian random variable with variance 2.m/N^2  
        
        if(fileOutput==1)  fprintf(output, "%e \t %e \t %d \n", mean_modBin[j], Monopole[j], modes_perbin);   
        
        loIndex       = hiIndex;
    }
    
    if(fileOutput==1)  fclose(output);
    
    return 0;
}


int MultipoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    printf("\n\nPerforming multipole calculation. (to Quadrupole order)");
    
    // assign log k binning for pk and assign memory for binning. (logk_min, logk_max, # bins).
    // prep_pkbinning(-2., log10(modkMax), kbin_no);
    
    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);
    
    loIndex      = 0;
    hiIndex      = 0;

    for(k=0; k<modBinNumb-1; k++){
        mean_modBin[k]    = 0.0;
        Monopole[k]       = 0.0;
        Quadrupole[k]     = 0.0;
    }
    
    for(i=0; i<modeCount; i++){
        if(Array[i][0] >= logk_limits[0]){
            loIndex = i; 
            break;
        }
    }
    
    if(fileOutput==1)  output = fopen(filepath, "w");
    
    for(j=0; j<modBinNumb-1; j++){
        //  Linear regression against P_i(k, \mu_i)  
        //  L_i  = 0.5*(3.*\mu_i**2 -1.)
        //  P_i  = P_i(k, \mu_i)
    
        double Li          = 0.0;
        double Pi          = 0.0;
    
        double Sum_Li      = 0.0;
        double Sum_Li2     = 0.0;
    
        double Sum_Pi      = 0.0;
        double Sum_PiLi    = 0.0;
        
        modes_perbin = 0;
        
        // Find the range of indices corresponding to the modes in a given k interval. 
        for(i=loIndex; i<modeCount; i++){
            if(Array[i][0] > logk_limits[j+1]){
                hiIndex = i;
                break;
            } 
        }
        
        for(i=loIndex; i<hiIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){        
                mean_modBin[j] += Array[i][0];
                modes_perbin   += 1;
            
                         // L_i = 0.5*(3.*mu**2 -1.)
                     
                Li              = 0.5*(3.*pow(Array[i][1], 2.) - 1.);     
                Pi              = Array[i][2];
            
                // if(j==2)  printf("\n%e \t %e \t %e \t %e", Array[i][0], Array[i][1], Array[i][2], Li);
            
                Sum_Li         += Li; 
                Sum_Li2        += Li*Li;
            
                Sum_Pi         += Pi;
                Sum_PiLi       += Pi*Li;
            }
         }
    
         mean_modBin[j]        /= modes_perbin;
        	    
         // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
         // Required to invert AP  = B. 2x2 matrix inversion.
         	
         // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
         //            (c d)        

         // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li
     
         double detA;
         // det = ad - bc.
     
         detA                       = modes_perbin*Sum_Li2 - Sum_Li*Sum_Li;
     
         // if(detA<pow(10., -10.))      printf("\n\nCannot decompose into monopole and quadrupole, detA = 0.0");
     
         // if(j==2)  printf("\n\ndet A: %e", detA); 
     
         // (P_0, P_2)^T = (A^-1)B  = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
                     Monopole[j]    = (1./detA)*( Sum_Li2*Sum_Pi - Sum_Li*Sum_PiLi);
                   Quadrupole[j]    = (1./detA)*(-Sum_Li*Sum_Pi  + modes_perbin*Sum_PiLi);
                   
         // for(jj=0; jj<= j; jj++)  printf("\n%d \t %e \t %e \t %e \t %d", j, mean_modBin[jj], Monopole[jj], Quadrupole[jj], modesperbin[jj]);
         
         // printf("\n\n");
         
         
         //for(k=LowerBinIndex; k<UpperBinIndex; k++){
           //if((mu_lolimit<Array[k][1]) && (Array[k][1]<mu_hilimit)){
             //Li                       = 0.5*(3.*pow(Array[k][1], 2.) - 1.);
         
             //kMonopole_expError[j]   += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(Sum_Li2 - Li*Sum_Li), 2.);
           
             //kQuadrupole_expError[j] += pow((*pt2Pk)(meanKBin[j])*pow(1. + beta*pow(Array[k][1], 2.), 2.)*pow(1. + 0.5*pow(meanKBin[j]*velDispersion*Array[k][1], 2.), -1.)*(-1.*Sum_Li + + modesperbin[j]*Li), 2.);
           //}
         //}

           //kMonopole_expError[j]   *= pow(detA, -2.); 
         //kQuadrupole_expError[j]   *= pow(detA, -2.);
         
         loIndex   = hiIndex;

	 // printf("\n\n%s", filepath);
         
         if((fileOutput==1) && (detA>pow(10., -6.))){  
	   // printf("\n\n**Printing to file.**  %s", filepath);
	   
	   printf("\n%le \t %le \t %le \t %d", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);

	   fprintf(output, "%e \t %e \t %e \t %d \n", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);
	 }

	 // if(detA > pow(10., -6.))  printf("\n%le \t %le \t %le \t %d", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);
     } 
    
     if(fileOutput==1) fclose(output);
       
     // free logk_limits, mean_modk, binnedPk, modes_perbin
     // free_pkRegression();
     
     return 0;
}


int HexadecapoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double Hexadecapole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2 + P_4(k)*(35x^4 -30x^2 +3)/8.
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    printf("\n\nPerforming multipole calculation to Hexadecapole order.");
    
    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    // Order by mod k to ensure binning is the mean between LowerBinIndex and UpperBinIndex.
    qsort(Array, modeCount, sizeof(Array[0]), FirstColumnCompare);

    loIndex = 0;
    hiIndex = 0;
    
    for(k=0; k<modBinNumb-1; k++){
        mean_modBin[k]      = 0.0;
        Monopole[k]         = 0.0;
        Quadrupole[k]       = 0.0;
        Hexadecapole[k]     = 0.0;
    }
    
    for(i=0; i<modeCount; i++){
        if(Array[i][0] >= logk_limits[0]){
            loIndex = i; 
            
	    break;
        }
    }
    
    if(fileOutput==1)  output = fopen(filepath, "w");

    for(j=0; j<modBinNumb-1; j++){
        //  Linear regression against P_i(k, \mu_i)  
        //  L_i  = 0.5*(3.*\mu_i**2 -1.)
        //  P_i  = P_i(k, \mu_i)
        
        double Li          = 0.0;
        double Pi          = 0.0;
        double Ji          = 0.0;
    
        double Sum_Li      = 0.0;
        double Sum_Li2     = 0.0;
        double Sum_Ji      = 0.0;
        double Sum_Ji2     = 0.0;
    
        double Sum_Pi      = 0.0;
        double Sum_PiLi    = 0.0;
        double Sum_PiJi    = 0.0;
        double Sum_LiJi    = 0.0;
        
	modes_perbin       =   0;

        // Find the range of indices corresponding to the modes in a given k interval. 
        for(i=loIndex; i<modeCount; i++){
            if(Array[i][0] > logk_limits[j+1]){
                hiIndex = i;
                break;
            } 
        }
        
        for(i=loIndex; i<hiIndex; i++){
            if((mu_lolimit<Array[i][1]) && (Array[i][1]<mu_hilimit)){        
                mean_modBin[j]    += Array[i][0];
                modes_perbin      += 1;
                    
                Pi                 = Array[i][2];
                Li                 = 0.5*(3.*pow(Array[i][1], 2.) - 1.);    
                Ji                 = (35.*pow(Array[i][1], 4.) - 30.*pow(Array[i][1], 2.) + 3.)/8.; 
            
                Sum_Li            += Li; 
                Sum_Li2           += Li*Li;
                Sum_Ji            += Ji;
                Sum_Ji2           += Ji*Ji;
            
                Sum_Pi            += Pi;
                Sum_PiLi          += Pi*Li;
                Sum_PiJi          += Pi*Ji;
                Sum_LiJi          += Li*Ji;
            }
        }
    
        mean_modBin[j]            /= modes_perbin;
        	    
        // For a matrix A, paramater vector, (P_0, P_2, P_4)^T, P and vector B
        // Required to invert AP  = B. 3x3 matrix inversion.
         	
        // A reads as (a b c) for a = sum_Modes 1, b = sum L_i, c = sum J_i, d = sum_Modes L_i, e = sum_Modes Li**2, f = sum_Modes Ji*Li, g = sum_Modes Ji, h = sum_Modes Li*Ji, i = sum_Modes Ji*Ji 
        //            (d e f)
        //            (g h i)        

        // and B = (b1, b2, b3)^T for b1 = sum_Modes hat Pi, b2 = sum_Modes (hat Pi)*Li, b3 = sum_Modes (hat Pi)*Ji
     
        double detA;
        // det = ad - bc.
     
        detA = modes_perbin*(Sum_Li2*Sum_Ji2 - Sum_LiJi*Sum_LiJi) - Sum_Li*(Sum_Li*Sum_Ji2 - Sum_LiJi*Sum_Ji) + Sum_Ji*(Sum_Li*Sum_LiJi - Sum_Li2*Sum_Ji);
     
        // (P _0, P_2)^T = (A^-1)B = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
            Monopole[j]  = (1./detA)*(     (Sum_Li2*Sum_Ji2 - Sum_LiJi*Sum_LiJi)*Sum_Pi - (Sum_Li*Sum_Ji2          - Sum_Ji*Sum_LiJi)*Sum_PiLi + (Sum_Li*Sum_LiJi         - Sum_Ji*Sum_Li2)*Sum_PiJi);
          Quadrupole[j]  = (1./detA)*( -1.*(Sum_Li*Sum_Ji2  - Sum_LiJi*Sum_Ji  )*Sum_Pi + (modes_perbin*Sum_Ji2  - Sum_Ji*Sum_Ji  )*Sum_PiLi - (modes_perbin*Sum_LiJi - Sum_Ji*Sum_Li )*Sum_PiJi);
        Hexadecapole[j]  = (1./detA)*(     (Sum_Li*Sum_LiJi - Sum_Li2*Sum_Ji   )*Sum_Pi - (modes_perbin*Sum_LiJi - Sum_Li*Sum_Ji  )*Sum_PiLi + (modes_perbin*Sum_Li2  - Sum_Li*Sum_Li )*Sum_PiJi);
    
        loIndex = hiIndex;
 
	if((fileOutput == 1) && (detA>pow(10., -6.)))  fprintf(output, "%e \t %e \t %e \t %e \t %d \n", mean_modBin[j], Monopole[j], Quadrupole[j], Hexadecapole[j], modes_perbin);
    } 
    
    if(fileOutput == 1) fclose(output);
    
    return 0;
}


double invnbar_chisq(double chi){
    return chi*chi/interp_nz(chi);
}

double chisq(double chi){
    return chi*chi;
}

double chicubed(double chi){
    return chi*chi*chi;
}

double chicubed_nbar(double chi){
    return chi*chi*chi*interp_nz(chi);
}

double chisq_nbar(double chi){
    return chi*chi*interp_nz(chi);
}

double calc_volavg_invnbar(){
    // Shot noise correction for a varying background number density.  
    // P(k) is the volume avg. therefore so must be the shot noise correction. 

    double vol            =  qromb(&chisq, loChi, hiChi);
    
    double volavg_invnbar =  qromb(&invnbar_chisq, loChi, hiChi);
    
    // solid angle drops out if nbar is independent of direction. 
    return volavg_invnbar/vol;
}


double calc_volavg_chi(){
    double vol            =  qromb(&chisq, loChi, hiChi);
    
    double volavg_chi     =  qromb(&chicubed, loChi, hiChi);
    
    return volavg_chi/vol;
}


double calc_galavg_chi(){
    double vol            =  qromb(&chisq_nbar, loChi, hiChi);
    
    double volavg_chi     =  qromb(&chicubed_nbar, loChi, hiChi);
    
    return volavg_chi/vol;
}


double chiSq_fkpweight(double chi){
    return chi*chi*pow(interp_nz(chi)*fkpPk/(1. + interp_nz(chi)*fkpPk), 2.);
}


double calc_volavg_fkpweights(){
    double volavg_fkpweights     =  qromb(&chiSq_fkpweight, loChi, hiChi);
    
    return volavg_fkpweights;
}


int Gaussian_filter(double GaussianFilter_radius, int zero_mean){
    // Gaussian filter the array overdensity, filter radius set by global variable GaussianFilter_radius.
    prep_fftw();
    
    fftw_execute(p);
    
    // Gaussian smooth the counts. 
    double GaussianFilter;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;
                
                if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
                if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
                if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;      
       
                // Radius 1.         
                GaussianFilter                     = exp(-1.*kSq*0.5*pow(GaussianFilter_radius, 2.));
            
                H_k[Index][0]                     *= GaussianFilter;
                H_k[Index][1]                     *= GaussianFilter;
                
                H_k[Index][0]                     /= n0*n1*n2;
                H_k[Index][1]                     /= n0*n1*n2;
            }
        }
    }
    
    if(zero_mean == 1){
      // Return a zero mean field.
      H_k[0][0] = 0.0;
      H_k[0][1] = 0.0;
    }
    
    int negkIndex;
    
    // Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
    for(k=n0-1; k>=n0/2; k--){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                negkIndex          = k*n1*n2 + j*n2 + i;
                       
                Index              = 0;
                
                // zero maps to zero on reflection through the origin.
                if(i!=0)  Index   += (n2 - i);
                if(j!=0)  Index   += (n1 - j)*n2;
                          Index   += (n0 - k)*n1*n2;

                H_k[negkIndex][0]  =     H_k[Index][0];
                H_k[negkIndex][1]  = -1.*H_k[Index][1];
                
                if(negkIndex == Index)   H_k[Index][1] = 0.0; // purely real
            }
        }
    } 
    
    // in the k_z=0 plane one semi-circle is independent, k_y>0.         
    for(j=n1-1; j>=n1/2; j--){
        for(i=0; i<n2; i++){
            negkIndex          = j*n2 + i;
                       
            Index              = 0;
                
            // zero maps to zero on reflection through the origin.
            if(i!=0)  Index   += (n2 - i);
                      Index   += (n1 - j)*n2;

            H_k[negkIndex][0]  =     H_k[Index][0];
            H_k[negkIndex][1]  = -1.*H_k[Index][1];
            
            if(negkIndex == Index)   H_k[Index][1] = 0.0;
        }
    }
    
    // on the line k_z=k_y=0, one half is independent, k_x>=0.
    for(i=n2-1; i>=n2/2; i--){
        negkIndex          = i;
                       
        Index              = 0;
                
        // zero maps to zero on reflection through the origin.
        Index             += (n2 - i);

        H_k[negkIndex][0]  =      H_k[Index][0];
        H_k[negkIndex][1]  =  -1.*H_k[Index][1];
            
        if(negkIndex == Index)    H_k[Index][1] = 0.0;
    }

    iplan              = fftw_plan_dft_3d(n0, n1, n2, H_k, smooth_overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(iplan);

    // should be unnecessary.
    for(j=0; j<n0*n2*n1; j++)  smooth_overdensity[j][1] = 0.0;

    return 0;
}


int Gaussianfield(){
    int m0, m1, m2;
    
    double Power, amplitude, phase, expectation; 
    
    overdensity        = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
    H_k                = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(fftw_complex)); 
    
    iplan              = fftw_plan_dft_3d(n0, n1, n2, H_k, overdensity, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    
    bsigma8            =        1.0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                m0 = k;
                m1 = j;
                m2 = i;

                if(m2>n2/2)  m2                   -= n2;
                if(m1>n1/2)  m1                   -= n1;
                if(m0>n0/2)  m0                   -= n0;
                
                k_x                                = kIntervalx*m2;
                k_y                                = kIntervaly*m1;
                k_z                                = kIntervalz*m0;
                
                Index                              = k*n1*n2 + j*n2 + i;
                
                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);

                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;
                                                                                     // shot noise contribution. 
                // expectation                     = (*pt2Pk)(kmodulus)/TotalVolume; //  + (1./TotalVolume)*(*pt2shot)(1.);

                // Monopole and Quadrupole components, with z taken as the line of sight direction. 
                // expectation                     = (kaiser_multipole(kmodulus, beta, 0) + kaiser_multipole(kmodulus, beta, 2)*0.5*(3.*mu*mu -1.))*Pk_powerlaw(kmodulus, 5., 1.8)/TotalVolume;

                expectation                        = (*pt2Pk)(kmodulus)/TotalVolume;
                
                // expectation                    *= 1. + 0.5*pow(mu, 2.);
            
                // expectation                    *= pow(1. + 0.5*pow(mu, 2.), 2.);
                
                // fingers of God RSD.
                // expectation                    /= 1. + 0.5*pow(k_z*2., 2.);

                expectation                       *= (1. + LegendrePolynomials(mu, 2));

                // expectation                    *= spherical_tophat(kmodulus, 8.)*spherical_tophat(kmodulus, 8.)*pow(app_sigma8, -2.)*pow(   bsigma8,  2.);


                Power                              = -log(gsl_rng_uniform(gsl_ran_r))*expectation;
                
                // Note AMPLITUDE IS NOT the expectation. 
                // amplitude                          = sqrt(Power);
                amplitude                          = sqrt(expectation);
                                
                phase                              = 2.*pi*gsl_rng_uniform(gsl_ran_r);
                
                // Assuming Cic assignment scheme
                H_k[Index][0]                      = amplitude*cos(phase);
                H_k[Index][1]                      = amplitude*sin(phase);
                
                // WindowFunc                         = 1.;

                //if(k_x != 0.){
		          //  WindowFunc                    *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);}
                
                //if(k_y != 0.){
		          //  WindowFunc                    *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);}
                
                //if(k_z != 0.){
		          //  WindowFunc                    *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);}		      
	            
	            // H_k[Index][0]                  *= pow(WindowFunc, 2.);
	        	// H_k[Index][1]                  *= pow(WindowFunc, 2.);
	        }
        }
    }
    
    // Zero mean. Mean is always real for a real input fn.
    H_k[0][0] = 0.0;
    H_k[0][1] = 0.0;
    
    int negkIndex;
    
    // Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.
    for(k=n0-1; k>=n0/2; k--){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                negkIndex          = k*n1*n2 + j*n2 + i;
                       
                Index              = 0;
                
                // zero maps to zero on reflection through the origin.
                if(i!=0)  Index   += (n2 - i);
                if(j!=0)  Index   += (n1 - j)*n2;
                          Index   += (n0 - k)*n1*n2;

                H_k[negkIndex][0]  =     H_k[Index][0];
                H_k[negkIndex][1]  = -1.*H_k[Index][1];
                
                if(negkIndex == Index)   H_k[Index][1] = 0.0; // purely real
            }
        }
    } 
    
    // in the k_z=0 plane one semi-circle is independent, k_y>0.         
    for(j=n1-1; j>=n1/2; j--){
        for(i=0; i<n2; i++){
            negkIndex          = j*n2 + i;
                       
            Index              = 0;
                
            // zero maps to zero on reflection through the origin.
            if(i!=0)  Index   += (n2 - i);
                      Index   += (n1 - j)*n2;

            H_k[negkIndex][0]  =     H_k[Index][0];
            H_k[negkIndex][1]  = -1.*H_k[Index][1];
            
            if(negkIndex == Index)   H_k[Index][1] = 0.0;
        }
    }
    
    // on the line k_z=k_y=0, one half is independent, k_x>=0.
    for(i=n2-1; i>=n2/2; i--){
        negkIndex          = i;
                       
        Index              = 0;
                
        // zero maps to zero on reflection through the origin.
        Index             += (n2 - i);

        H_k[negkIndex][0]  =      H_k[Index][0];
        H_k[negkIndex][1]  =  -1.*H_k[Index][1];
            
        if(negkIndex == Index)    H_k[Index][1] = 0.0;
    }

    fftw_execute(iplan);

    // should be unnecessary.
    for(j=0; j<n0*n2*n1; j++)  overdensity[j][1] = 0.0;
    
    fftw_free(H_k);
    
    fftw_destroy_plan(iplan);
    
    
    sprintf(filepath,"%s/Data/SpectralDistortion/GRF_mask_MonoAndQuad_CellSize_2.00_papercheck.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(k=0; k<n0*n1*n2; k++)  fprintf(output, "%le \n", 1.);  // overdensity[k][0]

    fclose(output);
    
    return 0;
}


int lnNormfield(){
    Gaussianfield();
            
    double GRF_var = 0.0;
            
    // Already a zero mean field. 
    // for(j=0; j<n0*n1*n2; j++) GRF_var         += pow(overdensity[j][0], 2.);

    // GRF_var                                   /= n0*n1*n2;

    GRF_var = 1.011;

    for(j=0; j<n0*n2*n1; j++)  overdensity[j][0] = exp(overdensity[j][0] - 0.5*GRF_var) - 1.;

    return 0; 
}

/*
int printGaussianfield(){
    sprintf(filepath, "%s/Data/SpectralDistortion/GRF_mask_MonoAndQuad_CellSize_%.2f.dat", root_dir, CellSize);

    output = fopen(filepath, "w");
    
    for(j=0; j<n0*n1*n2; j++)  fprintf(output, "%e \n", densityArray[j]);
    
    fclose(output);

    return 0;
}*/

/*
int BinnedPkForMuInterval(double lowerMuLimit, double upperMuLimit, char filepath[], int modeCount){
    printf("\nPerforming P(k) binning for given mu interval, %f < mu < %f", lowerMuLimit, upperMuLimit);

    int muInterval_modeCount = 0; 

    for(k=0; k<modeCount; k++){
        if((lowerMuLimit<polar2Dpk[k][1]) && (polar2Dpk[k][1] < upperMuLimit)){
            muIntervalPk[k][0]    = polar2Dpk[k][0];
            muIntervalPk[k][1]    = polar2Dpk[k][2];  
        
            muInterval_modeCount += 1;
        }
    }
    
    printf("\nNumber of modes in interval:  %d", muInterval_modeCount);
    
    MonopoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polar_pkcount, filepath, kbinInterval, 0.0, 1.0, 1);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++)     fprintf(output, "%e \t %e \n", meanKBin[j], TotalVolume*binnedPk[j]);
    
    fclose(output);
    
    return 0;
}*/



int DualBinning(int NumberModes, double** DualParamArray, double** BinnedDualParamArray){
    int m;

    int bin_no = 50;

    double firstBinLimits[bin_no];
    double secndBinLimits[bin_no];

    int modesPerBin[bin_no][bin_no];

    int firstColLowerBinIndex     = 0;
    int firstColUpperBinIndex     = 0;

    int secndColLowerBinIndex     = 0;
    int secndColUpperBinIndex     = 0;
    

    qsort(DualParamArray, NumberModes, sizeof(DualParamArray[0]), FirstColumnCompare);

    // Order by first then second column.
    printf("\nDual param array sorted.");
    
    for(j=0; j<bin_no; j++) firstBinLimits[j] = 0.01 + j*1./bin_no;
    for(j=0; j<bin_no; j++) secndBinLimits[j] = 0.01 + j*1./bin_no;


    // for(j=0; j<20; j++)  printf("\n%.3lf \t %.3lf \t %.3lf", DualParamArray[j][0], DualParamArray[j][1], DualParamArray[j][2]);
    
    for(j=0; j<bin_no; j++){
        for(k=0; k<bin_no; k++){
            BinnedDualParamArray[j][k] = 0.0;
            // mean_firstCol[j][k]        = 0.0;
            // mean_secndCol[j][k]        = 0.0;
            modesPerBin[j][k]          =   0;    
        }
    }
    
    for(j=0; j<NumberModes; j++){
        if(DualParamArray[j][0]   >= firstBinLimits[0]){
            firstColLowerBinIndex = j; 
        
            break;
        }
    }
    
    for(j=0; j<bin_no; j++){
        for(i=firstColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][0] > firstBinLimits[j+1]){
                firstColUpperBinIndex = i;
                break;
            } 
        }
        
        secndColLowerBinIndex = firstColLowerBinIndex;
    
	    qsort(&DualParamArray[firstColLowerBinIndex], firstColUpperBinIndex - firstColLowerBinIndex, sizeof(DualParamArray[0]), SecondColumnCompare);
    
        // meets the lower bin limit
        for(i=secndColLowerBinIndex; i<NumberModes; i++){
            if(DualParamArray[i][1]   >= secndBinLimits[0]){
                secndColLowerBinIndex = i; 
        
                break;
            }
        }      
       
        // meets the upper bin limit. 
        for(k=0; k<bin_no; k++){
            for(m=secndColLowerBinIndex; m<firstColUpperBinIndex; m++){
                if(DualParamArray[m][1] > secndBinLimits[k+1]){
                    secndColUpperBinIndex = m;
                    break;
                } 
            }
        
            for(m=secndColLowerBinIndex; m<secndColUpperBinIndex; m++){
                BinnedDualParamArray[j][k]    += DualParamArray[m][2];
                // mean_firstCol[j][k]           += DualParamArray[m][0];
                // mean_secndCol[j][k]           += DualParamArray[m][1];
                modesPerBin[j][k]             += 1;    
            }
            
            if(modesPerBin[j][k] != 0)  BinnedDualParamArray[j][k] /= modesPerBin[j][k];
            // if(modesPerBin[j][k] != 0)  mean_firstCol[j][k]        /= modesPerBin[j][k];
            // if(modesPerBin[j][k] != 0)  mean_secndCol[j][k]        /= modesPerBin[j][k];
            
            secndColLowerBinIndex = secndColUpperBinIndex;
        }
        
        firstColLowerBinIndex = firstColUpperBinIndex;
    }
    
    return 0;
}
