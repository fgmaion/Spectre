#include "qSortCompare.c"

int PkCalc(){
    fftw_execute(plan);
    /*
    oldprep_pkRegression();
    
    // correct_ind_modes();
    correct_all_modes();

    oldnosort_MultipoleCalc(kbin_no, mean_modk, Monopole, Quadrupole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 0);
    */
    // Correct_modes();
    
    // observedQuadrupole();

    nosort_MultipoleCalc();
    
    // MultipoleCalc(kbin_no, mean_modk, Monopole, Quadrupole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 0);
    
    return 0;
}


int prep_r2c_modes(){
  int dummy;
  
  // r2c returns half the modes on the direction in which overdensity changes first, i.e. x. 
  // #pragma omp parallel for private(Index, k, j, i, k_x, k_y, k_z, pk, kmodulus, mu)
  for(k=0; k<n0; k++){
    k_z = kIntervalz*k;

    if(k_z > zNyquistWaveNumber)  k_z   -= n0*kIntervalz;
    
    for(j=0; j<n1; j++){
      k_y = kIntervaly*j;

      if(k_y > yNyquistWaveNumber)  k_y   -= n1*kIntervaly;

      // for(i=0; i<(n2/2 + 1); i++){
      for(i=0; i<n2; i++){ // limit is (n2/2 + 1)*n1*n0 for r2c.
        k_x = kIntervalx*i;

        if(k_x > xNyquistWaveNumber)  k_x   -= n2*kIntervalx; //  Remove for r2c. int rather than double condition.  
        
        // Index                                  = k*n1*(n2/2+1) + j*(n2/2+1) + i;
        Index                                  = k*n1*n2 + j*n2 + i;
        
        kSq                                    = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

        kmodulus                               = pow(kSq, 0.5);

        mu                                     = k_z/kmodulus; // Assume mean is not passed.

        kLi[Index]                             = gsl_sf_legendre_P2(mu);  // L_2 = 0.5*(3.*mu**2 -1.) 

        kM2[Index]                             = gsl_sf_sinc(k_x*0.5/xNyquistWaveNumber); // Computes \sinc(x) = \sin(\pi x) / (\pi x).
        kM2[Index]                            *= gsl_sf_sinc(k_y*0.5/yNyquistWaveNumber);
        kM2[Index]                            *= gsl_sf_sinc(k_z*0.5/zNyquistWaveNumber);

        kM2[Index]                             = pow(kM2[Index], 2.0);  // Correct mass assignment of randoms; cic = 2, ngp = 1.

        dummy                                  = (int)  floor((log10(kmodulus) - logk_min)/logk_interval);
        
        // Needs properly fixed.  Discarding zero index info.
        //                                     if                                            then    else
        kind[Index]                            = ((dummy >= 0) && (dummy < kbin_no)) ? dummy : 0;
        
        // Each available mode has an index in the binning array. 
        Sum_Li[kind[Index]]                   += kLi[Index];
        Sum_Li2[kind[Index]]                  += kLi[Index]*kLi[Index];

        modes_perbin[kind[Index]]             += 1;

        mean_modk[kind[Index]]                += kmodulus;

        // printf("\n%d \t %.4lf \t %.4lf \t %d", Index, kLi[Index], kM2[Index], kind[Index]);
      }
    }
  }

  int sum_modes = 0;
  
  for(j=0; j<kbin_no; j++){
    mean_modk[j]  /= modes_perbin[j];

    detA[j]        = modes_perbin[j]*Sum_Li2[j] - Sum_Li[j]*Sum_Li[j];

    // printf("\n%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf", modes_perbin[j], mean_modk[j], Sum_Li[j], Sum_Li2[j], detA[j]);

    sum_modes     += modes_perbin[j];
  }

  printf("\n%d \t %d", sum_modes, n0*n1*n2);
  
  return 0;
}


int correct_all_modes(){
  int start;
  double pk;
  double rand_shot = 0.0, gal_shot = 0.0;

  for(j=0; j<rand_number; j++)  rand_shot += pow(rand_weight[j], 2.);  // shot noise from randoms.

  rand_shot *= alpha*alpha;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j]/sampling[j], 2.);  // galaxy shot noise, inc. sampling.
  }

  printf("\n\nShot noise: randoms %.4lf, galaxies %.4lf", rand_shot, gal_shot);
  
  for(k=0; k<n0; k++){  // k_z > 0.0; no restriction on k_x or k_y.
    k_z = kIntervalz*k;

    if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;
    
    for(j=0; j<n1; j++){
      k_y = kIntervaly*j;

      if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
      
      for(i=0; i<n2; i++){
        k_x = kIntervalx*i;

        if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;  // Way to remove this?
      	
        Index                              = k*n1*n2 + j*n2 + i;

        PkCorrections(Index, k_x, k_y, k_z, rand_shot, gal_shot, &pk, &kmodulus, &mu);

        polar_pk[Index][0] = kmodulus;   
        polar_pk[Index][1] = fabs(mu);   
        polar_pk[Index][2] = pk;
      }
    }
  }

  polar_pkcount = n0*n1*n2;
  
  return 0;
}


int correct_ind_modes(){
  int start;
  double pk;
  double rand_shot = 0.0, gal_shot = 0.0;

  for(j=0; j<rand_number; j++)  rand_shot += pow(rand_weight[j], 2.);  // shot noise from randoms.

  rand_shot *= alpha*alpha;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j]/sampling[j], 2.);  // galaxy shot noise, inc. sampling.
  }

  printf("\n\nShot noise: randoms %.4lf, galaxies %.4lf", rand_shot, gal_shot);

  //** NOTE:  Clipping is assumed to not alter the galaxy shot noise
  // as (a)   The volume affected is very small (~ 1 %)
  //    (b)   Galaxies removed account for linear fluctuations to be
  //          smaller, as opposed to a change in number density.
  
  for(i=1; i<n2/2+1; i++){
    k_x = kIntervalx*i;  // k_z = 0.0 && k_y = 0.0 && k_x > 0.0; drop the mean (k=0) mode. 

    PkCorrections(i, k_x, 0.0, 0.0, rand_shot, gal_shot, &pk, &kmodulus, &mu);
    
    polar_pk[i-1][0] = kmodulus;         // One hemi-sphere is independent, e.g. i) k_z >= 0.0
    polar_pk[i-1][1] = fabs(mu);         // ii) in the k_z=0 plane one semi-circle is independent, k_y>0.
    polar_pk[i-1][2] =       pk;         // iii) on the line k_z=k_y=0, one half is independent, k_x>=0.
  }
  
  start = -n2/2;
  
  for(j=1; j<n1/2+1; j++){  // k_z = 0.0; k_y > 0.0;
    k_y = kIntervaly*j;
    for(i=0; i<n2; i++){
      k_x = kIntervalx*i;

      if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;

      Index                              = j*n2 + i;   // k=0 by assumption;   

      PkCorrections(Index, k_x, k_y, 0.0, rand_shot, gal_shot, &pk, &kmodulus, &mu);

      // Starts at n2/2 from first loop, but k_y starts at 1. 
      polar_pk[start + Index][0]         = kmodulus;   // One hemi-sphere is independent, e.g. i) k_z >= 0.0
      polar_pk[start + Index][1]         = fabs(mu);   // ii) in the k_z=0 plane one semi-circle is independent, k_y>0.
      polar_pk[start + Index][2]         = pk;         // iii) on the line k_z=k_y=0, one half is independent, k_x>=0.
    }
  }
  
  start = -n1*n2/2 + n2/2;
  
  // #pragma omp parallel for private(Index, k, j, i, k_x, k_y, k_z, pk, kmodulus, mu)
  for(k=1; k<=n0/2; k++){  // k_z > 0.0; no restriction on k_x or k_y. 
     k_z = kIntervalz*k;

     for(j=0; j<n1; j++){
       k_y = kIntervaly*j;

       for(i=0; i<n2; i++){
   	     k_x = kIntervalx*i;

         if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;  // Way to remove this?
         if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;

         Index                              = k*n1*n2 + j*n2 + i;

         PkCorrections(Index, k_x, k_y, k_z, rand_shot, gal_shot, &pk, &kmodulus, &mu);
	 
         polar_pk[start + Index][0] = kmodulus;   // One hemi-sphere is independent, e.g. i) k_z >= 0.0
         polar_pk[start + Index][1] = fabs(mu);   // ii) in the k_z=0 plane one semi-circle is independent, k_y>0.
         polar_pk[start + Index][2] = pk;         // iii) on the line k_z=k_y=0, one half is independent, k_x>=0.
	                                              // in the k_z=0 plane one semi-circle is independent, k_y>0.
       }
     }
  }
  
  polar_pkcount = n0*n1*n2/2 + n1*n2/2 + n2/2;  // Number of modes, Index is (polar_pkcount -1).
  
  return 0;
}


int observedQuadrupole(){
  if(data_mock_flag == 0){
    // sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/pk/d0_%d/W%d/mock_%d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, (int) appliedClippingThreshold, fieldFlag, loopCount, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);  // Clipping. 10 June 2016

    sprintf(filepath, "%s/W1_Spectro_V7_4/pk.dat", root_dir); 
  }

  else if(data_mock_flag == 1){
    sprintf(filepath, "%s/John/data_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);
  }

  nosort_MultipoleCalc(filepath, 0);
  
  return 0;
}


int MultipoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
    // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
    // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

    int          loIndex;
    int          hiIndex;
    int     modes_perbin;

    printf("\n\nPerforming multipole calculation. (to Quadrupole order)");
    
    // assign log(k) binning for P(k) and assign memory for binning. (logk_min, logk_max, # bins).
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
        
        for(i=loIndex; i<modeCount; i++){
            if(Array[i][0] > logk_limits[j+1]){
              hiIndex = i;  // Find the range of indices corresponding to the modes in a given k interval.

              break;
            } 
        }
        
        for(i=loIndex; i<hiIndex; i++){
             mean_modBin[j] += Array[i][0];
             modes_perbin   += 1;
            
             Li              = 0.5*(3.*pow(Array[i][1], 2.) - 1.);  // L_i = 0.5*(3.*mu**2 -1.)     
             Pi              = Array[i][2];
            
             Sum_Li         += Li; 
             Sum_Li2        += Li*Li;
            
             Sum_Pi         += Pi;
             Sum_PiLi       += Pi*Li;
         }
    
         mean_modBin[j]        /= modes_perbin;
        	    
         // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
         // Required to invert AP  = B. 2x2 matrix inversion.
         	
         // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
         //            (c d)        

         // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li
     
         double detA;  // det = ad - bc.
     
         detA                       = modes_perbin*Sum_Li2 - Sum_Li*Sum_Li;
     
           Monopole[j]              = (1./detA)*( Sum_Li2*Sum_Pi - Sum_Li*Sum_PiLi);
         Quadrupole[j]              = (1./detA)*(-Sum_Li*Sum_Pi  + modes_perbin*Sum_PiLi);
                            
         loIndex   = hiIndex;

         if(detA>pow(10., -6.))  printf("\n%le \t %le \t %le \t %d", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);
         
         if((fileOutput==1) && (detA>pow(10., -6.)))  fprintf(output, "%e \t %e \t %e \t %d \n", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin);
     } 
    
     if(fileOutput==1) fclose(output);
       
     // free logk_limits, mean_modk, binnedPk, modes_perbin
     // free_pkRegression();
     
     return 0;
}


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


int Cartesian2Dpk(int modeCount){
  // for(j=0; j<modeCount; j++)  printf("\n%.3lf \t %.3lf \t %.3lf", twodim_pk[j][0], twodim_pk[j][1], twodim_pk[j][2]);
  // DualBinning(modeCount, twodim_pk, d2_binnedpk);

  printf("\n\nWriting 2D pk.");

  sprintf(filepath,"%s/W1_Spectro_V7_2/data_v1.7/pk_2d/pk_2d_W%d_%.1lf_%.1lf_Jf_%d.dat", root_dir, fieldFlag, lo_zlim, hi_zlim, (int) Jenkins_foldfactor);

  output = fopen(filepath, "w");

  // corresponds to k=0.5 for **current** binning.
  for(k=0; k<kbin_no; k++){
    for(j=0; j<kbin_no; j++){
      // fprintf(output, "%e \t", d2_binnedpk[k][j]);
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
        fprintf(output, "%e \t %e \t %e \t %e \n", mean_mu[j][k], mean_modk[j][k], TotalVolume*polar2DBinnedPk[j][k], (*pt2Pk)(mean_modk[j][k])*pow(1. + beta*pow(mean_mu[j][k], \
    2.), 2.)/(1. + 0.5*pow(mean_modk[j][k]*mean_mu[j][k]*velDispersion, 2.)));                                                                                                               
      }                                                                                                                                                                                     
    }                                                                                                                                                                                                                                                                                                                                                                                     
    fclose(output);                                                                                                                                                                          
                                                                                                                                                                                             
    return 0;                                                                                                                                                                                 
    }
*/
