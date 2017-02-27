int Gaussianfilter(double array[], int len, double Interval, double sigma){
    int        newlen;
    int    kernelSize;
  
    double        Sum;
    double Norm = 0.0;

    double* NewArray;

    kernelSize  = (int) ceil(4.*sigma/Interval); // (int) ceil(len/8);

    newlen      = len + 2*kernelSize;


    NewArray    =  (double *) malloc(newlen*sizeof(*NewArray));


    for(j=0; j<newlen; j++) NewArray[j] = 0.0;

    // padded array.                                                                                                                                          
    for(j=kernelSize; j<newlen - kernelSize; j++)   NewArray[j] = array[j-kernelSize];

    for(i=-kernelSize; i<kernelSize+1; i++)    Norm += exp(-0.5*pow(i*Interval/sigma, 2.0));

    for(j=kernelSize; j<newlen - kernelSize; j++){
        Sum    = 0.0;

        for(i=-kernelSize; i<kernelSize + 1; i++){
            Sum +=  pow(Norm, -1.)*exp(-0.5*pow(i*Interval/sigma, 2.0))*NewArray[j + i];
        }

        array[j - kernelSize] = Sum;
    }

    return 0;
}


int test_Gaussianfilter(){
  double test[100]; 

  prep_nbar();

  printf("\n\nChecking divisible by 8, len(chibins): %d", chibin_no);

  for(j=0; j<100; j++)  test[j] = 0.0;

  // Dirac delta. 
  test[50] = 1.0;

  Gaussianfilter(test, 100, 1., 10.);

  sprintf(filepath, "%s/W1_Spectro_V7_2/Gaussian_filterCheck.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<100; j++)  fprintf(output, "\n%.4lf \t %.4lf", j*1., test[j]);

  return 0;
}


int reflect_up(double maxr, int j, int Index){
    double dist, reflected_chi, vol_factor; 
    
    int rIndex;
    
    dist            = maxr - rDist[j];

    reflected_chi   = rDist[j] + 2.*dist;
    
    vol_factor      = pow(reflected_chi/rDist[j], 2.);

    Index           = (int) floor(reflected_chi/chi_interval);
            
    // chibins[Index] += vol_factor*reflected_chi/sampling[j];
            
    Nchi[Index]    += vol_factor*clip_galweight[j]/sampling[j]; 

    return 0;
}

int reflect_down(double minr, int j, int Index){
    double dist, reflected_chi, vol_factor; 
    
    int rIndex;

    dist            = rDist[j] - minr;

    reflected_chi   = rDist[j] - 2.*dist;

    // Expected number of galaxies scales with volume. 
    vol_factor      = pow(reflected_chi/rDist[j], 2.);

    rIndex          = (int) floor(reflected_chi/chi_interval);
            
    // chibins[rIndex] += vol_factor*reflected_chi/sampling[j];
            
    Nchi[rIndex]    += vol_factor*clip_galweight[j]/sampling[j]; 

    return 0;
}


int smoothed_nbar_calc(double kernel_width, int reflect){    
    // equal intervals in comoving distance, for both W1 and W4 fields.  
    prep_nbar();
    
    double  chi, norm = 0.0;

    double minr;
    double maxr;

    int global_fieldFlag;

    // Store initial fieldFlag. 
    global_fieldFlag = fieldFlag;
    
    // Two field average.
    for(ii=1; ii<5; ii=ii+10){ // fix to 3 for joint-field
      fieldFlag = ii;

      // analysis on **mock** catalogues.                                                                            
      if(data_mock_flag == 0){
	    // Nagoya v6. (v7.) mocks //                                                                      
	    // if(loopCount<10)       sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_00%d_spec_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, loopCount);
	    // else if(loopCount<100) sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_0%d_spec_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
	    // else                   sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_%d_spec_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);

            if(loopCount<10)        sprintf(filepath, "%s/mocks_v1.7/W%d/mock_00%d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
	    else if(loopCount<100)  sprintf(filepath, "%s/mocks_v1.7/W%d/mock_0%d_VAC_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);
	    else                    sprintf(filepath, "%s/mocks_v1.7/W%d/mock_%d_VAC_Nagoya_v6_Samhain.dat",    vipersHOD_dir, fieldFlag, loopCount);
    
	    CatalogueInput_500s(filepath);
      
  	    // Choice of redshift from zcos, zpec, zphot, zobs.                                               
	    gal_z = &zobs[0];

	    // load sampling according to local TSR.                                                                       
	    spec_weights();

	    // clipping weights.
	    if(appliedClippingThreshold < 1000.)  load_clippingweights();

            // nbar over full redshift range, no z cut. 
            assignAcceptance();
      
	    // cut at z boundaries, with reflection. 
	    // if(reflect==1) assignAcceptance();
      }
      
      // analysis on **data** catalogues.                                                                                   
      if(data_mock_flag == 1){
	    // Also assigns sampling: ESR = TSR x SSR.
  	    DataInput();

	    // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                                                                                   
            gal_z = &zobs[0];

	    // do not impose redshift cuts just yet, obtain full n(z). Some redshifts >2, remove these objects.
	    for(j=0; j<Vipers_Num; j++){  
	      if((0.1<zobs[j]) && (zobs[j]<1.7)) Acceptanceflag[j] = true;	    
	      
	      else{Acceptanceflag[j] = false;}
	    }
      }
      
      for(j=0; j<Vipers_Num; j++)  rDist[j] = interp_comovingDistance(zobs[j]);

      minr = AcceptedMin(rDist, Acceptanceflag, Vipers_Num);
      maxr = AcceptedMax(rDist, Acceptanceflag, Vipers_Num);
      
      printf("\nStationary distances accepted: %le \t %le", minr, maxr);

      for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
            Index             = (int) floor(rDist[j]/chi_interval);
            
            Nchi[Index]      +=  clip_galweight[j]/sampling[j]; 
            
	    if((loChi <= rDist[j]) && (rDist[j] <= hiChi))  norm +=  clip_galweight[j]/sampling[j];

            if((reflect == 1) &&  ((rDist[j] - minr) < kernel_width))  reflect_down(minr, j, Index);
            if((reflect == 1) &&  ((maxr - rDist[j]) < kernel_width))    reflect_up(maxr, j, Index);
	}  
      }
    }

    // Filter counts. 
    Gaussianfilter(Nchi, chibin_no, chi_interval, kernel_width);

    for(j=0; j<chibin_no; j++)  comovVol[j] = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]  = Nchi[j]/comovVol[j]; 

    // or filter <n>
    // Gaussianfilter(nbar, chibin_no, chi_interval, kernel_width);
    
    // Remove reflected galaxies. 
    for(j=0; j<chibin_no; j++){
      if((reflect == 1) && (chibins[j] > maxr) || (chibins[j] < minr)){  
	nbar[j] = 0.0;
	
	Nchi[j] = 0.0;
      }
    }

    // assuming the surveyed volume is a fair sample -> \sum_g 1/ESR in the volume should be \int <n> d3x. Rescale <n> such that this is the case,                                                                                          
    // this imposes the integral constraint. Calculate renormalisation, \int <n> dV = N_W1 + N_W4, ESR corrected.  1% correction. 
    double result = 0.0;

    for(j=0; j<chibin_no; j++){
      // redshift limits.                                                                                                                                                                                                                   
      if((chibins[j] >= loChi) && (chibins[j] <= hiChi)){
        result += pow(chibins[j], 2.)*chi_interval*nbar[j];
      }
    }

    result *= sqdegs2steradians(TotalW1W4area);

    printf("\n\n<n> normalisation ratio: %.4lf", norm/result);

    // renormalise.                                                                                                                                                                                                                         
    for(j=0; j<chibin_no; j++)     nbar[j] *= norm/result;
    
    // if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/smoothed_nbar/nbar_smoothed_reflected_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg.dat", root_dir, kernel_width, loopCount);
    // if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_0/nbar_smoothed_reflected_%.1lf_%.2lf_%.2lf_Nagoya_v7_WX_SPECTRO_V7_twofield_avg.dat", root_dir, kernel_width, lo_MBlim, hi_MBlim);
    
    // if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg_v2.dat", root_dir, kernel_width, loopCount);
    // if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_twofield_avg.dat", root_dir, kernel_width);
    
    if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg_clippedWeights.dat", root_dir, kernel_width, loopCount);
    // if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_twofield_avg.dat", root_dir, kernel_width);

    // printf("\n\n%s", filepath);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)  fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);

    // Restore global field Flag.
    fieldFlag = global_fieldFlag;
    
    return 0;
}
