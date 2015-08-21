int Gaussianfilter(double array[], int len, double Interval, double sigma){
    int       newlen;
    int   kernelSize;
  
    double        Sum;
    double Norm = 0.0;

    double* NewArray;

    kernelSize = (int) ceil(len/8);

    newlen     = len + 2*kernelSize;


    NewArray   =  (double *) malloc(newlen*sizeof(*NewArray));


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

/*
int NormedGaussianfilter(float (*transformArray)(float, int), float (*inversetransform)(float, int), float array[], float xarray[], int len, double Interval, float sigma, double externNorm, float output[]){

    float Interim[len];
  
    double InterimNorm = 0.0;
    
    for(j=0; j<len; j++)  Interim[j]   = (*transformArray)(array[j], j);
    
    Gaussianfilter(Interim, xarray, len,  Interval, sigma);
    
    for(j=0; j<len; j++)  Interim[j]   = (*inversetransform)(Interim[j], j);
      
    for(j=0; j<len; j++)  InterimNorm += (double) Interim[j];
    
    for(j=0; j<len; j++)  output[j]    = (externNorm/InterimNorm)*Interim[j];

    return 0;
}*/


int reflect_up(double chi, int j, int Index){
    double dist, reflected_chi, vol_factor; 
    
    int rIndex;
    
    dist            = hiChi - chi;

    reflected_chi   = chi + 2.*dist;
    
    vol_factor      = pow(chi/reflected_chi, 2.);

    Index           = (int) floor(reflected_chi/chi_interval);
            
    chibins[Index] += vol_factor*reflected_chi/sampling[j];
            
    Nchi[Index]    += vol_factor*1./sampling[j]; 

    return 0;
}


int reflect_down(double chi, int j, int Index){
    double dist, reflected_chi, vol_factor; 
    
    int rIndex;
    
    dist            = chi - loChi;

    reflected_chi   = chi - 2.*dist;

    vol_factor      = pow(reflected_chi/chi, 2.);

    rIndex          = (int) floor(reflected_chi/chi_interval);
            
    chibins[rIndex] += vol_factor*reflected_chi/sampling[j];
            
    Nchi[rIndex]    += vol_factor*1./sampling[j]; 

    return 0;
}


int smoothed_nbar_calc(double kernel_width, int reflect){    
    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    prep_nbar();
    
    double chi;

    int global_fieldFlag;

    // Store initial fieldFlag. 
    global_fieldFlag = fieldFlag;

    // Two field average.
    for(ii=1; ii<5; ii=ii+3){
      fieldFlag = ii;

      // analysis on **mock** catalogues.                                                                            
      if(data_mock_flag == 0){
	    //** Nagoya v6. (v7.) mocks **//                                                                      
	    if(loopCount<10)       sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_00%d_spec_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, loopCount);
	    else if(loopCount<100) sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_0%d_spec_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
	    else                   sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_%d_spec_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);
    
	    CatalogueInput_500s(filepath);
      
  	    // Choice of redshift from zcos, zpec, zphot, zobs.                                               
	    gal_z = &zobs[0];

	    // load sampling according to local TSR.                                                                       
	    spec_weights();

        assignAcceptance();
      }

      // analysis on **data** catalogues.                                                                                   
      if(data_mock_flag == 1){
  	    // W1 catalogue.
	    sprintf(filepath, "%s/W1_Spectro_V7_0/W%d_SPECTRO_V7_0.txt", root_dir, fieldFlag);
    
  	    DataInput(filepath);
    
	    spec_weights();
	
	    // Acceptance criteria to be applied: 
	    //   i) Redshift cut
	    //  ii) zFlag cut, 2 to 9 inclusive. 
	    // iii) photoMask == 1 
	    assignAcceptance_WX_SPECTRO_V7();
      }
      
      for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
            chi               = interp_comovingDistance(zobs[j]);
	    
            Index             = (int) floor(chi/chi_interval);
            
            chibins[Index]   += chi/sampling[j];
            
            Nchi[Index]      +=  1./sampling[j]; 
            
            if((reflect==1) &&  (chi-loChi< 2.*kernel_width))   reflect_down(chi, j, Index);
            
            if((reflect==1) &&  (hiChi-chi< 2.*kernel_width))     reflect_up(chi, j, Index);
	    }  
	  }
    }
    
    for(j=0; j<chibin_no; j++){  
        chibins[j]  /= Nchi[j];
    
        if(Nchi[j] == 0) chibins[j] = (j + 0.5)*chi_interval;
    }
    
    for(j=0; j<chibin_no; j++)  comovVol[j] = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]  = Nchi[j]/comovVol[j]; 

    Gaussianfilter(nbar, chibin_no, chi_interval, kernel_width);
    
    if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/smoothed_nbar/nbar_smoothed_reflected_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg.dat", root_dir, kernel_width, loopCount);
    if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_0/nbar_smoothed_reflected_%.1lf_%.2lf_%.2lf_Nagoya_v7_WX_SPECTRO_V7_twofield_avg.dat", root_dir, kernel_width, lo_MBlim, hi_MBlim);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)  fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);

    // Restore global field Flag.
    fieldFlag = global_fieldFlag;
     
    return 0;
}


double nbar_dV(double chi){
    return pow(chi, 2.)*interp_nz(chi);
}


int prep_inverseCumulative_nbar(){
    // Calculate cumulative distribution of nbar and splint it's inverse.
    double                           norm;
    
    norm = qromb(&nbar_dV, 700., 2500.);
    
    for(i=0; i<100; i++){
        chi_cumulative_nbar[i] = 1250. + i*(2500. - 700.)/100.;
        
        cumulative_nbar[i]     = qromb(&nbar_dV, 700., chi_cumulative_nbar[i]);
        
        cumulative_nbar[i]    /= norm;
    }
    
    // for(i=0; i<58; i++)  printf("\n%e \t %e", chi_cumulative_nbar[i], cumulative_nbar[i]);
    
    spline(cumulative_nbar, chi_cumulative_nbar, 58, 1.0e31, 1.0e31, cumulative_nbar_2d);
    
    // printf("\n\n%e", inverse_cumulative_nbar(0.4));
    
    return 0;
}


double inverse_cumulative_nbar(double arg){
    // if F is the cumulative distribution of nbar.
    // return the inverse of F, for random generation.
    
    double result;

    splint(cumulative_nbar, chi_cumulative_nbar, cumulative_nbar_2d,  58, arg, &result);

    return result;
}
