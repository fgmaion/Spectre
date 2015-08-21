int prep_nbar(){
  chibin_no =  (int)       ceil((interp_comovingDistance(1.4) - interp_comovingDistance(0.0))/chi_interval);   
                                                                                                                             
  zbins     =  (double *)  realloc(zbins,    chibin_no*sizeof(*zbins));
  chibins   =  (double *)  realloc(chibins,  chibin_no*sizeof(*chibins));
  
  Nchi      =  (double *)  realloc(Nchi,     chibin_no*sizeof(*Nchi));
  nbar      =  (double *)  realloc(nbar,     chibin_no*sizeof(*nbar));

  comovVol  =  (double *)  realloc(comovVol, chibin_no*sizeof(*comovVol));
  
  // Second derivatives.   
  nbar_2d   =  (double *)  realloc(nbar_2d,  chibin_no*sizeof(*nbar_2d));

  for(j=0; j<chibin_no;   j++){ 
      zbins[j]    = 0.0;
      chibins[j]  = 0.0;
      Nchi[j]     = 0.0;
      nbar[j]     = 0.0;  
      comovVol[j] = 0.0;
  }

  return 0;
}


int nbar_calc(int mocks){    
    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    prep_nbar();
    
    double chi;
    
    int global_fieldFlag;

    // Store initial fieldFlag. 
    global_fieldFlag = fieldFlag;
    
    for(loopCount=1; loopCount<mocks+1; loopCount++){
        printf("\n\n%d", loopCount);
        
        for(ii=1; ii<5; ii=ii+3){
          fieldFlag = ii;
      
          if(loopCount<10)        sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_00%d_spec_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, loopCount);
          else if(loopCount<100)  sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_0%d_spec_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
          else                    sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_%d_spec_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);
      
          // Choice of redshift from zcos, zpec, zphot, zobs.
          CatalogueInput_500s(filepath);
          // CatalogueInput_mockTSR(filepath);
        
          gal_z = &zobs[0];
        
          // currently loads local TSR (Granett) weights for mock TSR catalogue.
          spec_weights();
        
          assignAcceptance();
        
          for(j=0; j<Vipers_Num; j++){
            chi             = interp_comovingDistance(zobs[j]);
            
            Index           = (int) floor(chi/chi_interval);
            
            chibins[Index] += chi/sampling[j];
            

            Nchi[Index]    +=  1./sampling[j];
          }
        } 
    }
    
    for(j=0; j<chibin_no; j++){  
        chibins[j]  /= Nchi[j];
    
        if(Nchi[j] == 0) chibins[j] = (j + 0.5)*chi_interval;
    }
    
    // avg. number of gals in bin dChi in one mock. 
    for(j=0; j<chibin_no; j++)  Nchi[j]     /= mocks;
    
    for(j=0; j<chibin_no; j++)  comovVol[j]  = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]   = Nchi[j]/comovVol[j]; 
    
    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mockavg_specweight_nbar_%.2lf_%.2lf_Nagoya_v6_Samhain.dat", root_dir, lo_MBlim, hi_MBlim); 
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)   fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);
    
    // Restore global field Flag.
    fieldFlag = global_fieldFlag;
    
    return 0;
}


int randoms_nbar_calc(){    
    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    
    // chi_interval = 1.;
    
    prep_nbar();
    
    double chi;
        
    for(j=0; j<rand_number; j++){
        Index           = (int) floor(rand_chi[j]/chi_interval);
            
        chibins[Index] += rand_chi[j];
            
        Nchi[Index]    +=  1.;

    } 
    
    for(j=0; j<chibin_no; j++){  
        chibins[j]  /= Nchi[j];
    
        if(Nchi[j] == 0) chibins[j] = (j + 0.5)*chi_interval;
    }

    for(j=0; j<chibin_no; j++)  comovVol[j]  = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]   = Nchi[j]/comovVol[j]; 



    sprintf(filepath, "%s/W1_Spectro_V7_0/randoms_chi_reassignment_nbar_-90.50_-0.00_Nagoya_v7_Samhain.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)   fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);
    
    return 0;
}


int spline_nbar(int data_mock_flag, int twoField, int truth){
    prep_nbar();

    // analysis on **mock** catalogues. 
    if(data_mock_flag == 0){
        //** Nagoya v4 mask **//
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_parent.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specmask.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_spec.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specweight.dat", root_dir, lo_MBlim, hi_MBlim);
    
        //** Nagoya v4, mock TSR, downweight parent on a per quadrant basis. **//
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockTSR.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockTSR_weighted.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockTSR_localweighted.dat", root_dir, lo_MBlim, hi_MBlim);

        //** Nagoya v4, mock spec cats, overlay grid across survey.  (i) Poisson sample on a per cell basis, mean 40 %. (ii) Non-Poisson, one per cell.**//
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockSpec.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockSpec_nonpoisson.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockSpec_nonpoisson_weighted.dat", root_dir, lo_MBlim, hi_MBlim);

        //** Nagoya v4, non-Poisson sampling. One galaxy per cell.  
        // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_mockSpec_nonpoisson_Granettweighted.dat", root_dir, lo_MBlim, hi_MBlim);
    
        if(twoField == 0){
	        //** Nagoya v6., local TSR weighted. **//
	        sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/nbar/specweight_nbar_%.2lf_%.2lf_Nagoya_v6_Samhain.dat", root_dir, lo_MBlim, hi_MBlim); 
	    }
	    
	    if(twoField == 1){
	        //** Nagoya v7., Gaussian smoothed, two field estimate, as for the data. local TSR weighted. **//
	        // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/smoothed_nbar/nbar_smoothed_Nagoya_v7_Samhain_mock_%d_twofield_avg.dat", root_dir, loopCount);
	    
    	    if(truth==0){
	            // Includes reflection.
	            sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/smoothed_nbar/nbar_smoothed_reflected_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg.dat", root_dir, nz_smoothRadius, loopCount);
	        }
	        
	        if(truth==1){
	            sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/mockavg_specweight_nbar_-90.50_-0.00_Nagoya_v6_Samhain.dat", root_dir);
	        }
	    }
    }
    
    // analysis on **VIPERS DATA** catalogues. 
    if(data_mock_flag == 1){
        //** Nagoya v5. and Samhain **//
        // sprintf(filepath, "%s/W1_Spectro_V5_0/smoothed_nbar_-90.50_-0.00_Nagoya_v5_W1W4_realdata.dat", root_dir, lo_MBlim, hi_MBlim);
        // sprintf(filepath, "%s/W1_Spectro_V5_0/mocks_nbar_%.2lf_%.2lf_Nagoya_v5_Samhain_spec.dat", root_dir, lo_MBlim, hi_MBlim);   
    
        //** Nagoya v7 and Samhain **//
        // sprintf(filepath, "%s/W1_Spectro_V7_0/nbar_smoothed_%.2lf_%.2lf_Nagoya_v7_WX_SPECTRO_V7_twofield_avg.dat", root_dir, lo_MBlim, hi_MBlim);    
    
        //** Nagoya v7 and Samhain with reflection **//
        sprintf(filepath, "%s/W1_Spectro_V7_0/nbar_smoothed_reflected_%.1lf_%.2lf_%.2lf_Nagoya_v7_WX_SPECTRO_V7_twofield_avg.dat", root_dir, nz_smoothRadius, lo_MBlim, hi_MBlim);    
    }
    
    printf("\n\nn bar: %s\n\n", filepath);
    
    inputfile = fopen(filepath, "r");

    for(j=0; j<chibin_no; j++)  fscanf(inputfile, "%le \t %le \n", &chibins[j], &nbar[j]);
    
    fclose(inputfile);

    spline(chibins, nbar, chibin_no, 1.0e31, 1.0e31, nbar_2d);
        
    prep_inverseCumulative_nbar();

    return 0;
}


double interp_nz(double chi){
    splint(chibins, nbar, nbar_2d,  chibin_no, chi, &Interim);
    
    return Interim;
}
