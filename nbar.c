int prep_nbar(){
  chibin_no =  (int)       ceil((interp_comovingDistance(2.0) - interp_comovingDistance(0.0))/chi_interval);   
                                                                                                                             
  zbins     =  (double *)  realloc(zbins,    chibin_no*sizeof(*zbins));
  chibins   =  (double *)  realloc(chibins,  chibin_no*sizeof(*chibins));
  Nchi      =  (double *)  realloc(Nchi,     chibin_no*sizeof(*Nchi));
  nbar      =  (double *)  realloc(nbar,     chibin_no*sizeof(*nbar));
  comovVol  =  (double *)  realloc(comovVol, chibin_no*sizeof(*comovVol));
  nbar_2d   =  (double *)  realloc(nbar_2d,  chibin_no*sizeof(*nbar_2d));

  for(j=0; j<chibin_no; j++){ 
      Nchi[j]     = 0.0;
      nbar[j]     = 0.0;  
      comovVol[j] = 0.0;
  
      chibins[j] = (j + 0.5)*chi_interval;  // initialise chi vals.
    
      zbins[j] = interp_inverseComovingDistance(chibins[j]);
  }

  return 0;
}


int spline_nbar(int truth){
    prep_nbar();

    if(data_mock_flag == 0){	 // analysis on mocks.   
      if(truth==0){
        // smoothed counts, 1% renormalisation to \sum E^-1 (sum over W1 and W4).
        sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_mock_%d_twofield_avg_v2.dat", root_dir, nz_smoothRadius, loopCount);
      }
	        
      if(truth==1){
        // ** PARENT **
        sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/nbar_100_smoothedCounts/nbar_smooth_0.0_Nagoya_v7_Samhain_parent_mocks_avg_twofield_avg.dat", root_dir);
      }
     }
    
    if(data_mock_flag == 1){  // analysis on data.
      // Post Stefano comparison. Smoothed counts, normalised to joint field galaxy count.
      sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/nbar_100_smoothedCounts/nbar_smooth_%.1lf_Nagoya_v7_Samhain_twofield_avg.dat", root_dir, nz_smoothRadius); 
    }

    printf("\n\n%s", filepath);
    
    inputfile = fopen(filepath, "r");

    for(j=0; j<chibin_no; j++)  fscanf(inputfile, "%le \t %le \n", &chibins[j], &nbar[j]);
        
    fclose(inputfile);
    
    spline(chibins, nbar, chibin_no, 1.0e31, 1.0e31, nbar_2d);

    prep_inverseCumulative_nbar();
        
    return 0;
}

int nbar_calc(int mocks){    
    double chi;
    
    int global_fieldFlag;

    global_fieldFlag = fieldFlag;  // Store initial fieldFlag.

    prep_nbar();  // Equal intervals in comoving distance, for both W1 and W4 fields.

    for(loopCount=1; loopCount<mocks+1; loopCount++){
        printf("\n\n%d", loopCount);
        
        for(ii=1; ii<5; ii=ii+3){
          fieldFlag = ii;
      
	  sprintf(filepath, "%s/mocks_v1.7/W%d/mock_%03d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
      
          CatalogueInput_500s(filepath);
        
          spec_weights();

	  assignAcceptance_true();
        
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
   
      Nchi[j]     /= mocks;
    
      comovVol[j]  = sqdegs2steradians(TotalW1W4area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
      nbar[j]   = Nchi[j]/comovVol[j]; 
    }

    sprintf(filepath, "%s/W1_Spectro_V7_2/mocks_v1.7/nbar_100_smoothedCounts/nbar_smooth_0.0_Nagoya_v7_Samhain_mock_avg_twofield_avg.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)   fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);

    fclose(output);
    
    // Restore global field Flag.
    fieldFlag = global_fieldFlag;
    
    return 0;
}


double interp_nz(double chi){
    splint(chibins, nbar, nbar_2d,  chibin_no, chi, &Interim);
    
    return Interim;
}
