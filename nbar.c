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
    
    for(loopCount=1; loopCount<mocks+1; loopCount++){
        printf("\n\n%d", loopCount);
    
        // new 500s. limit to 0.7<z<0.8, limit to linear bias of 1.495903 corresponding to ~ -20.0 mag galaxies. at this mag. volume limited to z=0.85
        // if(loopCount<10)        sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_lm_00%d_gal.dat", vipersHOD_dir, loopCount);
        // else if(loopCount<100)  sprintf(filepath, "%s/mocks_W1_v9.0_500/mock_lm_0%d_gal.dat",  vipersHOD_dir, loopCount);
        // else                    sprintf(filepath, "%s/mocks_W1_v9.0_500/mock_lm_%d_gal.dat",   vipersHOD_dir, loopCount);
      
        if(loopCount<10)        sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_00%d_spec.dat", vipersHOD_dir, loopCount);
        else if(loopCount<100)  sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_0%d_spec.dat",  vipersHOD_dir, loopCount);
        else                    sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_%d_spec.dat",   vipersHOD_dir, loopCount);
      
        // Choice of redshift from zcos, zpec, zphot, zobs.
        CatalogueInput_500s(filepath);
        
        spec_weights();
        
        for(j=0; j<Vipers_Num; j++){
            if((lo_MBlim<M_B[j]) && (M_B[j]<hi_MBlim)){        
                chi             = interp_comovingDistance(zobs[j]);
            
                Index           = (int) floor(chi/chi_interval);
            
                // chibins[Index] += chi;
                chibins[Index] += chi/sampling[j];
            
                // Nchi[Index]    +=  1.;
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
    
    for(j=0; j<chibin_no; j++)  comovVol[j]  = sqdegs2steradians(W1area)*(pow((j+1)*chi_interval, 3.) - pow(j*chi_interval, 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]   = Nchi[j]/comovVol[j]; 
    
    
    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_parent.dat", root_dir, lo_MBlim, hi_MBlim);
    sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_v4_Nagoya_v4_spec.dat", root_dir, lo_MBlim, hi_MBlim);
    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specmask.dat", root_dir, lo_MBlim, hi_MBlim);
    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specweight.dat", root_dir, lo_MBlim, hi_MBlim);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)   fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);
    
    return 0;
}


int spline_nbar(){
    prep_nbar();

    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_parent.dat", root_dir, lo_MBlim, hi_MBlim);
    sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_spec.dat", root_dir, lo_MBlim, hi_MBlim);
    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specmask.dat", root_dir, lo_MBlim, hi_MBlim);
    // sprintf(filepath, "%s/Data/500s/nbar_%.2lf_%.2lf_Nagoya_v4_specweight.dat", root_dir, lo_MBlim, hi_MBlim);

    inputfile = fopen(filepath, "r");

    for(j=0; j<chibin_no; j++)  fscanf(inputfile, "%le \t %le \n", &chibins[j], &nbar[j]);
    
    fclose(inputfile);

    spline(chibins, nbar, chibin_no, 1.0e31, 1.0e31, nbar_2d);
    
    // for(j=0; j<chibin_no; j++)  printf("%le \t %le \t %le \n", chibins[j], nbar[j], interp_nz(chibins[j]));

    return 0;
}


double interp_nz(double chi){
    splint(chibins, nbar, nbar_2d,  chibin_no, chi, &Interim);
    
    return Interim;
}
