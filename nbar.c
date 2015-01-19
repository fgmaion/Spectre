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
    
    double chi, chi_lo;
    
    chi_lo = interp_comovingDistance(0.0);
    
    for(loopCount=1; loopCount<mocks; loopCount++){
        printf("\n\n%d", loopCount);
    
        // new 500s. limit to 0.7<z<0.8, limit to linear bias of 1.495903 corresponding to ~ -20.0 mag galaxies. at this mag. volume limited to z=0.85
        if(loopCount<10)        sprintf(filepath, "%s/mocks_W1_v9.0_500/mock_00%d_parent.dat", vipersHOD_dir, loopCount);
        else if(loopCount<100)  sprintf(filepath, "%s/mocks_W1_v9.0_500/mock_0%d_parent.dat",  vipersHOD_dir, loopCount);
        else                    sprintf(filepath, "%s/mocks_W1_v9.0_500/mock_%d_parent.dat",  vipersHOD_dir, loopCount);
      
        // Choice of redshift from zcos, zpec, zphot, zobs.
        CatalogueInput_500s(filepath);
        
        for(k=0; k<Vipers_Num; k++){
            if((-20.2<M_B[k]) && (M_B[k] < -19.8)){        
                chi             = interp_comovingDistance(zobs[k]);
            
                Index           = (int) floor((chi - chi_lo)/chi_interval);
            
                Nchi[Index]    += 1.;
            }
        } 
    }
    
    for(j=0; j<chibin_no; j++)     Nchi[j] /= mocks;
    
    for(j=0; j<chibin_no; j++)  chibins[j]  = (j+0.5)*chi_interval;
    
    for(j=0; j<chibin_no; j++) comovVol[j]  = sqdegs2steradians(W1area)*(pow(chibins[j+1], 3.) - pow(chibins[j], 3.))/3.;
    
    for(j=0; j<chibin_no; j++)     nbar[j]  = Nchi[j]/comovVol[j]; 
    
    
    sprintf(filepath, "%s/Data/likelihood/chi_nbar_0.7_z_0.8_20.0_mags.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<chibin_no; j++)   fprintf(output, "%e \t %e \n", chibins[j], nbar[j]);
    
    fclose(output);
    
    return 0;
}


int spline_nbar(){
    prep_nbar();

    sprintf(filepath, "%s/Data/likelihood/chi_nbar_0.7_z_0.8_20.0_mags.dat", root_dir);

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
