int set_mem(void){
  grander = (double *)  calloc(mono_order, sizeof(*grander));
   tenner = (double *)  calloc(mono_order, sizeof(*tenner));
    sixer = (double *)  calloc(mono_order, sizeof(*sixer));
   fourer = (double *)  calloc(mono_order, sizeof(*fourer));
  
  return 0;
}


int get_a11kindices(int start, char filepath[]){
  char   firstfilepath[200];
  char  foldedfilepath[200];
  
  sprintf(firstfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, start, lo_zlim, hi_zlim);
    
  printf("\n\n%s", firstfilepath);
  
  inputfile  = fopen(firstfilepath, "r");
    
  line_count(inputfile, &lineNo);

  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim<ChiSq_kmin){
      chiSq_kminIndex = i + 1;
    }

    if(Interim > jenkins_fold_kjoin){
      jenkins_foldIndex_unfoldedfile = i;

      break;
    }
  }

  fclose(inputfile);

  printf("\n\nk for switching to folded measurement: %.3lf (%d)", jenkins_fold_kjoin, jenkins_foldIndex_unfoldedfile);

  // again, but for folded measurement.
  sprintf(foldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);

  inputfile  = fopen(foldedfilepath, "r");

  line_count(inputfile, &lineNo);

  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim < jenkins_fold_kjoin){
      jenkins_foldIndex_foldedfile = i + 1;
    }

    if(Interim > 0.8){
      chiSq_kmaxIndex = i;

      break;
    }
  }

  fclose(inputfile);

  mono_order = (jenkins_foldIndex_unfoldedfile - chiSq_kminIndex) + (chiSq_kmaxIndex - jenkins_foldIndex_foldedfile);
  order      = mono_order*hiMultipoleOrder;
  
  return 0;
}


int set_MeanMultipoles(double* array, int ld0, int mocks, int start){
  char     Nthfilepath[200];
  char Nfoldedfilepath[200];

  sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock", covariance_mocks_path, ld0, fieldFlag);
  
  for(k=0; k<mocks; k++){
    sprintf(Nthfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, k + start, lo_zlim, hi_zlim);

    inputfile = fopen(Nthfilepath, "r");

    for(i=0; i<jenkins_foldIndex_unfoldedfile; i++){
      if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

      else{
        fscanf(inputfile, "%*le \t %le \t %*le \t %*d \n", &Interim);

        array[i - chiSq_kminIndex] += Interim;
      }
    }

    fclose(inputfile);
    
    sprintf(Nfoldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, k + start, lo_zlim, hi_zlim);  // add in folded measurements, e.g. at k_join = 0.2;

    inputfile = fopen(Nfoldedfilepath, "r");

    for(i=0; i<chiSq_kmaxIndex; i++){
      if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

      else{
        fscanf(inputfile, "%*le \t %le \t %*le \t %*d \n", &Interim);

        array[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex] += Interim;
      }
    }

    fclose(inputfile);
  }

  for(k=0; k<mono_order; k++) array[k] /= mocks;
  
  return 0;
}


int amplitude_correction(double* array, double* sq_amp){
  // remember p(k) scales with squared factor.
  double mean_grander = 0.0;
  double meaner       = 0.0;
  double sq_factor    = 0.0;

  for(j=0; j<mono_order; j++){
    mean_grander += grander[j]; // only need ratio, number of modes averaged drops out. 
          meaner +=   array[j];

    if(all_kVals[j] > 0.09)  break;
  }
  
  sq_factor = mean_grander/meaner;

  *sq_amp = (1./sq_factor); // inverse factor. 

  printf("\n\nAmplitude Correction:  P_xx(k) = (%.6lf)^2  P_1000(k)", sqrt(*sq_amp));

  return 0;
}

  
int get_shotnoise(int ld0, int start, int mocks, double klo, double khi, double* shot_amp){
  // Needs multiply folded measurements to get to sufficiently high k-max without aliasing. 

  double k, P0;
  int    nmodes = 0;
  char   Nthfilepath[200];
  
  
  *shot_amp = 0.0;

  for(i=start; i<mocks; i++){
    sprintf(   filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock", covariance_mocks_path, ld0, fieldFlag);
    sprintf(Nthfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_4.dat", filepath, i + start, lo_zlim, hi_zlim);

    // printf("\n\nHERE: %s", Nthfilepath);
    
    inputfile = fopen(Nthfilepath, "r");

    line_count(inputfile, &lineNo);
   
    for(j=0; j<lineNo; j++){
      fscanf(inputfile, "%le \t %le \t %*le \t %*d \n", &k, &P0);
      
      if((klo <= k) && (k <= khi)){
        *shot_amp += P0;

        nmodes    += 1;
      }
    }

    fclose(inputfile);
  }
  
  *shot_amp /= nmodes; // mean power in the interval k_lo < k < k_hi.

  printf("\n\nShot noise correction:  P_shot = %.6lf", *shot_amp);
  
  return 0;
}


double prep_a11(int mocks, int start, double* sq_amp, double* shot_amp){
  double klo = 1.5;
  double khi = 2.0;

  sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock", covariance_mocks_path, 1000, fieldFlag); // can use any d0 to establish kVals. 
  
  get_a11kindices(start, filepath);

  set_mem(); // uses mono_order, as set by get_kindices. 
  
  set_MeanMultipoles(grander, 1000, mocks, start);
  set_MeanMultipoles( tenner,   10, mocks, start);
  set_MeanMultipoles(  sixer,    6, mocks, start);
  set_MeanMultipoles( fourer,    4, mocks, start);
  
  if(d0 == 1000){
    amplitude_correction(grander, sq_amp);
    
    get_shotnoise(1000, start, mocks, klo, khi, shot_amp);
  }
  
  else if(d0 == 10){
    set_MeanMultipoles(tenner, 10, mocks, start);

    amplitude_correction(tenner, sq_amp);

    get_shotnoise(10, start, mocks, klo, khi, shot_amp);
  }

  else if(d0 == 6){
    set_MeanMultipoles(sixer, 6, mocks, start);
    
    amplitude_correction(sixer, sq_amp);

    get_shotnoise(6, start, mocks, klo, khi, shot_amp);
  }

  else if(d0 == 4){
    set_MeanMultipoles(fourer, 4, mocks, start);
    
    amplitude_correction(fourer, sq_amp);

    get_shotnoise(4, start, mocks, klo, khi, shot_amp);
  }

  else{
    printf("\n\nMocks for given threshold are not available.");

    exit(EXIT_FAILURE);
  }


  u0 = inverse_erf(2.*sqrt(*sq_amp) - 1.);  // Estimate u0 from amplitude suppression.

  clip_distcoeff = C_n(u0, 1);

  return 0;
}
