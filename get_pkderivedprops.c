int set_clippingvars(void){
  u0             = inverse_erf(2.*sqrt(1./mean_suppression) - 1.);  // Estimate u0 from amplitude suppression.

  clipmono_amp   = 0.25*pow(1. + gsl_sf_erf(u0), 2.);
  
  clip_distcoeff = C_n(u0, 1);

  printf("\n\nHere, u0 is %.3lf", u0);
  
  return 0;
}

int get_mocksshotnoise(){  
  mean_shot = 0.0;

  sprintf(filepath, "%s/mocks_v1.7/pk_derivedprops/d0_%d/W%d/shotnoise_zlim_%.1lf_%.1lf.dat", covariance_mocks_path, d0, fieldFlag, lo_zlim, hi_zlim);
  
  inputfile = fopen(filepath, "r");

  linecount_header(inputfile, 1, &shot_ninstance);

  shotnoise_instances = malloc(shot_ninstance*sizeof(*shotnoise_instances));

  fscanf(inputfile, "%*[^\n]\n", NULL); // skip one line (hashed comment)
   
  for(j=0; j<shot_ninstance; j++){
    fscanf(inputfile, "%*d    %*lf    %lf\n", &shotnoise_instances[j]); // Without clipping, <n> estimate. 
    // fscanf(inputfile, "%*d    %lf    %*lf\n", &shotnoise_instances[j]); // With clipping, high-k estimate.
    
    // printf("\n%d \t %lf", j, shotnoise_instances[j]);
    
    mean_shot += shotnoise_instances[j];
  }

  mean_shot /= shot_ninstance;
  
  fclose(inputfile);

  if(shot_ninstance != CatalogNumber){
    printf("\n\nProblem with shot noise instances: number of mocks = %d, shot noise instances = %d", CatalogNumber, shot_ninstance);

    exit(EXIT_FAILURE);
  }

  else{
    printf("\n\nMean shot noise: %.4lf", mean_shot);
  }
  
  return 0;
}

int get_mocksclippedamplitudes(){
  mean_suppression = 0.0;

  sprintf(filepath, "%s/mocks_v1.7/pk_derivedprops/d0_%d/W%d/suppression_zlim_%.1lf_%.1lf.dat", covariance_mocks_path, d0, fieldFlag, lo_zlim, hi_zlim);
  
  inputfile = fopen(filepath, "r");
  
  line_count(inputfile, &suppression_ninstance);
  
  suppression_instances = malloc(suppression_ninstance*sizeof(*suppression_instances));
  
  for(j=0; j<suppression_ninstance; j++){
    fscanf(inputfile, "%*d \t %lf \n", &suppression_instances[j]);

    mean_suppression += suppression_instances[j];
  }
  
  fclose(inputfile);
  
  mean_suppression /= suppression_ninstance;
  
  if(suppression_ninstance != CatalogNumber){
    printf("\n\nProblem with shot noise instances: number of mocks = %d, suppression instances = %d \n\n", CatalogNumber, suppression_ninstance);

    exit(EXIT_FAILURE);
  }

  else{
    printf("\n\nMean amplitude correction: %.4lf", mean_suppression);
  }
  
  return 0;
}

int get_datashotnoise(){  
  sprintf(filepath, "%s/data_v1.7/pk_derivedprops/d0_%d/W%d/shotnoise_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim);
  
  inputfile = fopen(filepath, "r");

  if(inputfile == NULL){
    printf("\n\nError retrieving data shot noise.");

    exit(EXIT_FAILURE);
  }

  fscanf(inputfile, "%*[^\n]\n", NULL); // skip one line (hashed comment)
  
  fscanf(inputfile, "%*lf \t %lf \n", &mean_shot);  // Without clipping, <n> estimate.
  // fscanf(inputfile, "%lf \t %*lf \n", &mean_shot);  // With clipping, high-k estimate.
  
  fclose(inputfile);

  printf("\n\nData mean shot noise (%s): %.4lf", filepath, mean_shot);
  
  return 0;
}

int get_dataclippedamplitude(){  
  sprintf(filepath, "%s/data_v1.7/pk_derivedprops/d0_%d/W%d/suppression_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim);

  inputfile = fopen(filepath, "r");

  fscanf(inputfile, "%lf \n", mean_suppression);

  fclose(inputfile);

  printf("\n\nMean amplitude correction: %.4lf", mean_suppression);
  
  return 0;
}

int set_oldclippingvars(){
  if(data_mock_flag == 0){ // mocks
    if(d0      == 1000.)  mean_suppression    =  1.00;
    else if(d0 ==   10.)  mean_suppression    =  1.30;
    else if(d0 ==    6.)  mean_suppression    =  1.85;
    else                  mean_suppression    =  3.00;
  }

  if(data_mock_flag == 1){ // data
    if(fieldFlag == 1){
      if(d0      == 1000.)  mean_suppression    =  1.00;
      else if(d0 ==   10.)  mean_suppression    =  1.30;
      else if(d0 ==    6.)  mean_suppression    =  1.70;
      else                  mean_suppression    =  2.70;
    }

    if(fieldFlag == 4){
      if(d0      == 1000.)  mean_suppression    =  1.00;
      else if(d0 ==   10.)  mean_suppression    =  1.30;
      else if(d0 ==    6.)  mean_suppression    =  1.90;
      else                  mean_suppression    =  2.70;
    }
  }
      
  set_clippingvars();

  printf("\n\nHere, u0 is %.3lf", u0);
  
  return 0;
}

int set_oldshotnoise(){ // Clipping. 
  if(data_mock_flag == 0){ // mocks
    if(fieldFlag == 1){
      if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 277.16;   // Correct monopole.
      else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 243.94; 
      else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 208.99;
      else                  for(j=0; j<mono_order; j++) xdata[j] -= 167.90;
    }
  
    if(fieldFlag == 4){
      if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 289.16;   // Correct monopole.
      else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 249.26;  
      else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 208.00;
      else                  for(j=0; j<mono_order; j++) xdata[j] -= 168.59;
    }
  }

  else if(data_mock_flag == 1){ // data
    if(fieldFlag == 1){
      if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 335.538;  // Correct monopole.
      else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 283.993;  
      else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 229.602;
      else                  for(j=0; j<mono_order; j++) xdata[j] -= 180.127;
    }

    if(fieldFlag == 4){
      if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 331.20;   // Correct monopole. 
      else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 278.29;   
      else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 227.51;
      else                  for(j=0; j<mono_order; j++) xdata[j] -= 177.75;
    }
  }

  else{
    printf("\n\nAnother option required.");

    exit(EXIT_FAILURE);
  }
  
  return 0;
}
