int load_withoutfolding(char filepath[]){
  char  firstfilepath[200];

  sprintf( firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, lo_zlim, hi_zlim);

  inputfile = fopen(firstfilepath, "r");

  for(i=0; i<chiSq_kmaxIndex; i++){
    if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%*le \t %le \t %le \t %*d \n", &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);
    }
  }

  fclose(inputfile);

  return 0;
}


int load_withfolding(char filepath[]){
  double loc_k;

  char  firstfilepath[200];
  char foldedfilepath[200];

  sprintf( firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, lo_zlim, hi_zlim);
  sprintf(foldedfilepath, "%s_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, lo_zlim, hi_zlim);

  inputfile = fopen(firstfilepath, "r");

  for(i=0; i<jenkins_foldIndex_unfoldedfile; i++){
    if(i<chiSq_kminIndex) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &loc_k, &xdata[i - chiSq_kminIndex], &xdata[mono_order + i - chiSq_kminIndex]);
    }
  }

  fclose(inputfile);

  inputfile = fopen(foldedfilepath, "r");  // Add in folded measurements from e.g. jenkins_fold_kjoin = 0.1;

  for(i=0; i<chiSq_kmaxIndex; i++){
    if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%le \t %le \t %le \t %*d \n", &loc_k, &xdata[i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex], &xdata[mono_order + i - jenkins_foldIndex_foldedfile + jenkins_foldIndex_unfoldedfile - chiSq_kminIndex]);
    }
  }

  fclose(inputfile);

  return 0;
}


int load_mock(int mockNumber){
  //  Chapter 7, snipping results.  Snipped results are raw power, need shot noise corrected.
  // d0 = {1000, 10, 6}, Ps = {277.16, 243.944925, 208.99925}.  Determined from Jf_8 results.

  sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock_%03d", covariance_mocks_path, d0, fieldFlag, mockNumber);

  printf("\n\n%s", filepath);
  
  if(ChiSq_kmax <= jenkins_fold_kjoin){
    // printf("\n\nLoading mock without folding");

    load_withoutfolding(filepath);  // amplitude rescaling in load_withoutfolding.
  }

  else{
    // printf("\n\nLoading mock with folding");

    load_withfolding(filepath);  // amplitude recsaling in load_withfolding.
  }

  // snipping_amplitudeCorrection(xdata, mono_order);  //

  // snipping_shotnoise(xdata, mono_order);  // corrects monopole only;

  // printf("\n\nCorrelated data.");

  // for(j=0; j<mono_order; j++)  printf("\n%le \t %le \t %le", kVals[j], xdata[j], xdata[j + mono_order]);

  for(j=0; j<mono_order; j++){
    kdata[j]              = kVals[j];
    kdata[j + mono_order] = kVals[j];
  }
  
  return 0;
}


int load_data(){
  // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/d0_%d/W%d/data", root_dir, d0, fieldFlag);
  // sprintf(filepath, "%s/W1_Spectro_V7_2/data_v1.7/pk/old_d0_%d/W%d/data", root_dir, d0, fieldFlag);
  sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/pk/d0_%d/W%d/data", root_dir, d0, fieldFlag);

  if(ChiSq_kmax <= jenkins_fold_kjoin){
    printf("\n\nLoading data without folding");

    load_withoutfolding(filepath);
  }

  else{
    printf("\n\nLoading data with folding");

    load_withfolding(filepath);  // Ignore reference to mock, should work just fine.
  }

  if(fieldFlag == 1){
    if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 335.538;  // Correct monopole.
    else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 283.993;  // Note: as a difference estimator, covariance does not have to explicity corrected.
    else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 229.602;
    else if(d0 ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 180.127;
  }

  if(fieldFlag == 4){
    if(d0      == 1000.)  for(j=0; j<mono_order; j++) xdata[j] -= 331.20;  // Correct monopole.
    else if(d0 ==   10.)  for(j=0; j<mono_order; j++) xdata[j] -= 278.29;  // Note: as a difference estimator, covariance does not have to explicity corrected.
    else if(d0 ==    6.)  for(j=0; j<mono_order; j++) xdata[j] -= 227.51;
    else if(d0 ==    4.)  for(j=0; j<mono_order; j++) xdata[j] -= 177.75;
  }
  
  return 0;
}
  

int set_meanmultipoles(){
  for(j=0; j<order; j++) xdata[j] = MeanMultipoles[j];  
  
  return 0;
}
