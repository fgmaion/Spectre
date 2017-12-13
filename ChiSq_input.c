int load_withoutfolding(char filepath[]){
  char  firstfilepath[200];

  if(mull == 0){
    sprintf(firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, lo_zlim, hi_zlim);
  }

  else{
    sprintf(firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, lo_zlim, hi_zlim);
  }
  
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
  double             loc_k;

  char  firstfilepath[200];
  char foldedfilepath[200];

  if(mull == 0){
    sprintf(firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, lo_zlim, hi_zlim);
  }

  else{
    sprintf(firstfilepath, "%s_zlim_%.1lf_%.1lf_Jf_1.dat", filepath, lo_zlim, hi_zlim);
  }
  
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
  if(mull == 0){
    sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock_%03d", covariance_mocks_path, d0, fieldFlag, mockNumber);
  }

  else{
    sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock_%d", covariance_mocks_path, d0, fieldFlag, mockNumber);
  }
    
  if(ChiSq_kmax <= jenkins_fold_kjoin){
    load_withoutfolding(filepath);  // amplitude rescaling no longer in load_withoutfolding.
  }

  else{
    load_withfolding(filepath);  // amplitude recsaling no longer in load_withfolding.
  }

  // printf("\n\nCorrelated data.");
  // for(j=0; j<mono_order; j++)  printf("\n%le \t %le \t %le", kVals[j], xdata[j], xdata[j + mono_order]);

  for(j=0; j<mono_order; j++){
    kdata[j]              = kVals[j];
    kdata[j + mono_order] = kVals[j];
  }
  
  return 0;
}

int load_data(void){
  sprintf(filepath, "%s/data_v1.7/pk/d0_%d/W%d/data", covariance_mocks_path, d0, fieldFlag);  // covariance_mocks_path updated according to mull flag.
    
  if(ChiSq_kmax <= jenkins_fold_kjoin){
    load_withoutfolding(filepath);
  }

  else{
    load_withfolding(filepath);
  }

  for(j=0; j<mono_order; j++){
    kdata[j]              = kVals[j];
    kdata[j + mono_order] = kVals[j];
  }
  
  return 0;
}
