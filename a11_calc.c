int set_mem(void){
  grander = (double *)  calloc(mono_order, sizeof(*grander));
   tenner = (double *)  calloc(mono_order, sizeof(*tenner));
    sixer = (double *)  calloc(mono_order, sizeof(*sixer));
   fourer = (double *)  calloc(mono_order, sizeof(*fourer));
  
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


double a11_scale(int mocks, int start){
  set_mem();
  
  set_MeanMultipoles(grander, 1000, mocks, start);
  set_MeanMultipoles( tenner,   10, mocks, start);
  set_MeanMultipoles(  sixer,    6, mocks, start);
  set_MeanMultipoles( fourer,    4, mocks, start);

  for(k=0; k<mono_order; k++)  printf("\nHERE: %.6lf", tenner[k]);
  
  return 0;
}
