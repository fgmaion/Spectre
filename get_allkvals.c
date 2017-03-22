int get_allkvals(int start, char filepath[]){
  char   firstfilepath[200];
  char  foldedfilepath[200];

  sprintf(firstfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, start, lo_zlim, hi_zlim);

  // printf("\n\n%s", firstfilepath);
  
  inputfile  = fopen(firstfilepath, "r");
  
  line_count(inputfile, &lineNo);
  
  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim > jenkins_fold_kjoin){
      jenkins_foldIndex_unfoldedfile = i;

      break;
    }
  }
  
  fclose(inputfile);
  
  // file with folded measurements.
  sprintf(foldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);

  inputfile  = fopen(foldedfilepath, "r");

  line_count(inputfile, &lineNo);

  for(i=0; i<lineNo; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &Interim);

    if(Interim < jenkins_fold_kjoin){
      jenkins_foldIndex_foldedfile = i + 1;
    }
  }

  fclose(inputfile);
  
  mono_order = jenkins_foldIndex_unfoldedfile + lineNo - jenkins_foldIndex_foldedfile;

  kVals = (double  *) realloc(kVals, mono_order*sizeof(*kVals));
  
  inputfile  = fopen(firstfilepath, "r");

  for(i=0; i<jenkins_foldIndex_unfoldedfile; i++){
    fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &kVals[i]);

    // printf("\n%.6lf", kVals[i]);
  }

  fclose(inputfile);
  
  sprintf(foldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);  // add in folded measurements, e.g. at k_join = 0.2;
  
  inputfile = fopen(foldedfilepath, "r");
  
  for(i=0; i<lineNo; i++){
    if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &kVals[jenkins_foldIndex_unfoldedfile + i - jenkins_foldIndex_foldedfile]);

      // printf("\n%.6lf", kVals[jenkins_foldIndex_unfoldedfile + i - jenkins_foldIndex_foldedfile]);
    }
  }
  
  fclose(inputfile);

  // for(i=0; i<mono_order; i++)  printf("\n%.9lf", kVals[i]);
  
  return 0;
}
