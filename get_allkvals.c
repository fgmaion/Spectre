int get_allkvals(int start){
  char   firstfilepath[200];
  char  foldedfilepath[200];

  sprintf(filepath,      "%s/mocks_v1.7/pk/d0_%d/W%d/mock",   covariance_mocks_path, d0, fieldFlag);
  
  sprintf(firstfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_0.dat", filepath, start, lo_zlim, hi_zlim);
  
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

  allmono_order = jenkins_foldIndex_unfoldedfile + lineNo - jenkins_foldIndex_foldedfile;
    
  all_kVals     = (double  *) malloc(allmono_order*sizeof(*all_kVals));
  
  inputfile  = fopen(firstfilepath, "r");

  for(i=0; i<jenkins_foldIndex_unfoldedfile; i++)  fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &all_kVals[i]);

  fclose(inputfile);
  
  sprintf(foldedfilepath, "%s_%03d_zlim_%.1lf_%.1lf_Jf_2.dat", filepath, start, lo_zlim, hi_zlim);  // add in folded measurements, e.g. at k_join = 0.2;
  
  inputfile = fopen(foldedfilepath, "r");
  
  for(i=0; i<lineNo; i++){
    if(i<jenkins_foldIndex_foldedfile) fscanf(inputfile, "%*le \t %*le \t %*le \t %*d \n");

    else{
      fscanf(inputfile, "%le \t %*le \t %*le \t %*d \n", &all_kVals[jenkins_foldIndex_unfoldedfile + i - jenkins_foldIndex_foldedfile]);
    }
  }
  
  fclose(inputfile);
  
  return 0;
}


int allkvals_matchup(){
  double     diff;
  double min_diff;

  fftlog_indices = realloc(fftlog_indices, allmono_order*sizeof(*fftlog_indices));

  for(i=0; i<allmono_order; i++){
    min_diff = pow(10., 99.);

    for(j=0; j<FFTlogRes; j++){
      diff = fabs(mono_config->krvals[j][0] - all_kVals[i]);

      if(diff<min_diff){
        min_diff = diff;

        fftlog_indices[i]  = j;
      }
    }
  }

  printf("\n\nAll k vals matchup:");
  
  for(j=0; j<allmono_order; j++)  printf("\n%.6le \t %.6le", all_kVals[j], mono_config->krvals[fftlog_indices[j]][0]);
  
  return 0;
}
