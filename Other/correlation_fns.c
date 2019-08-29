double ra_criterion(double x, double y, double max){
	if(fabs(x-y) > max){
		return 180. - fmod(fabs(x-y), 180.);
	}

	else{
		return fmod(fabs(x-y), 180.);
	}
}


double dec_criterion(double x, double y){
	double dec360x;
	double dec360y;

	if(x<0.)  dec360x = 360. + x;
	else	  dec360x = x;

	if(y<0.)  dec360y = 360. + y;
	else	  dec360y = y;

	return ra_criterion(dec360x, dec360y, 180.);
}

int initialise_angularCorrelation(){
  rabinInterval  = (1./240.)/50.;

  decbinInterval = (1./20.)/50.; 

  ra_bins  = malloc(radecbinNumb*sizeof(double));
  dec_bins = malloc(radecbinNumb*sizeof(double));

  for(j=0; j<radecbinNumb;  j++)        ra_bins[j]  =   rabinInterval*(j+0.5);
  for(j=0; j<radecbinNumb;  j++)       dec_bins[j]  =  decbinInterval*(j+0.5);

  binned_pairs = realloc(binned_pairs, radecbinNumb*sizeof(double*));

  for(j=0; j<radecbinNumb; j++)  binned_pairs[j] = malloc(radecbinNumb*sizeof(double));

  return 0;
}

int correlationfn(char filepath[], double xCoor[], double yCoor[], double zCoor[], int Nobjects){
  double r, mu;

  double* rbins;

  int rbin;

  rbins = malloc(300*sizeof(double));

  double rbinInterval = 1.;

  output = fopen(filepath, "w");

  for(i=0; i<Nobjects; i++){
    printf("\n%d \t %d", loopCount, i);

    for(j=i+1; j<Nobjects; j++){
      //      if((flag_SSPOC[i] == 1) && (flag_SSPOC[j] ==1)){

      r = pow(xCoor[i] - xCoor[j], 2.) + pow(yCoor[i] - yCoor[j], 2.) + pow(zCoor[i] - zCoor[j], 2.);

      r = pow(r, 0.5);
      
      // mu = fabs(zCoor[i] - zCoor[j])/r;

      if(r<300.){
	rbin  = (int) floor(r/rbinInterval);
      
	rbins[rbin] += 1;
      }
      // }
    }
 

  if(i%10000 == 0){
    output = fopen(filepath, "w");

    for(j=0; j<300; j++)  fprintf(output, "%e \t %e \n", 0.5*j, rbins[j]);

    fclose(output);
  }
  }

  output = fopen(filepath, "w");

  for(j=0; j<300; j++)  fprintf(output, "%e \t %e \n", 0.5*j, rbins[j]);

  fclose(output);

    return 0;
}

int randoms_angular_correlationfn(){
  // printf("\n\nHighest dec: %e", (dec_bins[radecbinNumb-1]+0.5*decbinInterval)*60.*60.);

  for(i=0; i<rand_number; i++){
    printf("\n%d \t %d", loopCount, i);

    for(j=i+1; j<rand_number; j++){
      //                              if((flag_SSPOC[i] == 1) && (flag_SSPOC[j] ==1)){               \
                                                                                                                  
      dra  = fabs(rand_ra[i] - rand_ra[j]);   // ra_criterion(ra[i],   ra[j], 180.)/180.;                     

      ddec = fabs(rand_dec[i] - rand_dec[j]); // dec_criterion(dec[i], dec[j])/180.;                          

      // printf("\n%e \t %e", dra, ddec);

      if((dra<ra_bins[radecbinNumb-1]+0.5*rabinInterval) && (ddec<dec_bins[radecbinNumb-1]+0.5*decbinInterval)){
	rabin  = (int) floor(dra/rabinInterval);

	decbin = (int) floor(ddec/decbinInterval);

	// printf("\n%d \t %d", rabin, decbin);                                           

	// if((dra < 3./(60.*60.)) && (ddec < 72./(60.*60.))){
	  binned_pairs[rabin][decbin] += 1.0;
	  //}
      }
      //      }                                                                                 
    }

    if(i%10000 == 0){

      printf("\nprinting to file");
      printMockAvg_correlation();
    }
  }

  printMockAvg_correlation();

  return 0;
}



int angular_correlationfn(){
	for(i=0; i<Vipers_Num; i++){
	  printf("\n%d \t %d", loopCount, i);
	  
		for(j=i+1; j<Vipers_Num; j++){
		  //		  		  if((flag_SSPOC[i] == 1) && (flag_SSPOC[j] ==1)){		  
		    dra  = fabs(ra[i] - ra[j]);   // ra_criterion(ra[i],   ra[j], 180.)/180.;
			 
		    ddec = fabs(dec[i] - dec[j]); // dec_criterion(dec[i], dec[j])/180.; 

			if((dra<ra_bins[radecbinNumb-1]+0.5*rabinInterval) && (ddec<dec_bins[radecbinNumb-1]+0.5*decbinInterval)){
	 			rabin  = (int) floor(dra/rabinInterval);

 				decbin = (int) floor(ddec/decbinInterval); 

 				// printf("\n%d \t %d", rabin, decbin);

 			
				if((dra < 3./(60.*60.)) && (ddec < 72./(60.*60.))){
				  binned_pairs[rabin][decbin] += 1.0;
				}
				}
			//	}
		}
	}

	printMockAvg_correlation();
	
	return 0;
}

int printMockAvg_correlation(){	
        sprintf(filepath, "%s/randoms20_W1_Nagoya_angularCorrelation.dat", root_dir);

	output = fopen(filepath, "w");

	for(j=0; j<radecbinNumb; j++){
	  for(k=0; k<radecbinNumb; k++)  fprintf(output, "%e \t", binned_pairs[j][k]);
		
		fprintf(output, "\n");
	}
	return 0;
}
