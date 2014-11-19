int initialisefiberCatGen(int objectNumber){
  fiber_flag   = malloc(objectNumber*sizeof(double));

  shuffle_rows = malloc(2*objectNumber*sizeof(int));

  for(j=0; j<objectNumber; j++){  
    fiber_flag[j]   = 1.;
    
    shuffle_rows[j] = 0;
  }

  return 0;
}


int fiberCollision_cat(int objectNumber, double redshifts[], double xCoors[], double yCoors[], double zCoors[]){
  double dx, dy, dz, r;
  
  initialisefiberCatGen(objectNumber);
  
  for(j=0; j<objectNumber; j++){  
      if((redshifts[j]<0.65) || (redshifts[j]>0.75)){
        //    if((redshifts[j]>0.65) && (redshifts[j]<0.75) && (flag_Nagoya[j] == 1)){
        fiber_flag[j] = 0.;
      }
  }
  
  double inredshiftrange = 0.0;
  
  for(j=0; j<objectNumber; j++) inredshiftrange += fiber_flag[j];  
  
  // drawrows(shuffle_rows, objectNumber);
  
  printf("\n\nGenerating fiber collision catalogue. input objects: %d", objectNumber);
  
  int count = 0;
  int shuffled_i;
  int shuffled_j;
   
  double phi;
  
  for(i=0; i<objectNumber; i++){    
    // shuffled_i = shuffle_rows[i];
  
    if(i%10000==0)  printf("\n%d", i);
    
      if(fiber_flag[i] > 0.1){      
          for(j=i+1; j<objectNumber; j++){
              // shuffled_j = shuffle_rows[j];
          
	          if(fiber_flag[j] > 0.1){
   	              dx            = xCoors[i] - xCoors[j];

	              dy            = yCoors[i] - yCoors[j];
	  
	              dz            = zCoors[i] - zCoors[j]; 

                  // atan2 returns in the range -pi to pi
            	  // phi        = pi + atan2(dy, dx);

            	  r             = pow(pow(dx, 2.) + pow(dy, 2.) + pow(dz, 2.), 0.5);

                  // printf("\n%e \t %e \t %e \t %e", dx, dy, dz, r);

	              //  if((dx<0.0254) && (dy<0.6087)){  
            	  //  if((r>5.) && (r<20.) && (phi < (0.5*pi)*exp(-(r - 5.)/2.))){
	              if(r<10.){
	                 count += 1;
	   
               	     fiber_flag[j] = 0.;

               	     // printf("\n%d \t %e", count, r);
	              }
	         }
        }  
    }
  }
  
  printf("\n\nfraction removed: %f", 1.*count/inredshiftrange);
  
  sprintf(filepath, "%s/mocks_W1_Nagoya_0.65_0.75_modrfiber_v1.2/Selectedrandoms_W1_001_ALLINFO.cat", vipersHOD_dir, loopCount);

  // if(loopCount<10)  sprintf(filepath, "%s/mocks_W1_Nagoya_0.65_0.75_modrfiber_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
  // else              sprintf(filepath, "%s/mocks_W1_Nagoya_0.65_0.75_modrfiber_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

  output = fopen(filepath, "w");

  for(j=0; j<objectNumber; j++){  
      if((fiber_flag[j] == 1.)){
  	      fprintf(output, "%e \t %e \t %e \n", rand_ra[j], rand_dec[j], redshifts[j]);

	// fprintf(output, "%d \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %d \t %lf \t %lf \ \t %lf \t %d \t %d \t %d \t %lf \n", id[j], ra[j], dec[j], zcos[j], zpec[j], zobs[j], zphot[j], M_B[j], type[j], csr[j], sampling[j], sampling35[j], flag_Nagoya[j], flag_SSPOC[j], flag_SSPOC35[j], rand_sel[j]);
      }
  }
  
  fclose(output);
  
  return 0;
}
