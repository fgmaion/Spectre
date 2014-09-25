int randomGeneration(){
    assign_randmemory(0.001, Vipers_Num);
        
    for(j=0; j<rand_number; j++){
         rand_ra[j] = LowerRAlimit + (UpperRAlimit - LowerRAlimit)*gsl_rng_uniform(gsl_ran_r);    
         
        // uniform in sin(dec). 
        rand_dec[j] = sin(LowerDecLimit*pi/180.) + (sin(UpperDecLimit*pi/180.) - sin(LowerDecLimit*pi/180.))*gsl_rng_uniform(gsl_ran_r);
        rand_dec[j] = (180./pi)*asin(rand_dec[j]);
    
        rand_chi[j] = pow(LowerChiLimit, 3.) + (pow(UpperChiLimit, 3.) - pow(LowerChiLimit, 3.))*gsl_rng_uniform(gsl_ran_r);
        rand_chi[j] = pow(rand_chi[j], 1./3.);
        
        StefanoBasis(rand_number, rand_ra, rand_dec, rand_chi, rand_x, rand_y, rand_z);
        
        StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
        
        for(j=0; j<rand_number; j++){        
            xlabel                  = (int) floor((rand_x[j] - AxisLimsArray[0][0])/CellSize);    
            ylabel                  = (int) floor((rand_y[j] - AxisLimsArray[0][1])/CellSize);
            zlabel                  = (int) floor((rand_z[j] - AxisLimsArray[0][2])/CellSize);
    
            boxlabel                = (int)                 xlabel + n2*ylabel + n2*n1*zlabel;
    
            /*
            if(Cell_SurveyLimitsMask[boxlabel]  == 0){
                Cell_SurveyLimitsMask[boxlabel] += 1;
            }
            */
            
        }
    }
    
    return 0;
}


int assign_randmemory(double alpha, int Vipers_Num){
    rand_ra    = (double *) malloc(rand_number*sizeof(*rand_ra));
    rand_dec   = (double *) malloc(rand_number*sizeof(*rand_dec));
    rand_chi   = (double *) malloc(rand_number*sizeof(*rand_chi));

    rand_x     = (double *) malloc(rand_number*sizeof(*rand_x));
    rand_y     = (double *) malloc(rand_number*sizeof(*rand_y));
    rand_z     = (double *) malloc(rand_number*sizeof(*rand_z));

    rand_redshift = (double *) malloc(rand_number*sizeof(*rand_redshift));
    
    return 0;
}


int loadNagoya_rands(){
  sprintf(filepath, "/disk1/mjw/VIPERS_ValueAddedHOD/randoms20_W1_Nagoya_ra_dec_z.cat");

  inputfile = fopen(filepath, "r");

  ch          = 0;
  rand_number = 0;

  do{
    ch = fgetc(inputfile);
    if(ch == '\n')
      rand_number += 1;
  } while (ch != EOF);

  printf("\n\n%d randoms. up to date", rand_number);

  rewind(inputfile);

  assign_randmemory(0.1, rand_number);
  
  for(j=0; j<rand_number; j++){  
    fscanf(inputfile, "%le \t %le \t %le", &rand_ra[j], &rand_dec[j], &rand_redshift[j]);

    rand_chi[j] = interp_comovingDistance(rand_redshift[j]);
  }

  fclose(inputfile);
  
  for(j=0; j<rand_number; j++){
    rand_ra[j]                *= (pi/180.0);                                 // Converted to radians.                      
    rand_dec[j]               *= (pi/180.0);                                 // Converted to radians.                    

    rand_x[j]              =     rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);
    rand_y[j]              =     rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
    rand_z[j]              = -1.*rand_chi[j]*sin(rand_dec[j]);

    rand_ra[j]                /= (pi/180.0);                                 // Converted to degrees.                      
    rand_dec[j]               /= (pi/180.0);                                 // Converted to degrees.                    
  }

  printf("\n\nRandoms on input");

  printf("\nx max:  %f \t x min:  %f", arrayMax(rand_x, rand_number), arrayMin(rand_x, rand_number));
  printf("\ny max:  %f \t y min:  %f", arrayMax(rand_y, rand_number), arrayMin(rand_y, rand_number));
  printf("\nz max:  %f \t z min:  %f", arrayMax(rand_z, rand_number), arrayMin(rand_z, rand_number));
  printf("\nr max:  %f \t r min %f",   arrayMax(rand_chi, rand_number), arrayMin(rand_chi, rand_number));
  /*
  for(j=0; j<rand_number; j++){
    if((rand_ra[j] == 0.) || (rand_dec[j] == 0.) || (rand_chi[j]==0.))  printf("\nrand: %d", j);
  }  
  */
  StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
  
  printf("\n\nRandoms successfully loaded. + translated");

  printf("\nx max:  %f \t x min:  %f", arrayMax(rand_x, rand_number), arrayMin(rand_x, rand_number));
  printf("\ny max:  %f \t y min:  %f", arrayMax(rand_y, rand_number), arrayMin(rand_y, rand_number));
  printf("\nz max:  %f \t z min:  %f", arrayMax(rand_z, rand_number), arrayMin(rand_z, rand_number));
  // printf("\nr max:  %f \t r min %f",   arrayMax(rand_chi, rand_number), arrayMin(rand_chi, rand_number));
  
  return 0;
}


int loadRand(){
    sprintf(filepath, "%s/Stefano/xyzRan.dat", root_dir);
    
    inputfile     = fopen(filepath, "r");  
    
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

    ch          = 0;
    rand_number = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  rand_number += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    printf("\n\nRandom number: %d", rand_number);
    
    rand_x      = (double *) malloc(rand_number*sizeof(*rand_x));
    rand_y      = (double *) malloc(rand_number*sizeof(*rand_y));
    rand_z      = (double *) malloc(rand_number*sizeof(*rand_z));
    rand_chi    = (double *) malloc(rand_number*sizeof(*rand_chi));
    
    rand_accept = (bool   *) malloc(rand_number*sizeof(*rand_accept));
    
    for(j=0; j<rand_number; j++) rand_accept[j] = false;
    
    
    for(j=0; j<rand_number; j++){
        fscanf(inputfile, "%lf \t %lf \t %lf \n", &rand_x[j], &rand_y[j], &rand_z[j]);
        
        rand_chi[j] = pow(pow(rand_x[j], 2.) + pow(rand_y[j], 2.) + pow(rand_z[j], 2.), 0.5); 
        
        rand_x[j]   = rand_x[j] + stefano_trans_x;
        rand_y[j]   = rand_y[j] + stefano_trans_y;
        rand_z[j]   = rand_z[j] + stefano_trans_z;   
    
        if((LowerChiLimit < rand_chi[j]) && (rand_chi[j] < UpperChiLimit)){
            rand_accept[j] = true;
        }
    
    }
    
    fclose(inputfile);
    
    printf("\n\nRandoms successfully loaded. + translated");

    printf("\nx max:  %f \t x min:  %f", arrayMax(rand_x, rand_number), arrayMin(rand_x, rand_number));
    printf("\ny max:  %f \t y min:  %f", arrayMax(rand_y, rand_number), arrayMin(rand_y, rand_number));
    printf("\nz max:  %f \t z min:  %f", arrayMax(rand_z, rand_number), arrayMin(rand_z, rand_number));
    printf("\nr max:  %f \t r min %f",   arrayMax(rand_chi, rand_number), arrayMin(rand_chi, rand_number));
    
    for(j=0; j<rand_number; j++){        
        xlabel                  = (int) floor((rand_x[j] - AxisLimsArray[0][0])/CellSize);    
        ylabel                  = (int) floor((rand_y[j] - AxisLimsArray[0][1])/CellSize);
        zlabel                  = (int) floor((rand_z[j] - AxisLimsArray[0][2])/CellSize);
    
        boxlabel                = (int)                 xlabel + n2*ylabel + n2*n1*zlabel;
    
        if((LowerChiLimit < rand_chi[j]) && (rand_chi[j] < UpperChiLimit)){
            Cell_SurveyLimitsMask[boxlabel]  +=  1.;
	        accepted_rand                    +=  1.;
	    }
    }
    /*
    for(j=0; j<n0*n1*n2; j++){
        if(Cell_SurveyLimitsMask[j] <=4.){
            Cell_SurveyLimitsMask[j] = 0.0;
        }
    }
    */
    sprintf(filepath, "%s/Data/SpectralDistortion/VIPERS_mask.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<n0*n1*n2; j++){ 
        fprintf(output, "%e \n", Cell_SurveyLimitsMask[j]);
    }
    
    fclose(output);

    /*
    for(j=0; j<rand_number; j++){ 
        if((LowerChiLimit < rand_chi[j]) && (rand_chi[j] < UpperChiLimit)){ 
            accepted_rand += 1;
    
            cic_assign(rand_x[j], rand_y[j], rand_z[j], CellSize, Cell_SurveyLimitsMask, 1.0);
        }
    }
    */
    
    /*
    sprintf(filepath, "%s/Data/RandomsAssigned2Cells_cic.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<n0*n1*n2; j++){  
        if(Cell_SurveyLimitsMask[j] > 0.0){
            fprintf(output, "\n %e", Cell_SurveyLimitsMask[j]);
        }
    }
    
    fclose(output);
    */
    
    printf("\n\nAccepted number of randoms:  %d", accepted_rand);
    
    printf("\n\nAccepted co-ordinates. ");
    printf("\nx max:  %f \t x min:  %f", AcceptedMax(rand_x, rand_accept, rand_number), AcceptedMin(rand_x, rand_accept, rand_number));
    printf("\ny max:  %f \t y min:  %f", AcceptedMax(rand_y, rand_accept, rand_number), AcceptedMin(rand_y, rand_accept, rand_number));
    printf("\nz max:  %f \t z min:  %f", AcceptedMax(rand_z, rand_accept, rand_number), AcceptedMin(rand_z, rand_accept, rand_number));
    printf("\nr max:  %f \t r min %f",   AcceptedMax(rand_chi, rand_accept, rand_number), AcceptedMin(rand_chi, rand_accept, rand_number));
    
    free(rand_x);
    free(rand_y);
    free(rand_z);
    free(rand_chi);
        
    return 0;
} 


int cube_ranGen(double ya, double za, double yb, double zb, double yc, double zc, double yd, double zd, double xlo, double xhi, int rand_number){
    double m1, m2;
    double c1, c2;
    
    m1 = (zd - za)/(yd - ya);
    m2 = (zc - zb)/(yc - yb);
    
    c1 = (za*yd - ya*zd)/(yd - ya);
    c2 = (zb*yc - yb*zc)/(yc - yb);
    
    rand_x     = (double *) malloc(rand_number*sizeof(*rand_x));
    rand_y     = (double *) malloc(rand_number*sizeof(*rand_y));
    rand_z     = (double *) malloc(rand_number*sizeof(*rand_z));

    for(j=0; j<rand_number; j++){
        for(;;){    
            rand_x[j] = (AxisLimsArray[1][0] - AxisLimsArray[0][0])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][0];
            rand_y[j] = (AxisLimsArray[1][1] - AxisLimsArray[0][1])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][1];
            rand_z[j] = (AxisLimsArray[1][2] - AxisLimsArray[0][2])*gsl_rng_uniform(gsl_ran_r) + AxisLimsArray[0][2];
        
            if((rand_z[j] > za) && (rand_z[j] < zd) && (rand_z[j] >= m1*rand_y[j] + c1) && (rand_z[j] >= m2*rand_y[j] + c2) && (rand_x[j]>xlo) && (rand_x[j]<xhi)){
                break;
            }
        }
    }
    
    for(j=0; j<rand_number; j++){        
        xlabel                  = (int) floor((rand_x[j] - AxisLimsArray[0][0])/CellSize);    
        ylabel                  = (int) floor((rand_y[j] - AxisLimsArray[0][1])/CellSize);
        zlabel                  = (int) floor((rand_z[j] - AxisLimsArray[0][2])/CellSize);
    
        boxlabel                = (int)                 xlabel + n2*ylabel + n2*n1*zlabel;

        Cell_SurveyLimitsMask[boxlabel]  = 1.;
    }
    
    free(rand_x);
    free(rand_y);
    free(rand_z);
    
    return 0;
}
