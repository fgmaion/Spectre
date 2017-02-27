// created 10/02/2017
int calc_clippingweights(){
  if(appliedClippingThreshold >= 1000){
    load_clippingweights();

    return 0;
  }
  

  // ** calculate clipping weights, must have no folding. **//  
    double                       chi;
    double*             cell_weights;
    double*             gal_occupied;
    double*            rand_occupied;
    // double*              cell_chi;

    double     number_occupied = 0.0;

    // cell_chi   = malloc(n0*n1*n2*sizeof(*cell_chi));
    cell_weights  = malloc(n0*n1*n2*sizeof(*cell_weights));

     gal_occupied = malloc(n0*n1*n2*sizeof(*rand_occupied));
    rand_occupied = malloc(n0*n1*n2*sizeof(*rand_occupied));

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;

    for(j=0; j<n0*n1*n2; j++)    gal_occupied[j] = 0.0;
    for(j=0; j<n0*n1*n2; j++)   rand_occupied[j] = 0.0;

    // initialise to no clipping. 
    for(j=0; j<n0*n1*n2; j++)    cell_weights[j] = 1.0;
    //for(j=0; j<n0*n1*n2; j++)    cell_chi[j]     = 0.0;

    // binary array, is surveyed or not.                                                                                                          
    for(j=0; j<rand_number; j++){
      if(rand_accept[j] == true){
        // bound survey volume.                                                                                                                     
        boxlabel                              =                 boxCoordinates(rand_x, rand_y, rand_z, j);

        //cell_chi[boxlabel]                   +=                                               rand_chi[j];                                        

        rand_occupied[boxlabel]               =                                                         1;
      }
    }

    // Number of cells over which we average <1+ delta> rescaling.                                                                                 
    for(j=0; j<n0*n1*n2; j++){
      if(rand_occupied[j] > 0.){
        number_occupied += 1.;

        //cell_chi[j] /= rand_occupied[j];                                                                                                          
      }
    }

    printf("\n\nvolume available: %.2lf (h^-1 Mpc)^3.", CellVolume*number_occupied);

    // assign galaxies using NGP.  overdensity contains the galaxy counts. 
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){ 
	    chi                               =                          interp_comovingDistance(zobs[j]);
    
            boxlabel                          =                    boxCoordinates(xCoor, yCoor, zCoor, j);

	    gal_occupied[boxlabel]            = 1; 

            // Add (1./sampling) in units of nbar. 
            overdensity[boxlabel][0]         +=                       pow((*pt2nz)(chi)*sampling[j], -1.);
      }
    }

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]      /=                                        CellVolume;                             

    // Smooth n/<n>.  True/false flag for enforcing zero mean for the filtered field.                                                                            // Tapering at survey boundaries, due to smoothing with zero padding.                                                        
    Gaussian_filter(clipping_smoothing_radius, 0);
    // for(j=0; j<n0*n1*n2; j++) smooth_overdensity[j][0] = overdensity[j][0];

    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]      -=                   1.0;  // -1 or 0 determined by randoms.    

    // Scale (1 + delta) such that <1+ delta> = 1.; i.e. homogeneous "zero mean constraint"; preserving delta_min = -1. 
    double norm = 0.0;

    for(j=0; j<n0*n1*n2; j++){
      // even if randoms don't fully sample the field then still a volume average over a reasonably representative sample.  
      if(rand_occupied[j] > 0.){
	norm += smooth_overdensity[j][0]; // currently "smooth_overdensity" is really smoothed (1 + delta) 
      }
    }

    norm /= number_occupied;

    printf("\n\nMean renormalisation in clipping weights: %.2lf", norm);

    // rescaling of <1+delta>.  by rescaling (1+d), don't slip below d=-1. 
    for(j=0; j<n0*n1*n2; j++)  smooth_overdensity[j][0]  /=                                         norm;  
    
    // now it's delta. 
    for(j=0; j<n0*n1*n2; j++)  smooth_overdensity[j][0]  -=                                          1.0;
    
    // restore mask after smoothing.                                                                                                               
    for(j=0; j<n0*n1*n2; j++){
      if(rand_occupied[j] == 0.){
        // spuriously empty cells (due to lack of randoms where there should be).  but surveyed volume not reimposed. 
                 overdensity[j][0] = 0.0;
	  smooth_overdensity[j][0] = 0.0;
      }
    }
    
    // double var_threshold;
    for(j=0; j<n0*n1*n2; j++){           
      // var_threshold = appliedClippingThreshold*exp(pow((1700. - cell_chi[j])/400., 2.));

      // if(overdensity[j][0] > appliedClippingThreshold){    
      if(smooth_overdensity[j][0] > appliedClippingThreshold){
	    //	if(smooth_overdensity[j][0] > var_threshold){
	    //** SMOOTH ->> OVERDENSITY ONLY ** //
	    //if(rand_occupied[j] == 0.)  printf("%.4lf \n", smooth_overdensity[j][0]);
	
	    cell_weights[j]                   = (1. + appliedClippingThreshold)/(1. + overdensity[j][0]);

            // clipped volume fraction.
            if(rand_occupied[j] > 0.)  fraction_clipped  += 1.0;
        }
    }
    
    fraction_clipped /= number_occupied;
    
    printf("\n\napplied clipping threshold: %lf, fraction of cells clipped: %lf", appliedClippingThreshold, fraction_clipped);
    
    double mass_frac = 0.0, frac_norm = 0.0;

    if(data_mock_flag == 0){
      // sprintf(filepath, "%s/W1_Spectro_V7_1/W%d_clip_weights_d0_5/mock_%d_d0_%.2lf_256_pk.dat", root_dir, fieldFlag, loopCount, appliedClippingThreshold); 
      sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount); 
    }

    if(data_mock_flag == 1){
      // sprintf(filepath, "%s/W1_Spectro_V7_1/W%d_clip_weights_d0_5/mock_%d_d0_%.2lf_256_pk.dat", root_dir, fieldFlag, loopCount, appliedClippingThreshold);                                                                                 
      sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
    }

    int count = 0;

    output = fopen(filepath, "w");
    
    //  Clipping weights should not be renormalised, they simply multiply the (fkp_weight/sampling) without further renormalisation. 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j]  == true){ 
            boxlabel           = boxCoordinates(xCoor, yCoor, zCoor, j);

            clip_galweight[j]  =                 cell_weights[boxlabel];
        
	    if(clip_galweight[j] > 1.)  printf("\nspurious weight: %.2lf", clip_galweight[j]);

	    frac_norm         +=                                    1.0;

	    mass_frac         +=                 cell_weights[boxlabel];

	    // if(rand_occupied[boxlabel] == 0.){  printf("%d \t %.3lf \n", count, clip_galweight[j]); count +=1;}
	}
        
        else{clip_galweight[j] = 0.0;}
    
        fprintf(output, "%le \n", clip_galweight[j]);
    }

    fclose(output);

    printf("\n\nMean weight: %.3lf", mass_frac/frac_norm);

    /*
    // print x-axis projected weight. 
    double projected_weight[n0][n1];

    sprintf(filepath, "%s/W1_Spectro_V7_3/mock_1_d0_4_256_projectedweight.dat", root_dir);

    output = fopen(filepath, "w");

    for(k=0; k<n0; k++){
      for(j=0; j<n1; j++){
	projected_weight[k][j] = 0.0;

	for(i=0; i<n2; i++){
	  boxlabel = i + j*n2 + k*n2*n1;

	  projected_weight[k][j] += smooth_overdensity[boxlabel][0];
	}

	fprintf(output, "%.3lf \t", projected_weight[k][j]);
      }

      fprintf(output, "\n");
    }
    */
    return 0;
}


int load_clippingweights(){
    int line_no;
    
    // sprintf(filepath, "%s/W1_Spectro_V7_1/W%d_clip_weights_d0_5/mock_%d_d0_%.2lf_256_pk.dat", root_dir, fieldFlag, loopCount, appliedClippingThreshold); 
    // sprintf(filepath, "%s/Data/500s/hod_cube/clipping_weights_d0_%.2lf_256_pk.dat", root_dir, appliedClippingThreshold); 
    
    if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
    if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);      

    inputfile = fopen(filepath, "r");  
    
    ch         = 0;
    line_no    = 0;
    
    do{
        ch = fgetc(inputfile);        
    
        if(ch == '\n')
       	  line_no += 1;
    } while(ch != EOF);
    
    
    rewind(inputfile);
   
    for(j=0; j<line_no; j++)  fscanf(inputfile, "%le \n", &clip_galweight[j]);
    
    fclose(inputfile);
    
    return 0;
}
