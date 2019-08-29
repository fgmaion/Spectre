int xiMonopole(double** rr, double** dd, double* Monopole){
    double norm_gg, norm_gr, norm_rr;
    
    norm_gg = computeNorm(Vipers_Num,   Vipers_Num);
    norm_gr = computeNorm(Vipers_Num,  rand_number);
    norm_rr = computeNorm(rand_number, rand_number);

    double ddCollapsed[nlogbins];
    double rrCollapsed[nlogbins];

    for(i=0; i<nlogbins; i++){
        logrbins[i]    = pow(10., zerolog + (i + 0.5)*logbinsz);

        ddCollapsed[i] = 0.0;
        rrCollapsed[i] = 0.0;
    
        for(j=0; j<nlinbins; j++){
            ddCollapsed[i] += dd[i][j];  
            rrCollapsed[i] += rr[i][j];
        }
     
        Monopole[i] = ddCollapsed[i]/(rrCollapsed[i]*norm_gg/norm_rr) - 1.;
    }    
    
    return 0;
}


int landy_szalay(double **gg, double **gr, double **rr, int Vipers_Num, int rand_number, double **xi){
    double norm_gg, norm_gr, norm_rr;
    
    norm_gg = computeNorm(Vipers_Num,   Vipers_Num);
    norm_gr = computeNorm(Vipers_Num,  rand_number);
    norm_rr = computeNorm(rand_number, rand_number);

    for (i=0; i<nlogbins;i++){
        for (j=0; j<nlinbins; j++){
            // Landy-Szalay
	        // xi[i][j] = gg[i][j]/(rr[i][j]*norm_gg/norm_rr) - 2.*(0.25*norm_gr/norm_rr)*gr[i][j]/(rr[i][j]*norm_gg/norm_rr) + 1.;
	        
	        // Basic DD/RR -1 estimator 
	        xi[i][j] = gg[i][j]/(rr[i][j]*norm_gg/norm_rr) - 1.;
        }
    }

    return 0;
}


int assignMemory_xi(){
    point_gals  = (Particle *) realloc(point_gals,  Vipers_Num*sizeof(Particle));
    point_rands = (Particle *) realloc(point_rands, rand_number*sizeof(Particle));

    gg          = malloc(nlogbins*sizeof(double*));
    gr          = malloc(nlogbins*sizeof(double*));
    rr          = malloc(nlogbins*sizeof(double*));
    landy_xi    = malloc(nlogbins*sizeof(double*));
    
    gg_meanr    = malloc(nlogbins*sizeof(double*));
    gg_meanmu   = malloc(nlogbins*sizeof(double*));
     
    gr_meanr    = malloc(nlogbins*sizeof(double*));
    gr_meanmu   = malloc(nlogbins*sizeof(double*));
    
    rr_meanr    = malloc(nlogbins*sizeof(double*));
    rr_meanmu   = malloc(nlogbins*sizeof(double*));
    
    for(j=0; j<nlogbins; j++){
        gg[j] = malloc(nlinbins*sizeof(double));

        gr[j] = malloc(nlinbins*sizeof(double));

        rr[j] = malloc(nlinbins*sizeof(double));
        
        gg_meanr[j]  = malloc(nlinbins*sizeof(double));
        gg_meanmu[j] = malloc(nlinbins*sizeof(double));
        
        gr_meanr[j]  = malloc(nlinbins*sizeof(double));
        gr_meanmu[j] = malloc(nlinbins*sizeof(double));
        
        rr_meanr[j]  = malloc(nlinbins*sizeof(double));
        rr_meanmu[j] = malloc(nlinbins*sizeof(double));
        
        landy_xi[j] = malloc(nlinbins*sizeof(double));

        for(i=0; i<nlinbins; i++){        
            gg[j][i]       = 0.0;
            rr[j][i]       = 0.0;
            gr[j][i]       = 0.0;
        
            gg_meanr[j][i]  = 0.0;
            gg_meanmu[j][i] = 0.0;
            
            gr_meanr[j][i]  = 0.0;
            gr_meanmu[j][i] = 0.0;
        
            rr_meanr[j][i]  = 0.0;
            rr_meanmu[j][i] = 0.0;
        
            landy_xi[j][i] = 0.0;
        }
    }

    xi0      = malloc(nlogbins*sizeof(*xi0));
    xi2      = malloc(nlogbins*sizeof(*xi2));
    xi4      = malloc(nlogbins*sizeof(*xi4));

    logrbins = malloc(nlogbins*sizeof(*logrbins));

    pairsperbin = malloc(nlogbins*sizeof(*pairsperbin));
    
    for(i=0; i<nlogbins; i++){
        xi0[i]  = 0.0; 
        xi2[i]  = 0.0; 
        xi4[i]  = 0.0; 
        
        logrbins[i] = 0.0;
        
        pairsperbin[i] = 0;
    }
    
    return 0;
}


int computeCorrelation_fns(){
    printf("\n\nComputing correlation fns...");
    
    sprintf(surveyType, "lnNormal_10000_basicEstimator_tree_%d", loopCount);
    
    // logarithmic binning in r. 
    zerolog  =                    0.0;
    maxlog   =             log10(50.);
    logbinsz =             log10(1.4);
  
    nlogbins =  (int) ceil((maxlog - zerolog)/logbinsz);

    printf("\nlog bins: %d", nlogbins);

    // linear binning in mu. 
    zerolin  =                   0.00;
    maxlin   =                   1.10;
    linbinsz =                   0.05;

    nlinbins =  (int) ceil((maxlin - zerolin)/linbinsz);
  
    printf("\nlin bins: %d", nlinbins);
    
    assignMemory_xi();
    
    growTrees();
    
    /*
    bruteforce_nonodes(gg, gg_meanr, gg_meanmu, point_gals,  point_gals,  Vipers_Num,  Vipers_Num,  1);
  
    bruteforce_nonodes(gr, gr_meanr, gr_meanmu, point_gals,  point_rands, Vipers_Num,  rand_number, 0);
    
    bruteforce_nonodes(rr, rr_meanr, rr_meanmu, point_rands, point_rands, rand_number, rand_number, 1);
    */
    
    // printf("\n\nCounting RR pairs.");                                                                                                             
    // CountPairs_rMu(rr, rr_meanr, rr_meanmu, randTree, randTree, 1);
    
    // print_rr();
    
    printf("\n\nCounting DD pairs.");
    CountPairs_rMu(gg, gg_meanr, gg_meanmu, galTree,  galTree,  1);

    print_dd();
    
    // printf("\n\nCounting DR pairs.");
    // CountPairs_rMu(gr, gr_meanr, gr_meanmu, galTree, randTree,  0);
    
    // print_dr();
    
    // sprintf(filename, "DD_%s.dat", surveyType);
    // load2d(filename, gg);
    
    // sprintf(filename, "DR_%s.dat", surveyType);
    // load2d(filename, gr);
    
    sprintf(filename, "RR_lnNormal_10000_basicEstimator_tree.dat", surveyType);
    load2d(filename, rr);
    
    // Compute correlation fn. multipoles.
    landy_szalay(gg, gr, rr, Vipers_Num, rand_number, landy_xi);
   
    // print_xi();
 
    // xiMultipoles(rr, landy_xi, xi0, xi2, pairsperbin);

    xiMonopole(rr, gg, xi0);

    print_xiMultipoles();
    
    return 0;
}


int load2d(char filename[], double** array){
  sprintf(filepath, "%s/Data/Today/%s", root_dir, filename);

  inputfile = fopen(filepath, "r");

  for(j=0; j<nlogbins; j++){
    for(i=0; i<nlinbins; i++)  fscanf(inputfile, "%le \t", &array[j][i]);

    fscanf(inputfile, "\n");
  }

  fclose(inputfile);

  return 0;
}


int print_xiMultipoles(){
  sprintf(filepath, "%s/Data/Today/xiMultipoles_%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++)  fprintf(output, "%e \t %e \t %e \t %e \t %d \n", logrbins[j], xi0[j], xi2[j], xi4[j], pairsperbin[j]);

  fclose(output);

  return 0;
}


int print_xi(){
  sprintf(filepath, "%s/Data/Today/xi_%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", landy_xi[j][i]);

    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}


int print_rr(){
  sprintf(filepath, "%s/Data/Today/RR_lnNormal_10000_basicEstimator_tree.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){  
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", rr[j][i]);
    
    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}


int print_dd(){
  sprintf(filepath, "%s/Data/Today/DD_%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", gg[j][i]);

    fprintf(output, "\n");
  }

  fclose(output);
  
return 0;
}


int print_dr(){
  sprintf(filepath, "%s/Data/Today/DR_%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", gr[j][i]);

    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}
