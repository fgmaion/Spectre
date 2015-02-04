int xiMonopole(double** xi, double** mu, double* Monopole){
    int quality_flag;

    for(i=0; i<nlogbins; i++){
        Monopole[i]         = 0.0;
    
        quality_flag        = 1;
    
        for(j=0; j<nlinbins; j++){
            // second equality tests for NaNs in mu[i][j] (due to zero DD[i][j] pairs). NaN is the only number not equal to itself, not respected by all compilers. 
            // if((gr[i][j] <10) || (mu[i][j] != mu[i][j])){
            //    quality_flag = 0;
            // }
    
            // mu[i][j]        = zerolin + (j+0.5)*linbinsz;
    
            Monopole[i]    += xi[i][j]; //*LegendrePolynomials(mu[i][j], 0)*linbinsz;
        }
     
        logrbins[i] = pow(10., zerolog + (i+0.5)*logbinsz);
    
        if(quality_flag == 1){
            meanMono[i]   += Monopole[i];
            
            Mono_suitable_mockCount[i] += 1.;
        }
        
        printf("\n%e \t %e", logrbins[i], Monopole[i]);
    }    
    
    return 0;
}


int xiMonopole_error(double** xi, double** mu, double* Monopole){
    int quality_flag;
    
    for(i=0; i<nlogbins; i++){
        Monopole[i]         = 0.0;
    
        quality_flag        = 1;
    
        for(j=0; j<nlinbins; j++){
            // second equality tests for NaNs in mu[i][j] (due to zero DD[i][j] pairs). NaN is the only number not equal to itself, not respected by all compilers. 
            if((rr[i][j] <10) || (mu[i][j] != mu[i][j])){
               quality_flag = 0;
            }
    
            Monopole[i]    += xi[i][j]*LegendrePolynomials(mu[i][j], 0)*linbinsz;
        }
     
        logrbins[i] = pow(10., zerolog + (i+0.5)*logbinsz);
    
        if(quality_flag == 1){
            meanMono_error[i] += pow(Monopole[i] - meanMono[i], 2.);
        }
    }    
    
    return 0;
}


int xiQuadrupole(double** xi, double** mu, double* Quadrupole){
    int quality_flag;

    for(i=0; i<nlogbins; i++){
        Quadrupole[i]       = 0.0;
        
        quality_flag        =   1;
    
        for(j=0; j<nlinbins; j++){
            // second equality tests for NaNs in mu[i][j] (due to zero pairs). NaN is the only number not equal to itself, not respected by all compilers. 
            // if((rr[i][j] <10) || (mu[i][j] != mu[i][j])){
            //  quality_flag  = 0;
            //}
            
            mu[i][j]        = zerolin + (j+0.5)*linbinsz;
        
            // Quadrupole[i]  += 5.*xi[i][j]*LegendrePolynomials(mu[i][j], 2)*linbinsz;
        
            Quadrupole[i] += 5.*xi[i][j]*LegendrePolynomials(mu[i][j], 2);
        }
     
        logrbins[i] = pow(10., zerolog + (i+0.5)*logbinsz);
    
        if(quality_flag ==1){        
            meanQuad[i]   += Quadrupole[i];
            Quad_suitable_mockCount[i] += 1.;
        }
        
        printf("\n\nQuadrupole.");
        
        printf("\n%e \t %e", logrbins[i], Quadrupole[i]); 
    }    
    
    return 0;
}


int xiQuadrupole_error(double** xi, double** mu, double* Quadrupole){
    int quality_flag;
    
    for(i=0; i<nlogbins; i++){
        Quadrupole[i]       = 0.0;
    
        quality_flag        = 1;
    
        for(j=0; j<nlinbins; j++){
            // second equality tests for NaNs in mu[i][j] (due to zero DD[i][j] pairs). NaN is the only number not equal to itself, not respected by all compilers. 
            if((rr[i][j] <10) || (mu[i][j] != mu[i][j])){
               quality_flag = 0;
            }
    
            Quadrupole[i]  += 5.*xi[i][j]*LegendrePolynomials(mu[i][j], 2)*linbinsz;
        }
     
        logrbins[i] = pow(10., zerolog + (i+0.5)*logbinsz);
    
        if(quality_flag == 1){
            meanQuad_error[i] += pow(Quadrupole[i] - meanQuad[i], 2.);
        }
    }    
    
    return 0;
}


int xiMultipoles(double** rr, double** dd, double** xi, double* Monopole, double* Quadrupole, int* pairsperbin){
    // Least squares fit between theory prediction, xi(r, mu_i) and measured.  (Linear regression).

    // For a matrix A, paramater vector, (x_0, x_2)^T, X and vector B
    // Required to invert AX  = B. 2x2 matrix inversion.
         	
    // A reads as (a b) for a = sum_bins 1, b = sum_bins L_j, c = sum_bins L_j, d = sum_bins L_j**2.
    //            (c d)        

    // and B = (b1, b2)^T for b1 = sum_bins measured Xj, b2 = sum_bins (measured Xj)*Lj

    // Input theory values for xi to test, all good. monopole and quadrupole recovered correctly. 
    
    double detA, mu_j, Lj;

    double sum_Lj, sum_Lj2, sum_Xj, sum_LjXj;
    
    for(i=0; i<nlogbins; i++){
        sum_Xj   = 0.0;
        sum_Lj   = 0.0;
        sum_Lj2  = 0.0;
        sum_LjXj = 0.0;
        
        logrbins[i] = 0.0;
         
        pairsperbin[i] = 0; 
        
        for(j=0; j<nlinbins; j++){
            pairsperbin[i]  +=                       1.;
            
            logrbins[i]     +=           gg_meanr[i][j];
            
            // mu_j          =          gg_meanmu[i][j];
            mu_j             = zerolin + j*linbinsz;
                
            sum_Lj          += LegendrePolynomials(mu_j, 2);
            
            sum_Lj2         += pow(LegendrePolynomials(mu_j, 2), 2.);
                
            sum_Xj          += xi[i][j];
                                
            sum_LjXj        += xi[i][j]*LegendrePolynomials(mu_j, 2);
        }
        
        // logrbins[i]         /= pairsperbin[i];
        logrbins[i]          = pow(10., zerolog + i*logbinsz);
    
        // det = ad - bc.
        detA                 = pairsperbin[i]*sum_Lj2 - sum_Lj*sum_Lj;
     
        // (x_0, x_2)^T = (A^-1)B  = (1./detA)*(d*b1 - b*b2, -c*b1 + a*b2)^T.
     
           Monopole[i]       = (1./detA)*( sum_Lj2*sum_Xj - sum_Lj*sum_LjXj);
         Quadrupole[i]       = (1./detA)*(-sum_Lj*sum_Xj  + pairsperbin[i]*sum_LjXj);
    }    
    
    return 0;
}


int landy_szalay(){
    double norm_gg, norm_gr, norm_rr;
    
    norm_gg = computeNorm(Vipers_Num,   Vipers_Num);
    norm_gr = computeNorm(Vipers_Num,  rand_number);
    norm_rr = computeNorm(rand_number, rand_number);

    printf("\n\n%d \t %d \t %e\n\n", Vipers_Num, rand_number, norm_gg/norm_rr);

    for(i=0; i<nlogbins;i++){
        for(j=0; j<nlinbins; j++){
            // Landy-Szalay
	        // landy_xi[i][j] = gg[i][j]/(rr[i][j]*norm_gg/norm_rr) - 2.*(0.25*norm_gr/norm_rr)*gr[i][j]/(rr[i][j]*norm_gg/norm_rr) + 1.;
	        
	        // Basic DD/RR -1 estimator 
	        // landy_xi[i][j] = gg[i][j]/(rr[i][j]*norm_gg/norm_rr) - 1.;
	        
	        // Basic DD/DR -1 estimator 
	        // landy_xi[i][j] = gg[i][j]/(gr[i][j]*norm_gg/norm_gr) - 1.;
        
            landy_xi[i][j] = gg[i][j];
        
            // printf("%e \t", landy_xi[i][j]);
        }
        
        // printf("\n");
    }

    return 0;
}


int assignGal_xiMemory(){
    point_gals  = (Particle *) realloc(point_gals,  Vipers_Num*sizeof(Particle));

    for(j=0; j<nlogbins; j++){
        for(i=0; i<nlinbins; i++){        
            gg[j][i]       = 0.0;
            gr[j][i]       = 0.0;
        
            gg_meanr[j][i]  = 0.0;
            gg_meanmu[j][i] = 0.0;
            
            gr_meanr[j][i]  = 0.0;
            gr_meanmu[j][i] = 0.0;
                
            landy_xi[j][i]  = 0.0;
            
            mean_xi[j][i]   = 0.0; 
        }
    }

    for(i=0; i<nlogbins; i++){
        xi0[i]  = 0.0; 
        xi2[i]  = 0.0; 
        xi4[i]  = 0.0; 
        
        logrbins[i] = 0.0;
        
        pairsperbin[i] = 0;
    
        meanMono[i] = 0.0;
        meanQuad[i] = 0.0;
        
        Mono_suitable_mockCount[i] = 0.0;
        Quad_suitable_mockCount[i] = 0.0;
        
        meanMono_error[i] = 0.0;
        
        meanQuad_error[i] = 0.0;
    }
    
    return 0;
}


int assignMemory_xi(){
    point_rands = (Particle *) realloc(point_rands, rand_number*sizeof(Particle));

    gg          = malloc(nlogbins*sizeof(double*));
    gr          = malloc(nlogbins*sizeof(double*));
    rr          = malloc(nlogbins*sizeof(double*));
    
    landy_xi    = malloc(nlogbins*sizeof(double*));
    
    // window
    rr_0        = malloc(nlogbins*sizeof(double*));
    rr_2        = malloc(nlogbins*sizeof(double*));
    rr_4        = malloc(nlogbins*sizeof(double*));
    
    dummy_gg    = malloc(nlogbins*sizeof(double*));
    
    gg_meanr    = malloc(nlogbins*sizeof(double*));
    gg_meanmu   = malloc(nlogbins*sizeof(double*));
     
    gr_meanr    = malloc(nlogbins*sizeof(double*));
    gr_meanmu   = malloc(nlogbins*sizeof(double*));
    
    rr_meanr    = malloc(nlogbins*sizeof(double*));
    rr_meanmu   = malloc(nlogbins*sizeof(double*));
    
    mean_xi     = malloc(nlogbins*sizeof(double*));
    
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
        
        landy_xi[j]  = malloc(nlinbins*sizeof(double));
        
        rr_0[j]      = malloc(nlinbins*sizeof(double));
        rr_2[j]      = malloc(nlinbins*sizeof(double));
        rr_4[j]      = malloc(nlinbins*sizeof(double));
        
        dummy_gg[j]  = malloc(nlinbins*sizeof(double));
        
        mean_xi[j]   = malloc(nlinbins*sizeof(double)); 

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
        
            landy_xi[j][i]  = 0.0;
            
            rr_0[j][i]      = 0.0;
            rr_2[j][i]      = 0.0;
            rr_4[j][i]      = 0.0;
            
            mean_xi[j][i]   = 0.0; 
        }
    }

    xi0      = malloc(nlogbins*sizeof(*xi0));
    xi2      = malloc(nlogbins*sizeof(*xi2));
    xi4      = malloc(nlogbins*sizeof(*xi4));

    logrbins = malloc(nlogbins*sizeof(*logrbins));

    pairsperbin = malloc(nlogbins*sizeof(*pairsperbin));
    
    meanMono           = malloc(nlogbins*sizeof(*meanMono));
    meanQuad           = malloc(nlogbins*sizeof(*meanQuad));
    
    Mono_suitable_mockCount = malloc(nlogbins*sizeof(*Mono_suitable_mockCount));
    Quad_suitable_mockCount = malloc(nlogbins*sizeof(*Quad_suitable_mockCount));
    
    meanMono_error     = malloc(nlogbins*sizeof(*meanMono_error));
    meanQuad_error     = malloc(nlogbins*sizeof(*meanQuad_error));
    
    for(i=0; i<nlogbins; i++){
        xi0[i]  = 0.0; 
        xi2[i]  = 0.0; 
        xi4[i]  = 0.0; 
        
        logrbins[i] = 0.0;
        
        pairsperbin[i] = 0;
    
        meanMono[i] = 0.0;
        meanQuad[i] = 0.0;
        
        Mono_suitable_mockCount[i] = 0.0;
        Quad_suitable_mockCount[i] = 0.0;
        
        meanMono_error[i] = 0.0;
        
        meanQuad_error[i] = 0.0;
    }
    
    return 0;
}


int prep_randpairs(double Nrands, int loadrands, int loadRR){
    // poissonSample_homogeneous(Nrands);

    // load_homogeneous_rands(Nrands, loadrands);
    
    assignMemory_xi();
    
    if(loadrands==1){
        grow_randTree();
    }
    
    if(loadRR ==1){
        sprintf(filename, "NFW_profile_RR_lessRands.dat");
        load2d(filename, rr);
    }
    
    printf("\nRand number: %d", rand_number);
    
    // printf("\n\nCounting RR pairs.");                                                                                                             
    // bruteforce_nonodes(rr, rr_meanr, rr_meanmu, point_rands, point_rands, rand_number, rand_number, 1);    
    // CountPairs_rMu(rr, rr_meanr, rr_meanmu, randTree, randTree, 1);
    
    // print_rr();
    
    return 0;
}


int prep_pairwisepdf(){
    pairwisepdf_pairCount = 0;

    // (10/2), 10 mocks, 1/2 from distinct pairs. 
    pairwisepdf_pairs     = realloc(pairwisepdf_pairs, (40/2)*30000*30000*sizeof(double));

    return 0;
}


int print_pairwisepdf(){
    sprintf(filepath, "%s/Data/stacpolly/pairwisepdf_DispersionModel.dat", root_dir);

    output = fopen(filepath, "w");
    
    for(j=0; j<pairwisepdf_pairCount; j++)  fprintf(output, "%e \n", pairwisepdf_pairs[j]);

    fclose(output);

    return 0;
}


int NFW_profile_pairCount(){
    sprintf(surveyType, "NFW_profile_xi_%d", loopCount);

    assignGal_xiMemory();
    
    grow_galTree();

    printf("\n\nCounting DD pairs.");  //
    CountPairs_rMu(rr_0, rr_2, rr_4, gg_meanr, gg_meanmu, galTree,  galTree,  1);

    // print_dd();

    // sprintf(filename, "DD_%s.dat", surveyType);
    // load2d(filename, gg);
    
    // landy_szalay();

    // xiMonopole(landy_xi,   gg_meanmu, xi0);
    
    // xiQuadrupole(landy_xi, gg_meanmu, xi2);
    
    // print_xiMultipoles();
    
    // printf("\n%d \t %d", Vipers_Num, rand_number);
    
    return 0;
}

/*
int print_W2_2D(){
  double     mu;
  double      r;
  double theory;

  sprintf(filepath, "%s/Data/VIPERS_window2/rand_VIPERS_W1_xi_500_mask_0.7_0.8_gridded_theoryExpectation_hex_2D.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j=j+10){  
    for(i=0; i<nlinbins; i++){  
        mu = zerolin + i*linbinsz;
    
         r = pow(10., zerolog + j*logbinsz);
    
         if((r>=1.) && (r<400.)){
            theory = splint_VIPERS_maskMono(r) + LegendrePolynomials(mu, 2)*splint_VIPERS_maskQuad(r);
         
            fprintf(output, "%e \t %e \t %e \n", log10(r), mu, theory);
            // fprintf(output, "%e \t %e \t %e \n", r, mu, rr_0[j][i]);
         }
    }
  }

  fclose(output);

  return 0;
}*/


int randWindow_pairCount(){
    sprintf(surveyType, "maskmultipoles_W1_500s_xi_%.1f_%.1f_hiRes_hex", lo_zlim, hi_zlim);

    assignMemory_xi();
    
    grow_randTree();
    
    printf("\n\nCounting RR pairs.");
    CountPairs_rMu(rr_0, rr_2, rr_4, rr_meanr, rr_meanmu, randTree,  randTree,  1);
    
    // printf("\n\nCounting DR pairs.");
    // CountPairs_rMu(gr, gr_meanr, gr_meanmu, galTree, randTree,  0);

    // print_W2_2D();
    
    // print_dr();
    
    // sprintf(filename, "DD_%s.dat", surveyType);
    // load2d(filename, rr);
    
    // sprintf(filename, "DR_tree_%s.dat", surveyType);
    // load2d(filename, gr);
    
    // landy_szalay();

    xiMonopole(rr_0,   gg_meanmu, xi0);

    // monopole calc. of *weighted* pairs. 
    xiMonopole(rr_2,   gg_meanmu, xi2);
    
    xiMonopole(rr_4,   gg_meanmu, xi4);
    
    // xiQuadrupole(rr_2, gg_meanmu, xi2);
    
    print_xiMultipoles();

    return 0;
}


int computeCorrelation_fns(int mockNumber){
  printf("\n\nComputing correlation fns for %d mocks...", mockNumber);

  for(loopCount=1; loopCount<mockNumber; loopCount++){    
    sprintf(surveyType, "NFW_profile_xi%d", loopCount);
    
    // load_clustered(1);
    
    Vipers_Num = 147637.;
    
    assignGal_xiMemory();
    
    // toyTrees();
    
    // grow_galTree();
    
    // grow_randTree();
    
    // printf("\n\nCounting DD pairs.");
    // bruteforce_nonodes(gg, gg_meanr, gg_meanmu, point_gals,  point_gals,  Vipers_Num,  Vipers_Num,  1);
    // CountPairs_rMu(gg, gg_meanr, gg_meanmu, galTree,  galTree,  1);
    // print_dd();
    
    // printf("\n\nCounting DR pairs.");
    // bruteforce_nonodes(gr, gr_meanr, gr_meanmu, point_gals,  point_rands, Vipers_Num,  rand_number, 0);
    // CountPairs_rMu(gr, gr_meanr, gr_meanmu, galTree, randTree,  0);
    
    // print_dr();
    
    sprintf(filename, "NFW_profile_DD.dat", loopCount);
    load2d(filename, gg);
    
    // sprintf(filename, "DD_%s_meanr.dat", surveyType);
    // load2d(filename, gg_meanr);
    
    // sprintf(filename, "DD_%d_meanmu.dat", loopCount);
    // load2d(filename, gg_meanmu);
    
    // sprintf(filename, "/poissonSampled_clustered_fog_500_DR/DR_poissonSampled_clustered_fog_500_%d.dat", loopCount);
    // load2d(filename, gr);
    
    // Compute correlation fn. multipoles.
    landy_szalay();

    xiMonopole(landy_xi, gg_meanmu, xi0);
    
    // xiQuadrupole(landy_xi, gg_meanmu, xi2);
    
    // xiMultipoles(rr, gg, landy_xi, xi0, xi2, pairsperbin);
  }
  /*  
  for(j=0; j<nlogbins; j++)  meanMono[j] /= Mono_suitable_mockCount[j];
  for(j=0; j<nlogbins; j++)  meanQuad[j] /= Quad_suitable_mockCount[j];
    
    
  printf("\n\nComputing correlation fns errors...");
  
  for(loopCount=1; loopCount<mockNumber; loopCount++){        
    load_clustered(0);
    
    sprintf(filename, "DD_%d.dat", loopCount);
    load2d(filename, gg);
        
    sprintf(filename, "DD_%d_meanmu.dat", loopCount);
    load2d(filename, gg_meanmu);
        
    // Compute correlation fn. multipoles.
    landy_szalay();
    
    xiMonopole_error(landy_xi, gg_meanmu, xi0);
    
    xiQuadrupole_error(landy_xi, gg_meanmu, xi2);
  }
  
  for(j=0; j<nlogbins; j++)  meanMono_error[j] /= (Mono_suitable_mockCount[j] - 1.);
  for(j=0; j<nlogbins; j++)  meanMono_error[j]  = sqrt(meanMono_error[j]);
  
  for(j=0; j<nlogbins; j++)  meanQuad_error[j] /= (Quad_suitable_mockCount[j] - 1.);
  for(j=0; j<nlogbins; j++)  meanQuad_error[j]  = sqrt(meanQuad_error[j]);
  */
    
  print_xiMultipoles();
  
  // testNode_separations();
    
  // test_particleSeparations(point_gals[10], point_gals[11]);
    
  return 0;
}


int testNode_separations(){
    Node* atoyNode;
    Node* btoyNode;
    
    atoyNode = Create_toyChildNode(75., 85., 95., 105., 100., 120.);
    
    btoyNode = Create_toyChildNode(15., 15., 15., 20., 20., 20.);
    
    printf("\n%e", pow(10., log10_minimum_modDisplacementBetweenNodes(atoyNode, btoyNode)));
    printf("\n%e", pow(10., log10_maximum_modDisplacementBetweenNodes(atoyNode, btoyNode)));

    return 0;
}


int test_particleSeparations(Particle a, Particle b){
    printf("\n\nParticle a: %.2e \t %.2e \t %.2e", a.x[0], a.x[1], a.x[2]);
    printf("\n\nParticle b: %.2e \t %.2e \t %.2e", b.x[0], b.x[1], b.x[2]);
    
    printf("\n\na to b separation: %.3e", pow(10., log10_particleSeparation(a, b)));

    printf("\n\na to b mu: %.3e", pair_zmu(a, b));

    return 0;
}


int load2d(char filename[], double** array){
    sprintf(filepath, "%s/Data/stacpolly/%s", root_dir, filename);

    inputfile = fopen(filepath, "r");

    for(j=0; j<nlogbins; j++){
      for(i=0; i<nlinbins; i++)  fscanf(inputfile, "%le \t", &array[j][i]);

      fscanf(inputfile, "\n");
    }

    fclose(inputfile);

    return 0;
}


int print_xiMultipoles(){
    sprintf(filepath, "%s/Data/500s/%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){  
      if(logrbins[j] > 0.1){
        fprintf(output, "%e \t %e \t %e \t %e \n", logrbins[j], xi0[j], xi2[j], xi4[j]);
      }
  }
  
  fclose(output);

  return 0;
}


int print_rr(){
  sprintf(filepath, "%s/Data/stacpolly/NFW_profile_RR_lessRands.dat", root_dir);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){  
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", rr[j][i]);
    
    fprintf(output, "\n");
  }

  fclose(output);
  
  // print_rr_meanr();
    
  // print_rr_meanmu();

  return 0;
}



int print_rr_meanr(){
  sprintf(filepath, "%s/Data/lnnorm_anisotropic/RR_homogeneous_basicEstimator_tree_meanr_newMany.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){  
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", rr_meanr[j][i]);
    
    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}



int print_rr_meanmu(){
  sprintf(filepath, "%s/Data/lnnorm_anisotropic/RR_homogeneous_basicEstimator_tree_meanmu_newMany.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){  
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", rr_meanmu[j][i]);
    
    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}


int print_dr(){
  sprintf(filepath, "%s/Data/stacpolly/DR_tree_%s.dat", root_dir, surveyType);

  output = fopen(filepath, "w");

  for(j=0; j<nlogbins; j++){
    for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", gr[j][i]);

    fprintf(output, "\n");
  }

  fclose(output);

  return 0;
}
