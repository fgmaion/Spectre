int set_fittomean_params(void){
  double maxL_fsig8, maxL_sigv, maxL_bsig8;

  data_mock_flag = 2;
  
  calc_ChiSqs(0, 0);
  
  set_minChiSq();
  
  // calc_onedposteriors(&maxL_fsig8, &maxL_bsig8, &maxL_sigv);
  
  fsigma8       = minX2_fsig8;
  bsigma8       = minX2_bsig8;
  velDispersion = minX2_sigp;

  // fsigma8       =   0.05;
  // bsigma8       = 0.6875;
  // velDispersion = 5.6250;

  // fsigma8       = 0.1437;
  // bsigma8       = 0.6875;
  // velDispersion = 5.6250;
  
  printf("\n\nFit-to-mean params:  %.6lf \t %.6lf \t %.6lf", fsigma8, velDispersion, bsigma8);

  epsilon_pad   = 0.000000;
  alpha_pad     = 1.000000;


  load_data();

  get_datashotnoise();         // assigns to mean_shot.
    
  for(j=0; j<mono_order; j++)  xdata[j] -= mean_shot;

  ctypeskvals_matchup();
  
  calc_ChiSq(fsigma8, bsigma8, velDispersion, epsilon_pad);
  
  return 0;
}

int getmockmean_params(int d0){
  int                            line_no;
  double   fs8, bs8, sp, fs8L, bs8L, spL;
  
  fsigma8       = 0.0;
  velDispersion = 0.0;
  bsigma8       = 0.0;  

  sprintf(filepath, "%s/mocks_v1.7/fsig8/d0_%d/W%d/kmax_%.1lf/mocks_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, ChiSq_kmax, lo_zlim, hi_zlim);

  inputfile = fopen(filepath, "r");

  line_count(inputfile, &line_no);

  // printf("\n\nGetting mock mean parameters.  File: %s.  Input line number: %d", filepath, line_no);
  
  for(j=0; j<line_no; j++){
    fscanf(inputfile, "%lg \t %lg \t %lg \t %lg \t %lg \t %lg \n", &fs8, &sp, &bs8, &fs8L, &spL, &bs8L);

    // printf("\n%d \t %.6lf \t %.6lf \t %.6lf", j, fs8, sp, bs8);
    
    fsigma8       += fs8;
    bsigma8       += bs8;
    velDispersion +=  sp;
  }

  fclose(inputfile);

  fsigma8        /= (double) line_no;
  bsigma8        /= (double) line_no;
  velDispersion  /= (double) line_no;
      
  // printf("\n\nMock Mean:  %.6lf \t %.6lf \t %.6lf", fsigma8, velDispersion, bsigma8);

  epsilon_pad   = 0.000000;
  alpha_pad     = 1.000000;
  
  return 0;
}

int default_params(void){
  fsigma8       = 0.500000;
  velDispersion = 5.625000;
  bsigma8       = 0.652500;
  epsilon_pad   = 0.000000;
  alpha_pad     = 1.000000;
  
  // set_u0();

  return 0;
}

int bestfit_params(void){
  fsigma8       = minX2_fsig8;
  velDispersion = minX2_sigp;
  bsigma8       = minX2_bsig8;
  alpha_pad     = minX2_alpha_pad;
  epsilon_pad   = minX2_epsilon_pad;
  
  return 0;
}

int set_chiSq_intervals(void){
    fsigma8Interval       = (max_fsigma8     - min_fsigma8)/dRes;
    bsigma8Interval       = (max_bsigma8     - min_bsigma8)/dRes;
    sigmaInterval         = (max_velDisperse - min_velDisperse)/dRes;
    A11SqInterval         = (max_A11Sq       - min_A11Sq)/dRes;

    alpha_padInterval     = (max_alpha_pad   - min_alpha_pad)/dRes_ap;
    epsilon_padInterval   = (max_epsilon_pad - min_epsilon_pad)/dRes_ap;
    
    return 0;
}

int kvals_matchup(void){
  double     diff;
  double min_diff;
    
  fftlog_indices = realloc(fftlog_indices, mono_order*sizeof(*fftlog_indices));
    
  for(i=0; i<mono_order; i++){  
    min_diff = pow(10., 99.);
        
    for(j=0; j<allmono_order; j++){  
      diff = fabs(all_kVals[j] - kVals[i]);
            
      if(diff<min_diff){
        min_diff = diff;
              
        fftlog_indices[i]  = j;
      }
    }
  }
  
  return 0;
}

int calc_models(){  
  sprintf(filepath, "%s/models/realspace_%s_sig8_%.3lf_d0_%d_W%d_%.1lf_%.1f_res_%d.cat", models_path, model_flag, camb_sig8, d0, fieldFlag, lo_zlim, hi_zlim, Res);

  printf("\n\nWriting models: %s", filepath);
  
  output = fopen(filepath, "wb");
  
  for(int aa=0; aa<Res; aa++){    
    fsigma8 = min_fsigma8 + fsigma8Interval*aa;

    for(int bb=0; bb<Res; bb++){
      bsigma8 = min_bsigma8 + bsigma8Interval*bb;

      for(int cc=0; cc<Res; cc++){
        velDispersion = min_velDisperse + sigmaInterval*cc;

        for(int dd=0; dd<Res_ap; dd++){
          alpha_pad = min_alpha_pad + alpha_padInterval*dd;
          
          for(int ee=0; ee<Res_ap; ee++){
            epsilon_pad = min_epsilon_pad + epsilon_padInterval*ee;
            
            alpha_pad   = 1.0;
            epsilon_pad = 0.0;

            // printf("\nfsig8, bsig8, sigv: %.4lf \t %.4lf \t %.4lf", fsigma8, bsigma8, velDispersion);
            
            model_compute(aa, bb, cc, dd, ee, 0);  // updates convlmonoCorr, convlquadCorr 
            
            for(j=0; j<allmono_order; j++)  fwrite(&convlmonoCorr->pk[allfftlog_indices[j]][0], sizeof(double), 1, output);
            for(j=0; j<allmono_order; j++)  fwrite(&convlquadCorr->pk[allfftlog_indices[j]][0], sizeof(double), 1, output);

            // printf("\n\nNew model.");

            // for(j=0; j<allmono_order; j++)  printf("\n%.6le \t %.6le", convlmonoCorr->pk[allfftlog_indices[j]][0], convlquadCorr->pk[allfftlog_indices[j]][0]);
          }
        }
      }
    }
  }
  
  fclose(output);
  
  return 0;
}

int calc_ChiSqs(int mockNumber, int print){    
    int ll, mm, nn;

    minChiSq = pow(10., 12.);
    
    if     (data_mock_flag == 0)  load_mock(mockNumber);
    else if(data_mock_flag == 1)  load_data();
    else if(data_mock_flag == 2)  load_meanMultipoles();
    
    else{
      printf("\n\nChi sq. input is invalid.");
    }
    
    // printf("\n\nChi sq. input.");

    if(mull == 0){  // mull results are already shot noise subtracted. 
      if(data_mock_flag == 0){
        for(j=0; j<mono_order; j++)  xdata[j] -= shotnoise_instances[mockNumber - 1]; 
      }
    
      if(data_mock_flag == 1){
        get_datashotnoise();         // assigns to mean_shot. 

        for(j=0; j<mono_order; j++)  xdata[j] -= mean_shot;
      }
    }
    
    // set_oldshotnoise();
    
    // scale_Cov(130);
    
    for(j=0; j<order; j++){
       ydata[j] = 0.0;  // new zero mean, unit variance, decorrelated variables.
      dkdata[j] = 0.0;
      
      gsl_matrix_get_col(col, evec, j);
        
      for(k=0; k<order; k++){
        ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
       dkdata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*kdata[k];	// normalisation problem? 
      }
    }

    // printf("\n\nDecorrelated data (normalisation?): ");
    
    // for(j=0; j<order; j++)  printf("\n%+.4le \t %+.4le \t %+.4le", dkdata[j], ydata[j], sqrt(gsl_vector_get(eval, j)));
    
    // printf("\n\nChi sq. calc:");

    // Data races currently: #pragma omp parallel for private(kk, jj, ii, ll, mm, nn, fsigma8, bsigma8, velDispersion)
    for(jj=0; jj<Res; jj++){    
      fsigma8 = min_fsigma8 + fsigma8Interval*jj;

      for(kk=0; kk<Res; kk++){
        bsigma8 = min_bsigma8 + bsigma8Interval*kk;

        for(ii=0; ii<Res; ii++){
          velDispersion = min_velDisperse + sigmaInterval*ii;
                
          for(ll=0;ll<Res_ap; ll++){
            // alpha_pad = min_alpha_pad + alpha_padInterval*ll;
                
            for(mm=0; mm<Res_ap; mm++){
              // epsilon_pad = min_epsilon_pad + epsilon_padInterval*mm;
                    
              alpha_pad   = 1.0;
              epsilon_pad = 0.0; 

              ytheory_compute(jj, kk, ii, ll, mm);

              ChiSqGrid[jj][kk][ii][ll][mm] = 0.0;
              
              for(nn=0; nn<order; nn++){
                // ChiSqGrid[jj][kk][ii][ll][mm] += pow(xdata[nn] - xtheory[jj][kk][ii][ll][mm][nn], 2.)*pow(gsl_matrix_get(sigma_norm, nn, nn), 2.);
                
                ChiSqGrid[jj][kk][ii][ll][mm] += pow(ydata[nn] - ytheory[jj][kk][ii][ll][mm][nn], 2.)/gsl_vector_get(eval, nn);                

                // ChiSqGrid[jj][kk][ii][ll][mm] *= (1. - (order + 1)/(CatalogNumber - 1.)); // Hartlap et al. correction.
              }

              if(print == 1)  printf("\n%.6lf \t %.6lf \t %.6lf \t %.6lf", fsigma8, bsigma8, velDispersion, ChiSqGrid[jj][kk][ii][ll][mm]);
            }
          }
        }
      }
    }
    
    return  0;
}

int get_ydata(int mockNumber, int data_mock_flag){
  // Set up zero mean, decorrelated, unit variance, variables.

  if     (data_mock_flag == 0)  load_mock(mockNumber);
  else if(data_mock_flag == 1)  load_data();
  else if(data_mock_flag == 2)  load_meanMultipoles();
  else{
    printf("\n\nChi sq. input is invalid.");
  }
  /*
  if(data_mock_flag == 0){
    for(j=0; j<mono_order; j++)  xdata[j] -= shotnoise_instances[mockNumber - 1];
  }

  if(data_mock_flag == 1){
    get_datashotnoise();

    for(j=0; j<mono_order; j++)  xdata[j] -= mean_shot;
  }
  */
  for(j=0; j<order; j++){
    ydata[j] = 0.0;  // new zero mean, unit variance, decorrelated variables.
     
    gsl_matrix_get_col(col, evec, j);

    for(k=0; k<order; k++)  ydata[j] += gsl_vector_get(col, k)*gsl_matrix_get(sigma_norm, k, k)*xdata[k];
  }
  
  return 0;
}

int ctypeskvals_matchup(){
  double     diff;
  double min_diff;

  fftlog_indices = realloc(fftlog_indices, mono_order*sizeof(*fftlog_indices));

  for(i=0; i<mono_order; i++){
    min_diff = pow(10., 99.);

    for(j=0; j<FFTlogRes; j++){
      diff = fabs(mono_config->krvals[j][0] - kVals[i]);

      if(diff<min_diff){
        min_diff = diff;

        fftlog_indices[i]  = j;
      }
    }
  }
  
  return 0;
}

double calc_ChiSq(double dfsigma8, double dbsigma8, double dvelDispersion, double depsilon){     
    int nn;
    
    fsigma8       =       dfsigma8;
    bsigma8       =       dbsigma8;
    velDispersion = dvelDispersion;
    epsilon_pad   =       depsilon;
    alpha_pad     =            1.0;

    // ctypeskvals_matchup();
    
    // get_ydata(10, 0);
    
    model_compute(0, 0, 0, 0, 0, 0);
    
    for(j=0; j<mono_order; j++)  xtheory[0][0][0][0][0][j]              = convlmonoCorr->pk[fftlog_indices[j]][0];
    for(j=0; j<mono_order; j++)  xtheory[0][0][0][0][0][j + mono_order] = convlquadCorr->pk[fftlog_indices[j]][0];
    
    // for(j=0; j<mono_order; j++)  printf("\n%.6le \t %.6le \t %.6le", kVals[j], xtheory[0][0][0][0][0][j], xtheory[0][0][0][0][0][j + mono_order]);
    
    // ytheory_compute(0, 0, 0, 0, 0);  // over ride zeroth model.

    double dChiSq;
    double  ChiSq  = 0.0;

    printf("\n\nChi sq. eval.");
    
    for(nn=0; nn<order; nn++){
      dChiSq = pow(xdata[nn] - xtheory[0][0][0][0][0][nn], 2.)*pow(gsl_matrix_get(sigma_norm, nn, nn), 2.);

      ChiSq += dChiSq;
      
      printf("\n%.6le \t %.6le \t %le \t %le \t %lf", kVals[nn%mono_order],xdata[nn],xtheory[0][0][0][0][0][nn],1./gsl_matrix_get(sigma_norm,nn,nn),dChiSq);

      if(nn+1 == mono_order)  printf("\n");
    }

    // for(nn=0; nn<order; nn++)  ChiSq += pow(ydata[nn] - ytheory[0][0][0][0][0][nn], 2.)/gsl_vector_get(eval, nn);
    
    // printf("\n\nIndependent eval. of chi sq.: %.6le \n\n", ChiSq);

    return     ChiSq;
    // return -ChiSq/2.;
}

int print_model(double dfsigma8, double dbsigma8, double dvelDispersion, double depsilon){
  int nn;

  fsigma8       =       dfsigma8;
  bsigma8       =       dbsigma8;
  velDispersion = dvelDispersion;
  epsilon_pad   =       depsilon;
  alpha_pad     =            1.0;

  model_compute(0, 0, 0, 0, 0, 0);

  
  sprintf(filepath, "%s/maxlikes/mean_maxlikes_model_W%d_%.1lf_%.1lf_d0_%d.dat", outputdir, fieldFlag, lo_zlim, hi_zlim, d0);

  output = fopen(filepath, "w");

  fprintf(output, "## %.6le \t %.6le \t %.6le \t %.6le \n", dfsigma8, dbsigma8, dvelDispersion, depsilon);
  
  for(j=0; j<FFTlogRes; j++){
    if((0.01 < convlmonoCorr->krvals[j][0]) && (convlmonoCorr->krvals[j][0] < 2.)){      
      fprintf(output, "%.6le \t %.6le \t %.6le \n", convlmonoCorr->krvals[j][0], convlmonoCorr->pk[j][0], convlquadCorr->pk[j][0]);
    }
  }

  fclose(output);

  printf("\n\n");
  
  return 0;
}

int set_models(){
  // Read in models precomputed up to 0.8 and assign xtheory to be these values us to ChiSq_kmax only. 
  int    ll, mm, nn;
  double store[Res][Res][Res][Res_ap][Res_ap][all_order];
  
  sprintf(filepath, "%s/models/realspace_%s_sig8_%.3lf_d0_%d_W%d_%.1lf_%.1f_res_%d.cat", models_path, model_flag, camb_sig8, d0, fieldFlag, lo_zlim, hi_zlim, Res);

  printf("\n\nReading models: %s", filepath);
  
  inputfile = fopen(filepath, "rb");

  while((inputfile = fopen(filepath, "r")) == NULL){
    calc_models();
  }
  
  for(jj=0; jj<Res; jj++){ // f sigma8
    for(kk=0; kk<Res; kk++){ // b sigma8
      for(ii=0; ii<Res; ii++){ // sigma_p
        for(ll=0;ll<Res_ap; ll++){ // alpha_pad
          for(mm=0; mm<Res_ap; mm++){ // epsilon_pad
            fread(store[jj][kk][ii][ll][mm], sizeof(double), all_order, inputfile); // load mono and quad; size of all mono order, i.e. no ChiSq_kmax cut.

            for(j=0; j<all_order; j++) xtheory[jj][kk][ii][ll][mm][j] = 0.0;        // Clean

            for(j=0; j<mono_order; j++){
              xtheory[jj][kk][ii][ll][mm][j]              = store[jj][kk][ii][ll][mm][fftlog_indices[j]];
              xtheory[jj][kk][ii][ll][mm][j + mono_order] = store[jj][kk][ii][ll][mm][fftlog_indices[j] + allmono_order];

              // printf("\n%.6le \t %.6le", xtheory[jj][kk][ii][ll][mm][j], xtheory[jj][kk][ii][ll][mm][j + mono_order]);
            }
            
            // printf("\n\nNew model.");           
            // for(j=0; j<all_order; j++)  printf("\n%.6le \t %.6le", xtheory[jj][kk][ii][ll][mm][j], xtheory[jj][kk][ii][ll][mm][j + allmono_order]);
          }
        }
      }
    }
  }

  fclose(inputfile);
    
  return 0;
}

int test_chiSq(){
  // default_params();
  bestfit_params();
    
  model_compute(0, 0, 0, 0, 0, 0);

  sprintf(filepath, "%s/W1_Spectro_V7_4/model.dat", root_dir);
  
  output = fopen(filepath, "w");

  for(j=0; j<mono_order; j++){
    fprintf(output, "%.4le \t %.4le \t %.4le \t %.4le \t %.4le \n", kVals[j], xdata[j], xdata[j + mono_order], xtheory[0][0][0][0][0][j], xtheory[0][0][0][0][0][j + mono_order]);
  }
  
  fclose(output);
  
  return 0;
}
