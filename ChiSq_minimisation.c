int default_params(){
  fsigma8       = 0.14;
  A11Sq         = 1.00;
  velDispersion = 2.62;
  bsigma8       = 0.82;
  epsilon_pad   = 0.00;
  alpha_pad     = 1.00;

  set_u0();

  return 0;
}


int bestfit_params(){
  fsigma8       = minX2_fsig8;
  velDispersion = minX2_sigp;
  bsigma8       = minX2_bsig8;
  alpha_pad     = minX2_alpha_pad;
  epsilon_pad   = minX2_epsilon_pad;
  
  return 0;
}


int set_chiSq_intervals(){
    fsigma8Interval       = (max_fsigma8     - min_fsigma8)/dRes;
    bsigma8Interval       = (max_bsigma8     - min_bsigma8)/dRes;
    sigmaInterval         = (max_velDisperse - min_velDisperse)/dRes;
    A11SqInterval         = (max_A11Sq       - min_A11Sq)/dRes;

    alpha_padInterval     = (max_alpha_pad   - min_alpha_pad)/dRes_ap;
    epsilon_padInterval   = (max_epsilon_pad - min_epsilon_pad)/dRes_ap;
    
    return 0;
}


int kvals_matchup(){
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
  int aa, bb, cc, dd, ee, klo;

  // input_check();
  // prep_dlnPR_dlnk();   

  // printf("\n\nCalculating models: \n");
  
  sprintf(filepath, "%s/W1_Spectro_V7_4/models/realspace_%s_d0_%d_W%d_%.1lf_%.1f_res_%d.cat", root_dir, model_flag, d0, fieldFlag, lo_zlim, hi_zlim, Res);

  output = fopen(filepath, "wb");
  
  for(aa=0; aa<Res; aa++){
    fsigma8 = min_fsigma8 + fsigma8Interval*aa;

    for(bb=0; bb<Res; bb++){
      bsigma8 = min_bsigma8 + bsigma8Interval*bb;

      for(cc=0; cc<Res; cc++){
        velDispersion = min_velDisperse + sigmaInterval*cc;

        for(dd=0; dd<Res_ap; dd++){
          alpha_pad = min_alpha_pad + alpha_padInterval*dd;
          
          for(ee=0; ee<Res_ap; ee++){
            epsilon_pad = min_epsilon_pad + epsilon_padInterval*ee;
            
            alpha_pad   = 1.0;
            epsilon_pad = 0.0;

            // printf("\nfsig8, bsig8, sigv: %.4lf \t %.4lf \t %.4lf", fsigma8, bsigma8, velDispersion);
            
            model_compute(aa, bb, cc, dd, ee);  // updates convlmonoCorr, convlquadCorr 

            for(j=0; j<allmono_order; j++)  fwrite(&convlmonoCorr->pk[fftlog_indices[j]][0], sizeof(double), 1, output);
            for(j=0; j<allmono_order; j++)  fwrite(&convlquadCorr->pk[fftlog_indices[j]][0], sizeof(double), 1, output);
            
            // fwrite(xtheory[aa][bb][cc][dd][ee], sizeof(double),  order,   output);
            // fwrite(ytheory[aa][bb][cc][dd][ee], sizeof(double),  order,   output);
          }
        }
      }
    }
  }
  
  fclose(output);
  
  return 0;
}


int calc_ChiSqs(int mockNumber){    
    int ll, mm, nn;

    minChiSq = pow(10., 12.);
    
    if(data_mock_flag == 0)  load_mock(mockNumber);  // Needs amplitude corrected. 
    else                     load_data();

    // set_meanmultipoles();
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
    
    printf("\n\nChi sq. calc:");

    for(jj=0; jj<Res; jj++){    
      fsigma8 = min_fsigma8 + fsigma8Interval*jj;

      for(kk=0; kk<Res; kk++){
        bsigma8 = min_bsigma8 + bsigma8Interval*kk;

        for(ii=0; ii<Res; ii++){
          velDispersion = min_velDisperse + sigmaInterval*ii;
                
          for(ll=0;ll<Res_ap; ll++){
            alpha_pad = min_alpha_pad + alpha_padInterval*ll;
                
            for(mm=0; mm<Res_ap; mm++){
              epsilon_pad = min_epsilon_pad + epsilon_padInterval*mm;
                    
              alpha_pad   = 1.0;
              epsilon_pad = 0.0; 

              ytheory_compute(jj, kk, ii, ll, mm);

              ChiSqGrid[jj][kk][ii][ll][mm] = 0.0;
              
              for(nn=0; nn<order; nn++){
                ChiSqGrid[jj][kk][ii][ll][mm] += pow(ydata[nn] - ytheory[jj][kk][ii][ll][mm][nn], 2.)/gsl_vector_get(eval, nn);
                // ChiSqGrid[jj][kk][ii][ll][mm] += pow(xdata[nn] - xtheory[jj][kk][ii][ll][mm][nn], 2.)*pow(gsl_matrix_get(sigma_norm, nn, nn), 2.);
              }

              // printf("\n%.2lf \t %.2lf \t %.2lf \t %.2lf", fsigma8, velDispersion, bsigma8,  ChiSqGrid[jj][kk][ii][ll][mm]);
	      
              if(ChiSqGrid[jj][kk][ii][ll][mm] < minChiSq){
                minChiSq = ChiSqGrid[jj][kk][ii][ll][mm];
                    
                minX2_fsig8       = fsigma8;
                minX2_A11Sq       = A11Sq;
                minX2_sigp        = velDispersion;
                minX2_bsig8       = bsigma8;
                minX2_alpha_pad   = alpha_pad;
                minX2_epsilon_pad = epsilon_pad;

                printf("\n%.2lf \t %.2lf \t %.2lf \t %.2lf", fsigma8, velDispersion, bsigma8,  minChiSq);
              }
            }
          }
        }
      }
    }
    
    return  0;
}


int calc_marginalisedposteriors(){
  // calc_fsigma8Posterior();
  // calc_bsigma8Posterior();
  // calc_velDispPosterior();

  // output = fopen(filepath, "w");

  // fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", maxlike_fsig8,      maxlike_sigv,    maxlike_bsig8, ChiSq_expected);
  // fprintf(output, "%.6lf \t %.6lf \t %.6lf \t %.6lf \n", minChiSq_fsigma8, minChiSq_sigma, minChiSq_bsigma8,       minChiSq);

  // fclose(output);

  return 0;
}


int set_models(){
  int ll, mm;

  sprintf(filepath, "%s/W1_Spectro_V7_4/models/realspace_%s_d0_%d_W%d_%.1lf_%.1f_res_%d.cat", root_dir, model_flag, d0, fieldFlag, lo_zlim, hi_zlim, Res);
  
  inputfile = fopen(filepath, "rb");

  if(inputfile == NULL){
    printf("\n\nCalculating models");
    
    calc_models();
  }

  else{
    for(jj=0; jj<Res; jj++){ // f sigma8
      for(kk=0; kk<Res; kk++){ // b sigma8
        for(ii=0; ii<Res; ii++){ // sigma_p
          for(ll=0;ll<Res_ap; ll++){ // alpha_pad
            for(mm=0; mm<Res_ap; mm++){ // epsilon_pad
              fread(xtheory[jj][kk][ii][ll][mm], sizeof(double), all_order,   inputfile); // load mono and quad; size of all mono order, i.e. no ChiSq_kmax cut.  
            }
          }
        }
      }
    }

    fclose(inputfile);
  }
  
  return 0;
}


int cut_xtheory_bykmax(){
  int ll, mm, nn;

  double spare[Res][Res][Res][Res_ap][Res_ap][all_order];
  
  for(jj=0; jj<Res; jj++){ // fsigma8
    for(kk=0; kk<Res; kk++){ // bsigma8
      for(ii=0; ii<Res; ii++){ // velDispersion
        for(ll=0;ll<Res_ap; ll++){ // alpha_pad
          for(mm=0; mm<Res_ap; mm++){ // epsilon_pad

            for(j=0; j<all_order; j++)  spare[jj][kk][ii][ll][mm][j] = 0.0;
            
            for(j=0; j<mono_order; j++){
              spare[jj][kk][ii][ll][mm][j]              = xtheory[jj][kk][ii][ll][mm][fftlog_indices[j]];
              spare[jj][kk][ii][ll][mm][j + mono_order] = xtheory[jj][kk][ii][ll][mm][fftlog_indices[j] + allmono_order]; 
            }
          }
        }
      }
    }
  }

  
  for(jj=0; jj<Res; jj++){ // fsigma8
    for(kk=0; kk<Res; kk++){ // bsigma8
      for(ii=0; ii<Res; ii++){ // velDispersion
        for(ll=0;ll<Res_ap; ll++){ // alpha_pad
          for(mm=0; mm<Res_ap; mm++){ // epsilon_pad

            for(j=0; j<all_order; j++) xtheory[jj][kk][ii][ll][mm][0] = 0.0; // NaN

            for(j=0; j<mono_order; j++){
              xtheory[jj][kk][ii][ll][mm][j]              = spare[jj][kk][ii][ll][mm][j];
              xtheory[jj][kk][ii][ll][mm][j + mono_order] = spare[jj][kk][ii][ll][mm][j + mono_order];
            }
          }
        }
      }
    }
  }
  
  return 0;
}


int test_chiSq(){
  // default_params();
  bestfit_params();
    
  model_compute(0, 0, 0, 0, 0);

  sprintf(filepath, "%s/W1_Spectro_V7_4/model.dat", root_dir);
  
  output = fopen(filepath, "w");

  for(j=0; j<mono_order; j++)  fprintf(output, "%.4le \t %.4le \t %.4le \t %.4le \t %.4le \n", kVals[j], xdata[j], xdata[j + mono_order], xtheory[0][0][0][0][0][j], xtheory[0][0][0][0][0][j + mono_order]);
  
  fclose(output);
  
  return 0;
}
