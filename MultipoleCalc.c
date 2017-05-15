int nosort_MultipoleCalc(regress* inst, int mock_start){
  // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
  // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

  // Linear regression against P_i(k, \mu_i)
  //  L_i  = 0.5*(3.*\mu_i**2 -1.)
  //  P_i  = P_i(k, \mu_i)
  
  double              pk;
  double  gal_shot = 0.0;
  double rand_shot = 0.0;

  double   Sum_Pi[KBIN_NO];
  double Sum_PiLi[KBIN_NO];
  
  // update with alpha factors. 
  gal_shot  =  bare_gal_shot/alpha;
  rand_shot = bare_rand_shot*alpha;
  
  printf("\n\nShot noise: randoms %.4lf, galaxies %.4lf", rand_shot, gal_shot);
  /*
  if(fold == 0){
    print_nbarshot(mock_start);
  }
  */
  // clear arrays. 
  for(k=0; k<KBIN_NO; k++){
    inst->Monopole[k]   = 0.0;
    inst->Quadrupole[k] = 0.0;
  
    inst->Sum_Pi[k]     = 0.0;
    inst->Sum_PiLi[k]   = 0.0;

    Sum_Pi[k]           = 0.0; // needed for omp reduction. 
    Sum_PiLi[k]         = 0.0;
  }

  walltime("\n\nPerforming multipole calculation to quadrupole order:");

  #pragma omp parallel for reduction(+: Sum_Pi[:KBIN_NO], Sum_PiLi[:KBIN_NO]) private(Index, k, j, i, pk) if(thread ==1)
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      for(i=0; i<nx; i++){
        Index                        = k*n1*nx + j*nx + i;
        
        H_k[Index][0]               /= inst->kM2[Index];  // Correct mass assignment of randoms; cic = 2, ngp = 1.
        H_k[Index][1]               /= inst->kM2[Index];

        pk                           = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

        // pk                          -= rand_shot;
        // pk                          -=  gal_shot;  // Fit for constant shotnoise of galaxies when clipping

        Sum_Pi[inst->kind[Index]]   += pk;
        Sum_PiLi[inst->kind[Index]] += pk*inst->kLi[Index];
      }
    }
  }
     
  for(j=0; j<KBIN_NO; j++){
    // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.

    // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
    //            (c d)

    // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li

    inst->Sum_Pi[j]      = Sum_Pi[j];
    inst->Sum_PiLi[j]    = Sum_PiLi[j];
    
    inst->Monopole[j]    = (1./inst->detA[j])*( inst->Sum_Li2[j]*inst->Sum_Pi[j] - inst->Sum_Li[j]*inst->Sum_PiLi[j]);
    inst->Quadrupole[j]  = (1./inst->detA[j])*(-inst->Sum_Li[j]*inst->Sum_Pi[j] + inst->modes_perbin[j]*inst->Sum_PiLi[j]);
  
    if(log10(inst->detA[j]) > -6.0)  printf("\n%le \t %le \t %le \t %d", inst->mean_modk[j], inst->Monopole[j], inst->Quadrupole[j], inst->modes_perbin[j]);
  }

  walltime("\n\nPerforming multipole calculation to quadrupole order (end):");
  
  return 0;
}

int print_multipoles(regress* inst){
  // sprintf(filepath, "%s/W1_Spectro_V7_5/mocks_v1.7/pk/mask_pk/W%d/mock_%03d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, fieldFlag, loopCount, lo_zlim, hi_zlim, 2*fold);

  if(data_mock_flag == 0) sprintf(filepath, "%s/mocks_v1.7/pk/d0_%d/W%d/mock_%03d_zlim_%.1lf_%.1lf_Jf_%d.dat", outputdir, d0, fieldFlag, loopCount, lo_zlim, hi_zlim, 2*fold);
  if(data_mock_flag == 1) sprintf(filepath, "%s/data_v1.7/pk/d0_%d/W%d/data_zlim_%.1lf_%.1lf_Jf_%d.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim, 2*fold);
  
  output = fopen(filepath, "w");

  for(j=0; j<KBIN_NO; j++){
    if(log10(inst->detA[j]) > -6.0)  fprintf(output, "%le \t %le \t %le \t %d \n", inst->mean_modk[j], inst->Monopole[j], inst->Quadrupole[j], inst->modes_perbin[j]); 
  }
  
  fclose(output);
  
  return 0;
}

int trash_nbarshot_file(int start){
  int    status = 0;
  char   sys_command[200];
  
  if(data_mock_flag == 0) sprintf(filepath, "%s/mocks_v1.7/pk_derivedprops/d0_%d/W%d/nbarshotnoise_mocks_%d_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, start, lo_zlim, hi_zlim);
  if(data_mock_flag == 1) sprintf(filepath, "%s/data_v1.7/pk_derivedprops/d0_%d/W%d/nbarshotnoise_data_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim);

  sprintf(sys_command, "rm -f %s;", filepath);

  status = system(sys_command);

  if(status == 0){
    printf("\n\nTrashed nbar shotnoise file.");
  }

  else{
    printf("\n\nFailed to trash nbar shotnoise file.");
  }
  
  return 0;
}

int print_nbarshot(int start){
  double  gal_shot = 0.0;
  double rand_shot = 0.0;

  gal_shot  =  bare_gal_shot/alpha;  // update with alpha factors.
  rand_shot = bare_rand_shot*alpha;
  
  if(data_mock_flag == 0) sprintf(filepath, "%s/mocks_v1.7/pk_derivedprops/d0_%d/W%d/nbarshotnoise_mocks_%d_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, start, lo_zlim, hi_zlim);
  if(data_mock_flag == 1) sprintf(filepath, "%s/data_v1.7/pk_derivedprops/d0_%d/W%d/nbarshotnoise_data_zlim_%.1lf_%.1lf.dat", outputdir, d0, fieldFlag, lo_zlim, hi_zlim);

  printf("\n\n%s", filepath);
  
  output = fopen(filepath, "a");

  fprintf(output, "%le \t %le \n", gal_shot, rand_shot);

  fclose(output);
  
  return 0;
}
