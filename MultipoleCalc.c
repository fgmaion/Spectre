int nosort_MultipoleCalc(regress* inst){
  // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
  // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

  // Linear regression against P_i(k, \mu_i)
  //  L_i  = 0.5*(3.*\mu_i**2 -1.)
  //  P_i  = P_i(k, \mu_i)
  
  double              pk;
  double  gal_shot = 0.0;
  double rand_shot = 0.0;

  // randoms have chi reassigned on per mock basis. 
  for(j=0; j<rand_number; j++)  rand_shot += pow(rand_weight[j], 2.);  // shot noise from randoms.

  rand_shot *= alpha*alpha;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j]/sampling[j], 2.);  // galaxy shot noise, inc. sampling.
  }

  printf("\n\nShot noise: randoms %.4lf, galaxies %.4lf", rand_shot, gal_shot);
  
  // clear arrays. 
  for(k=0; k<KBIN_NO; k++){
    inst->Monopole[k]   = 0.0;
    inst->Quadrupole[k] = 0.0;

    inst->Sum_Pi[k]     = 0.0;
    inst->Sum_PiLi[k]   = 0.0;
  }

  printf("\n\nPerforming multipole calculation to quadrupole order:");

  // To do: #pragma omp parallel for private(Index, dummy, k, j, i, k_x, k_y, k_z, kSq, kmodulus, mu)
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      for(i=0; i<nx; i++){
        Index                        = k*n1*nx + j*nx + i;
        
        H_k[Index][0]               /= inst->kM2[Index];  // Correct mass assignment of randoms; cic = 2, ngp = 1.
        H_k[Index][1]               /= inst->kM2[Index];

        pk                           = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

        // pk                       -= rand_shot;
        // pk                       -=  gal_shot;  // Fit for constant shotnoise of galaxies when clipping

        inst->Sum_Pi[inst->kind[Index]]   += pk;
        inst->Sum_PiLi[inst->kind[Index]] += pk*inst->kLi[Index];
      }
    }
  }
     
  for(j=0; j<KBIN_NO; j++){
    // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.

    // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
    //            (c d)

    // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li

    inst->Monopole[j]    = (1./inst->detA[j])*( inst->Sum_Li2[j]*inst->Sum_Pi[j] - inst->Sum_Li[j]*inst->Sum_PiLi[j]);
    inst->Quadrupole[j]  = (1./inst->detA[j])*(-inst->Sum_Li[j]*inst->Sum_Pi[j] + inst->modes_perbin[j]*inst->Sum_PiLi[j]);
  
    if(log10(inst->detA[j]) > -6.0)  printf("\n%le \t %le \t %le \t %d", inst->mean_modk[j], inst->Monopole[j], inst->Quadrupole[j], inst->modes_perbin[j]);
  }
    
  return 0;
}


int print_multipoles(regress* inst){
  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/pk/d0_1000/W%d/mock_%03d_zlim_%.1lf_%.1lf_Jf_%d.dat", root_dir, fieldFlag, loopCount, lo_zlim, hi_zlim, 2*fold);
  
  output = fopen(filepath, "w");

  for(j=0; j<KBIN_NO; j++){
    if(log10(inst->detA[j]) > -6.0)  fprintf(output, "%le \t %le \t %le \t %d \n", inst->mean_modk[j], inst->Monopole[j], inst->Quadrupole[j], inst->modes_perbin[j]); 
  }
  
  fclose(output);
  
  return 0;
}
