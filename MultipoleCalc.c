int nosort_MultipoleCalc(){
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
  for(k=0; k<kbin_no; k++){
    Monopole[k]       = 0.0;
    Quadrupole[k]     = 0.0;
    Sum_Pi[k]         = 0.0;
    Sum_PiLi[k]       = 0.0;
  }

  printf("\n\nPerforming multipole calculation to Quadrupole order:");
  
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      // for(i=0; i<(n2/2+1); i++){
      for(i=0; i<n2; i++){
        // Index                        = k*n1*(n2/2+1) + j*(n2/2+1) + i;
        Index                        = k*n1*n2 + j*n2 + i;
        
        H_k[Index][0]               /= kM2[Index];  // Correct mass assignment of randoms; cic = 2, ngp = 1.
        H_k[Index][1]               /= kM2[Index];

        pk                           = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

        // pk                       -= rand_shot;
        // pk                       -=  gal_shot;  // Fit for constant shotnoise of galaxies when clipping

        Sum_Pi[kind[Index]]         += pk;
        Sum_PiLi[kind[Index]]       += pk*kLi[Index];
      }
    }
  }
     
  for(j=0; j<kbin_no; j++){
    // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.

    // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
    //            (c d)

    // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li

    Monopole[j]    = (1./detA[j])*( Sum_Li2[j]*Sum_Pi[j] - Sum_Li[j]*Sum_PiLi[j]);
    Quadrupole[j]  = (1./detA[j])*( -Sum_Li[j]*Sum_Pi[j] + modes_perbin[j]*Sum_PiLi[j]);
  
    if(log10(detA[j]) > -6.0)  printf("\n%le \t %le \t %le \t %d", mean_modk[j], Monopole[j], Quadrupole[j], modes_perbin[j]);
  }
    
  return 0;
}


int print_multipoles(char filepath, int file_output){
  output = fopen(filepath, "w");

  // if(log10(detA[j]) > -6.0)  fprintf(output, "%le \t %le \t %le \t %d \n", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin[j]); 

  fclose(output);
  
  return 0;
}
