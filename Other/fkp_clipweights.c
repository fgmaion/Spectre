int alpha_clipcalc(){
  accepted       =   0;
  daccepted_gals = 0.0;

  //double daccepted_rand = 0.0;
  
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      accepted           += 1;

      daccepted_gals     += 1./sampling[j];

      //daccepted_gals     += clip_galweight[j]/sampling[j]; //  13/05/16.
    }
  }

  //for(j=0; j<rand_number; j++)  daccepted_rand += clip_randweight[j];
  
  alpha = 1.*daccepted_gals/rand_number;

  printf("\nInfo: d0=% 5d;  Total number of galaxies on input: %d, accepted: %d, accepted(weighted): % 6.2lf. (1./alpha): %.9lf;", d0, Vipers_Num, accepted, daccepted_gals, 1./alpha);

  return 0;
}


int calc_bare_fkpclipweights(){
  // stripped of (d0 dependent) alpha.   
  double nbar, chi;

  bare_fkp_norm = 0.0;

  for(j=0; j<rand_number; j++){
    nbar            = (*pt2nz)(rand_chi[j]);      // assumes randoms up to rand_number are all accepted.

    rand_weight[j] *= clip_randweight[j];

    bare_fkp_norm  += nbar*pow(rand_weight[j], 2.);   // FKP weights for randoms sets the normalisation
  }

  bare_fkp_norm  = sqrt(bare_fkp_norm);
  
  for(j=0; j<rand_number; j++)  rand_weight[j] /= bare_fkp_norm;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      chi               = interp_comovingDistance(zobs[j]);

      fkp_galweight[j]  = 1./(1. + fkpPk*(*pt2nz)(chi));

      fkp_galweight[j] /= bare_fkp_norm;
    }
  }
  
  printf("\n\nBare FKP norm: %.6lf", bare_fkp_norm);
  
  return 0;
}
