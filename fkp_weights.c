int calc_bare_fkpweights(){
  // stripped of (d0 dependent) alpha. 
  
  double nbar, chi;

  bare_fkp_norm = 0.0;

  for(j=0; j<rand_number; j++){
    nbar            =    interp_nz(rand_chi[j]);      // assumes randoms up to rand_number are all accepted.
    bare_fkp_norm  += nbar*pow(rand_weight[j], 2.);   // FKP weights for randoms sets the normalisation
  }

  bare_fkp_norm  = sqrt(bare_fkp_norm);
  
  for(j=0; j<rand_number; j++)  rand_weight[j] /= bare_fkp_norm;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      chi               = interp_comovingDistance(zobs[j]);

      fkp_galweight[j]  = 1./(1. + fkpPk*interp_nz(chi));

      fkp_galweight[j] /= bare_fkp_norm;
    }
  }
  
  printf("\n\nBare FKP norm: %.6lf", bare_fkp_norm);
  
  return 0;
}


double veff_integrand(double chi, void* p){
  double nbar;

  nbar  = interp_nz(chi);

  // return pow(chi, 2.);
  return pow(chi*nbar*fkpPk/(1. + nbar*fkpPk), 2.);
}


int calc_veff(double lo_z, double hi_z){
  // Calculate effective volume when fkp weighting. 
  gsl_function F;

  double veff, error, lor, hir;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  loopCount = 10;
  
  spline_nbar(0); // true <n>? no.
  
  F.function = &veff_integrand;

  lor = interp_comovingDistance(lo_z);
  hir = interp_comovingDistance(hi_z);

  gsl_integration_qags(&F, lor, hir, 0, 1e-7, 1000, w, &veff, &error);

  // veff *= sqdegs2steradians(W1area);

  // h^-1 Gpc -> Gpc
  // veff *= pow(0.673, -3.);

  veff /= pow(10., 9.);

  printf("\n\nBetween the limits: %.1lf < z < %.1lf, effective volume for W1 is: %.5lf (h^{-1} Gpc)^3", lo_z, hi_z, veff*sqdegs2steradians(W1area));
  printf("\n\nBetween the limits: %.1lf < z < %.1lf, effective volume for W4 is: %.5lf (h^{-1} Gpc)^3", lo_z, hi_z, veff*sqdegs2steradians(W4area));

  return 0;
}
