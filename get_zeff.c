int get_zeff(){
  double norm = 0.0;
  double rand_redshift;
  double weighted_redshift = 0.0; 
  
  spline_nbar(1); // set to parent <n(z)>

  pt2nz = &interp_nz;

  prep_inverseCumulative_nbar();
  
  #pragma omp parallel for reduction(+:norm, weighted_redshift)  private(j, rand_redshift) if(thread == 1)
  for(j=0; j<rand_number; j++){
       rand_chi[j]     = inverse_cumulative_nbar(rand_rng[j]);

    rand_weight[j]     = 1./(1. + (*pt2nz)(rand_chi[j])*fkpPk);

    rand_redshift      = interp_inverseComovingDistance(rand_chi[j]);

    norm              += rand_weight[j];
    
    weighted_redshift += rand_weight[j]*rand_redshift;
  }
  
  z_eff = weighted_redshift/norm;

  printf("\n\nEffective redshift: %.3lf", z_eff);
  
  return 0;
}
