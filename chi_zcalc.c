double HubbleCnst(double z){
  return 100.*pow(Om_v + Om_m*pow(1.+z, 3.), 0.5); // [h km s^-1 Mpc^-1]
}

double z_chi_integrand_series(double z){
  return pow(Om_v + Om_m, -1.5) - 4.5*Om_m*z*pow(Om_v + Om_m, -2.5) + Om_m*z*z*(16.875*Om_m*pow(Om_v + Om_m, -3.5) - 4.5*pow(Om_v + Om_m, -2.5)); // [lightSpeed_kmpersec/100.]
}

double z_chi_integrand(double z, void* p){
  (void) p;
  
  if(log10(z) < -11.0)  return z_chi_integrand_series(z);

  else{
    return 1./sqrt(Om_v + Om_m*pow(1.+ z, 3.) + Om_r*pow(1.+z, 4.) - (Om_tot - 1.)*pow(1.+z, 2.)); // [lightSpeed_kmpersec/100.]
  }
}

double interp_comovingDistance(double z){
  double result;

  splint(z_Array, ComovingDistance_z, z_ComovingDistance_2derivatives, 1000, z, &result);

  return result;  // Returns comoving distance at redshift z in h^-1 Mpc.
}

double interp_inverseComovingDistance(double r){
  double result;

  splint(ComovingDistance_z, z_Array, ComovingDistance_z_2derivatives, 1000, r, &result);

  return result;  // Returns z at comoving distance, [h^-1 Mpc].
}


int chi_zcalc(){
  printf("\nCalculating z - comoving distance.");

  gsl_integration_cquad_workspace* w = gsl_integration_cquad_workspace_alloc(1000); 

  size_t attempts;
  
  gsl_function  F;

  double result, error;

  F.function = &z_chi_integrand;
  
  // R_0*r = integral (c/H_0)*pow(Om_Lambda + Om_m*(1+z)^3 + Om_r*(1+z)^4 - (Om_tot -1)*(1+z)^2, -0.5) dz 
  for(i=1000; i>0; i--){
    z_Array[1000 - i]          =  2.0 + -2.*i/1000.0;
        
    gsl_integration_cquad(&F, 0.0, z_Array[1000-i], 1e-12, 0.0, w, &result, &error, &attempts);
    
    ComovingDistance_z[1000-i] = (lightSpeed_kmpersec/100.)*result;            
  }        

  // Working correctly, tested against Ned Wright Cosmology calculator.
  // First array must be a monotonically increasing function, start from redshift zero rather than redshift 2.0
  spline(z_Array, ComovingDistance_z, 1000, 1.0e31, 1.0e31, z_ComovingDistance_2derivatives);
  spline(ComovingDistance_z, z_Array, 1000, 1.0e31, 1.0e31, ComovingDistance_z_2derivatives);

  // -- Test -- //
  // gsl_integration_qags(&F, 0.0, 1.5, 0, 1e-7, 1000, w, &result, &error);
  // printf("\nComoving distance to redshift 1.5: %le [Mpc] \n\n", h*result);
  
  // -- Survey limits -- //
  // gsl_integration_qag(&F, 0.0, lo_zlim, 1e-10, 0, 2000, GSL_INTEG_GAUSS61, w, &loChi, &error);
  // gsl_integration_qag(&F, 0.0, hi_zlim, 1e-10, 0, 2000, GSL_INTEG_GAUSS61, w, &hiChi, &error);
  
  gsl_integration_cquad(&F, 0.0, lo_zlim, 1e-12, 0.0, w, &loChi, &error, &attempts);
  gsl_integration_cquad(&F, 0.0, hi_zlim, 1e-12, 0.0, w, &hiChi, &error, &attempts);

  loChi *= lightSpeed_kmpersec/100.;
  hiChi *= lightSpeed_kmpersec/100.;

  gsl_integration_cquad_workspace_free(w);

  printf("\nComoving distance to redshift 1.5: %le [Mpc] \n\n", h*interp_comovingDistance(1.5));
  
  // -- Derived -- //
  UniverseAge();

  linearGrowthRate();  // Must match declaration in cosmology_planck2015.h (cosmology_valueaddedmocks.h); This is NOT automatically ensured. //  
    
  return 0;
}
