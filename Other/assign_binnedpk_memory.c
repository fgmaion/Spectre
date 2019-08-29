int set_kbins_logkspacing(regress* inst){
  logk_interval = (logk_max - logk_min)/KBIN_NO;

  for(j=0; j<KBIN_NO; j++)  inst->logk_limits[j] = pow(10., logk_min + j*logk_interval);

  return 0;
}

int set_kbins_linearkspacing(regress* inst){
  double kmedia;

  kmedia = 2.*pi/800.;

  for(j=0; j<KBIN_NO; j++)  inst->logk_limits[j] = j*kmedia; // needs 800 bins.

  return 0;
}

int prep_pkRegression(void){  
  set_kbins_logkspacing(&flat);
  set_kbins_logkspacing(&half);
  set_kbins_logkspacing(&quart);  

  return 0;
}
