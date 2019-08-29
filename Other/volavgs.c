double unity(double chi){
  (void) chi;

  return 1.;
}

double invnbar_chisq(double chi){
  return chi*chi/(*pt2nz)(chi);
}

double chisq(double chi){
  return chi*chi;
}

double chicubed(double chi){
  return chi*chi*chi;
}

double chicubed_nbar(double chi){
  return chi*chi*chi*(*pt2nz)(chi);
}

double chisq_nbar(double chi){
  return chi*chi*(*pt2nz)(chi);
}


double calc_vol(){
  double result;

  result = (pow(hiChi, 3.) - pow(loChi, 3.))/3.; // Vol. of sphere when stripped of 4PI steradians.  

  // convert from h^-1 Mpc to Mpc. 
  // result *= pow(h, -3.);

  // full sky. 
  // result *= 4.*pi;

  // VIPERS W1 area. 
  if(fieldFlag == 1)  result *= sqdegs2steradians(W1area);  // printf("\n\nSTERADIANS: %.6lf", sqdegs2steradians(W1area));
  if(fieldFlag == 4)  result *= sqdegs2steradians(W4area);

  // Rota et al.: ~ parent bounds.
  // if(fieldFlag == 1)  result *= sqdegs2steradians(32.8125);  // printf("\n\nSTERADIANS: %.6lf", sqdegs2steradians(W1area));
  // if(fieldFlag == 4)  result *= sqdegs2steradians(17.875);

  // Rota et al.: ignoring gaps in a pointing, the area is 294.4 arcmin^2.  192 pointings to W1, 96 in W4.
  // if(fieldFlag == 1)  result *= sqdegs2steradians(15.701);  // printf("\n\nSTERADIANS: %.6lf", sqdegs2steradians(W1area));
  // if(fieldFlag == 4)  result *= sqdegs2steradians( 7.851);
  
  // Mpc to Gpc
  result *= pow(10., -9.);
  
  return result;
}


double chiSq_fkpweight2(double chi, void* p){
  (void) p;

  double nbar = (*pt2nz)(chi);

  return chi*chi*pow(nbar*fkpPk/(1. + nbar*fkpPk), 2.); // fkpPk normalisation. 
}


int calc_volavg_fkpweights2(){
  // Calculate cumulative nbar and splint its inverse.
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double result, error;

  gsl_function F;

  F.function = &chiSq_fkpweight2;

  // TEST:
  // loopCount = 10;
  // spline_nbar(0);  
  // END TEST. 

  printf("\n\nReassigning <n(z)> for vol. avg. FKP weights calc.");
  
  spline_nbar(1); //  set to parent <n(z)>, i.e. mock average. 
  
  gsl_integration_qags(&F, loChi, hiChi, 0, 1e-7, 1000, w, &result, &error);

  // VIPERS W1 area.
  if(fieldFlag == 1)  result *= sqdegs2steradians(W1area);  // printf("\n\nSTERADIANS: %.6lf", sqdegs2steradians(W1area));
  if(fieldFlag == 4)  result *= sqdegs2steradians(W4area);

  // Mpc to Gpc
  result *= pow(10., -9.);

  // [10^-3 (h^-1 Gpc)^3]
  result *= 1000.; 

  printf("\n\nTegmark volume: %.6lf", result);
  
  return 0;
}
