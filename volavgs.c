double unity(double chi){
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

  result = (pow(hiChi, 3.) - pow(loChi, 3.))/3.; // Vol. of sphere when stripped of 4pi steradians.  

  // convert from h^-1 Mpc to Mpc. 
  // result *= pow(h, -3.);

  // full sky. 
  // result *= 4.*pi;

  // VIPERS W1 area. 
  if(fieldFlag == 1)  result *= sqdegs2steradians(W1area);  // printf("\n\nSTERADIANS: %.6lf", sqdegs2steradians(W1area));
  if(fieldFlag == 4)  result *= sqdegs2steradians(W4area);
  
  // Mpc to Gpc
  result *= pow(10., -9.);
  
  return result;
}


double calc_volavg_chi(){
  double vol            =  qromb(&chisq, loChi, hiChi);

  double volavg_chi     =  qromb(&chicubed, loChi, hiChi);

  return volavg_chi/vol;
}


double calc_galavg_chi(){
  double vol            =  qromb(&chisq_nbar, loChi, hiChi);

  double volavg_chi     =  qromb(&chicubed_nbar, loChi, hiChi);

  return volavg_chi/vol;
}


double chiSq_fkpweight(double chi){
  return chi*chi*pow((*pt2nz)(chi)*fkpPk/(1. + (*pt2nz)(chi)*fkpPk), 2.);
}


double calc_volavg_fkpweights(){
  double volavg_fkpweights     =  qromb(&chiSq_fkpweight, loChi, hiChi);

  return volavg_fkpweights;
}
