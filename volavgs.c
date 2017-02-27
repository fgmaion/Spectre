double invnbar_chisq(double chi){
  return chi*chi/interp_nz(chi);
}

double chisq(double chi){
  return chi*chi;
}

double chicubed(double chi){
  return chi*chi*chi;
}

double chicubed_nbar(double chi){
  return chi*chi*chi*interp_nz(chi);
}

double chisq_nbar(double chi){
  return chi*chi*interp_nz(chi);
}

double calc_volavg_invnbar(){
  // Shot noise correction for a varying background number density.
  // P(k) is the volume avg. therefore so must be the shot noise correction.

  double vol            =  qromb(&chisq, loChi, hiChi);

  double volavg_invnbar =  qromb(&invnbar_chisq, loChi, hiChi);

  // solid angle drops out if nbar is independent of direction.
  return volavg_invnbar/vol;
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
  return chi*chi*pow(interp_nz(chi)*fkpPk/(1. + interp_nz(chi)*fkpPk), 2.);
}


double calc_volavg_fkpweights(){
  double volavg_fkpweights     =  qromb(&chiSq_fkpweight, loChi, hiChi);

  return volavg_fkpweights;
}
