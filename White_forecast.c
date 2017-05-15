// https://arxiv.org/pdf/0810.1518.pdf
double Veff(double nbar){
  double V0, factor;

  V0     =   calc_vol();
  V0    *= pow(10., 9.); // [(h^-1 Mpc)^3]

  factor = nbar*fkpPk/(1. + nbar*fkpPk);
  factor = pow(factor, 2.);
  
  return V0*factor;
}

double F_nm(int n, int m, double fid_f, double fid_b, double fid_sigz2){
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double result, error;

  gsl_function F;

  F.function = &nbar_dV;

  gsl_integration_qags(&F, loChi - 0.1, hiChi + 0.1, 0, 1e-7, 1000, w, &norm, &error);
}

double integrand(int n, int m, double k, double fid_f, double fid_b, double fid_sigz2){
  return 0.5*k*k/pow(2.*pi, 2.)*pair_derivs(n, m, double k, double mu, fid_f, fid_b, fid_sigz2) ; // [dk dmu]
}

double mu_integrand(int n, int m, double k, double fid_f, double fid_b, double fid_sigz2){
  

}

double pair_derivs(int n, int m, double k, double mu, double f, double b, double sigz2){
  return derivs(n, k, mu, f, b, sigz2)*derivs(m, k, mu, f, b, sigz2);
}

// Ordering [fsig8, bsig8, sigp]
double derivs(int n, double k, double mu, double f, double b, double sigz2){
  switch(n){
  case 0:
    return dlnP_dlnf(mu, f, b);
  case 1:
    return dlnP_dlnb(mu, f, b);
  case 2:
    return dlnP_dsigz2(k, mu);
  }
}

double dlnP_dlnf(double mu, double f, double b){
  return mu*mu* dlnP_dlnb(b, f, mu); // eqn. (3)
}

double dlnP_dlnb(double mu, double f, double b){
  return 2./(b + f*mu*mu); // eqn. (3)
}

double dlnP_dsigz2(double k, double mu){
  return -k*k*mu*mu;
}
