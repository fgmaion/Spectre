int ln_PR(double ln_k){
  // Given ln_k, return ln(P_R(k)), where P_R is the real space power spectrum. 
  double k;

  k = exp(ln_k);

  return log( (*pt2Pk)(k));
}


// Numerical evaluation of d ln(P_R) / d ln k with gsl.
int prep_dlnPR_dlnk(){
  gsl_function F;

  // Sanity check, d ln(P_R) / d ln k = 0.96 for large scale power law from inflation.

  double result, abserr, k, ln_k;

  spline_lnk    = malloc(2000*sizeof(*spline_lnk));
  dlnPR_dlnk    = malloc(2000*sizeof(*dlnPR_dlnk));
  dlnPR_dlnk_2D = malloc(2000*sizeof(dlnPR_dlnk_2D));

  F.function = &ln_PR;
  F.params   =      0;

  for(j=0; j<2000; j++){
    k = pow(10., -3. + 3.*j/2000.);

    ln_k = log(k);

    spline_lnk[j] = ln_k;

    gsl_deriv_central(&F, ln_k, 1e-8, &result, &abserr);
  
    dlnPR_dlnk[j] = result; 

    // printf ("%.10lf \t %.10lf +/- %.10lf\n", k, dlnPR_dlnk[j], abserr);
  }

  spline(spline_lnk, dlnPR_dlnk, 2000, 1.0e31, 1.0e31, dlnPR_dlnk_2D);

  return 0;
}


double splint_dlnPR_dlnk(double k){
  double Interim, lnk;
  
  lnk = log(k);

  splint(spline_lnk, dlnPR_dlnk, dlnPR_dlnk_2D, 2000, lnk, &Interim);

  return Interim;
}


// Assumes Kaiser-Lorentz model. 
double dP2_dlnk(double k, double beta, double sigma){
  double ks, PR, P0, P2;
  double M0, M2, M4, M6;

  // kappa in Cole et al. notation.                                                                          
  ks    = k*sigma;

  PR    = (*pt2Pk)(k);

  P0    = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2    = PR*kaiserLorentz_multipole(ks, beta, 2);

  M0    = muOrderZero(ks);
  M2    =  muOrderTwo(ks);
  M4    = muOrderFour(ks);  

  return P2*splint_dlnPR_dlnk(k) + (5./2.)*PR*(M0 -3.*(3.-2.*beta)*M2 -5.*beta*(6.-beta)*M4 -21.*beta*beta*M6 +2.*pow(1.+beta, 2.)/(1.+ks*ks/2.));
}


double dP0_dlnk(double k, double beta, double sigma){
  double ks, PR, P0, P2;
  double M0, M2, M4, M6;

  // kappa in Cole et al. notation.                                                                          
  ks    = k*sigma;

  PR    = (*pt2Pk)(k);

  P0    = PR*kaiserLorentz_multipole(ks, beta, 0);
  P2    = PR*kaiserLorentz_multipole(ks, beta, 2);

  M0    = muOrderZero(ks);
  M2    =  muOrderTwo(ks);
  M4    = muOrderFour(ks);

  return P0*splint_dlnPR_dlnk(k) +PR*(pow(1.+beta,2.)/(1.+ks*ks/2.) -M0 -6.*beta*M2 -5.*beta*beta*M4);
}


double AP_P0(double k_fid, double beta, double sigma, double epsilon, double local_alpha){
  // Alcock-Paczynski correction to the monopole, epsilon is the "warping" parameter.
  // See Padmanabhan and White. 

  double PR, P0, P2, kmod_tr, Jacobian;

  kmod_tr  = k_fid/local_alpha;

  Jacobian = pow(local_alpha, 3.);

  PR       = (*pt2Pk)(kmod_tr);

  P0       = PR*kaiserLorentz_multipole(kmod_tr*sigma, beta, 0);
  P2       = PR*kaiserLorentz_multipole(kmod_tr*sigma, beta, 2);

  return (P0 -(2./5.)*epsilon*dP2_dlnk(kmod_tr,beta,sigma) -(6./5.)*epsilon*P2)/Jacobian;
}


double AP_P2(double k_fid, double beta, double sigma, double epsilon, double local_alpha){
  double PR, P0, P2, kmod_tr, Jacobian;
  
  kmod_tr   = k_fid/local_alpha;

  Jacobian  = pow(local_alpha, 3.);

  PR        = (*pt2Pk)(kmod_tr);

  P0        = PR*kaiserLorentz_multipole(kmod_tr*sigma, beta, 0);
  P2        = PR*kaiserLorentz_multipole(kmod_tr*sigma, beta, 2);

  return ((1.-6.*epsilon/7.)*P2 -(4./7.)*epsilon*dP2_dlnk(kmod_tr,beta,sigma) -2.*epsilon*dP0_dlnk(kmod_tr,beta,sigma))/Jacobian; 
}


double KaiserLorentz_apMultipoles(double beta, double F, double f_perp, double k, double sigma, int MonoQuad){
  // n=4 limit of eqn. (A9) of Ballinger, reduces to an effective Kaiser factor -> can use Kaiser Lorentz                   
  // multipoles.                                                                                                             
  double beta_eff;

  beta_eff = (beta + 1.)/(F*F) - 1.;

  return kaiser_multipole(k, beta_eff, MonoQuad)*(*pt2Pk)(k)/(F*pow(f_perp, 3.+ 4.));
}


int AP_correction(){
  double f_perp, f_para, local_alpha, epsilon, sigma, beta, k, F;
   
  sigma  = 0.00;
  beta   = 0.00;

  f_para = 1.00;  
  f_perp = 1.10;
  
  // constraint for alpha = 1.
  // f_perp = pow(f_para, -0.5);

  // Flattening factor.
  F      = f_para/f_perp;

  printf("\n\nFlattening factor: %.4lf", F);

  local_alpha   = pow(pow(f_perp, 2.)*f_para, 1./3.);

  epsilon       = pow(f_para/f_perp, 1./3.) - 1.;
  
  printf("\n\nPadmanabhan & White parameters, alpha: %.4lf \t epsilon: %.4lf\n\n", local_alpha, epsilon);

  // Simple n=4 power law, corresponds to only an effecive Kaiser factor. See (A9) of Ballinger.  Otherwise, feeds 
  // off usual inputLinearPk() etc.                                                                                              
  // pt2Pk   = &powerlaw;

  prep_dlnPR_dlnk();
  
  int sampling = 200;
  
  sprintf(filepath, "%s/W1_Spectro_V7_0/quadrupole_padmanabhan.dat", root_dir);

  output = fopen(filepath, "w");

  for(jj=0; jj<sampling; jj++){
    k = pow(10., -2. + 3.*jj/sampling);    

    if(k>0.3)  break;

    fprintf(output, "%.4le \t %.4le \n", k, AP_P2(k,beta,f_para*sigma,epsilon,local_alpha));
  
    // printf("%.4le \t %.4le \n", k, KaiserLorentz_apMultipoles(beta,F,f_perp,k,sigma,2), AP_P2(k,beta,sigma,epsilon,local_alpha));
  }

  // fclose(output);
  
  // test_mulitpole_calc(1.0, 1.0, beta, sigma);
  
  test_mulitpole_calc(f_perp, f_para, beta, sigma);
  
  return 0;
}


double KaiserLorentz_AP(double f_perp, double f_para, double kmod, double mu, double beta, double sigma_prime){
  // eqn. (A8) of Ballinger. Returns P'(k) with appropriate mapping from (kmod, mu) - mode in assumed cosmology,
  // to that in the true cosmology.                             

  double F, J, arg, ap_aniso, eff_kaiser, damping;

  // Flattening factor.                                                                                          
  F = f_para/f_perp;

  // Jacobian = alpha (Xu et al., 2013)                                                                                       
  J = f_perp*f_perp*f_para;

  // Power-spectrum eval.                                                                                        
  arg = (kmod/f_perp)*sqrt(1. + mu*mu*(pow(F, -2.) -1.));

  // AP-anisotropy factor. No RSD dependence.                                                                    
  ap_aniso = pow(1. + mu*mu*(pow(F, -2.) -1.), -2.);

  // Effective Kaiser factor, RSD dependence.                                                                    
  eff_kaiser = pow(1. + mu*mu*((1.+beta)*pow(F, -2.) - 1.), 2.);

  // Effective damping.                                                                                          
  damping = 1./(1. + 0.5*pow(kmod*mu*sigma_prime, 2.));

  return (*pt2Pk)(arg)*ap_aniso*eff_kaiser*damping/J;
}


int test_mulitpole_calc(double f_perp, double f_para, double beta, double sigma){
  // Alcock-Paczynski effect yields a new P^{obs}(k_fid) which differs in amplitude, by (1/f_perp^2)*(1/f_para) factor, 
  // eqn (A5) of Ballinger and otherwise is the old power spectrum evaluated at a different mode, given by factors f_perp 
  // and f_para.
  
  double pk;

  // Ballinger amplitude term for P(k).
  double Jacobian;

  Jacobian = pow(f_perp*f_perp*f_para, -1.); 

  // Unprimed variables - corresponding to those values in the true cosmology.  
  double  noprime_kx, noprime_ky, noprime_kz, noprime_kSq, noprime_kmod, noprime_mu;

  polar_pkcount = 0;

  int m0, m1, m2;

  prep_pkRegression(-2., log10(modkMax), kbin_no);

  assign2DPkMemory();

  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      for(i=0; i<n2; i++){
	    m0 = k;
	    m1 = j;
	    m2 = i;

    	if(m2>n2/2)  m2                   -= n2;
	    if(m1>n1/2)  m1                   -= n1;
	    if(m0>n0/2)  m0                   -= n0;

    	// Assumed geometry, correspond to primed quantities in Ballinger.
	    k_x                                = kIntervalx*m2;
	    k_y                                = kIntervaly*m1;
	    k_z                                = kIntervalz*m0;

    	Index                              = k*n1*n2 + j*n2 + i;

    	// Primed, 
	    kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

    	kmodulus                           = pow(kSq, 0.5);
    
    	mu                                 = k_z/kmodulus;
	    if(kmodulus < 0.000001)       mu   = 0.0;

    	// For mode assumed cosmology k_x, calculate corresponding mode in the true cosmology. 
	    // true                            // fiducial
	    noprime_kx                         = k_x/f_perp;
	    noprime_ky                         = k_y/f_perp;
        noprime_kz                         = k_z/f_para;

    	noprime_kSq                        = pow(noprime_kx, 2.) + pow(noprime_ky, 2.) + pow(noprime_kz, 2.);

    	noprime_kmod                       = pow(noprime_kSq, 0.5);

	    noprime_mu                         = noprime_kz/noprime_kmod;

	    // Observed pk in the presence of the AP effect. 
	    // pk                              = (*pt2Pk)(noprime_kmod);

	    // eqn. (A8) of Ballinger. Analytic solution for Kaiser-Lorentz model p(k) in presence of AP effect.
	    pk                                 = KaiserLorentz_AP(f_perp, f_para, kmodulus, mu, beta, sigma);

	    // Negative quadrupole test.
	    // pk                              = (2. - 3.*0.5*(3*mu*mu - 1.));
  
	    // Multipoles calculated for primed (observed) variables
        polar_pk[polar_pkcount][0]         = kmodulus;
	    polar_pk[polar_pkcount][1]         = fabs(mu);
	    polar_pk[polar_pkcount][2]         =       pk;

	    twodim_pk[polar_pkcount][0]        = fabs(k_z);
	    twodim_pk[polar_pkcount][1]        = pow(k_y*k_y + k_x*k_x, 0.5);
	    twodim_pk[polar_pkcount][2]        = pk;

	    polar_pkcount                     += 1;
      }
    }
  }

  sprintf(filepath, "%s/W1_Spectro_V7_0/AP_Multipoles_Ballinger_pk_prediction_%.2lf_%.2lf_%.2lf.dat", root_dir, f_perp, f_para, beta);

  MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, polar_pkcount, filepath, 0.0, 1.0, 1);

  Cartesian2Dpk(polar_pkcount);
  
  return 0;
}

