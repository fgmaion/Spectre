int PkCorrections(int Index, double k_x, double k_y, double k_z, double rand_shot, double gal_shot, double* p_pk, double* p_kmod, double* p_mu){
  double kmodulus, mu, WindowFunc, kSq;

  kSq = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

  kmodulus = pow(kSq, 0.5);

  mu                                 = k_z/kmodulus;
  // if(kmodulus < 0.000001)       mu   = 0.0;     Assume ignoring mean.

  WindowFunc                         = gsl_sf_sinc(k_x*0.5/xNyquistWaveNumber); // Computes \sinc(x) = \sin(\pi x) / (\pi x).
  WindowFunc                        *= gsl_sf_sinc(k_y*0.5/yNyquistWaveNumber);
  WindowFunc                        *= gsl_sf_sinc(k_z*0.5/zNyquistWaveNumber);

  H_k[Index][0]                      = H_k[Index][0]/pow(WindowFunc, 2.); // gals - normalised rands.
  H_k[Index][1]                      = H_k[Index][1]/pow(WindowFunc, 2.); // If smoothing, do not correct mass assignment.

  // Correct mass assignment of randoms; cic = 2, ngp = 1.
  *p_pk                              = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

  // *p_pk                          -= rand_shot;
  // pk                             -=  gal_shot;  // Fit for constant shotnoise of galaxies when clipping

  *p_mu                              = mu;
  *p_kmod                            = kmodulus;

  return 0;
}


int oldprep_pkRegression(){
  polar_pk             = (double **) malloc(n0*n1*n2*sizeof(double*));

  for(j=0; j<n0*n1*n2; j++)  polar_pk[j]    = (double *)  malloc(3*sizeof(double));  // cols of mod k, mu, pk.

  return 0;
}


int oldnosort_MultipoleCalc(int modBinNumb, double mean_modBin[], double Monopole[], double Quadrupole[], double** Array, int modeCount, char filepath[], double mu_lolimit, double mu_hilimit, int fileOutput){
  
  // P(k, mu_i) = P_0(k) + P_2(k)*(3.*mu_i*mu_i - 1.)/2
  // Least squares fit between theory prediction, P(k, mu_i) and measured.  (Linear regression).

  printf("\n\nPerforming multipole calculation. (to Quadrupole order)");

  double Li, Pi;
  double Sum_Li[modBinNumb];
  double Sum_Li2[modBinNumb];

  double Sum_Pi[modBinNumb];
  double Sum_PiLi[modBinNumb];

  double detA[modBinNumb];  // det = ad - bc.

  int    modes_perbin[modBinNumb];

  double logk_interval;

  logk_interval        = (logk_max - logk_min)/modBinNumb;

  // assign log(k) binning for P(k) and assign memory for binning. (logk_min, logk_max, # bins).
  // prep_pkbinning(-2., log10(modkMax), kbin_no);

  for(k=0; k<modBinNumb; k++){
    mean_modBin[k]    = 0.0;

    Monopole[k]       = 0.0;
    Quadrupole[k]     = 0.0;
    modes_perbin[k]   =   0;

    Sum_Li[k]         = 0.0;
    Sum_Li2[k]        = 0.0;

    Sum_Pi[k]         = 0.0;
    Sum_PiLi[k]       = 0.0;
  }

  Li = Pi = 0.0;
  
  for(j=0; j<polar_pkcount; j++){
    // Linear regression against P_i(k, \mu_i)
    //  L_i  = 0.5*(3.*\mu_i**2 -1.)
    //  P_i  = P_i(k, \mu_i)

    Index = (int)  floor((log10(Array[j][0]) - logk_min)/logk_interval);

    if(Index < 0)  Index = 0;
    
    mean_modBin[Index]  += Array[j][0];
    modes_perbin[Index] += 1;

    Li                   = gsl_sf_legendre_P2(Array[j][1]);  // L_2 = 0.5*(3.*mu**2 -1.)
    Pi                   = Array[j][2];

    Sum_Li[Index]       += Li;
    Sum_Li2[Index]      += Li*Li;

    Sum_Pi[Index]       += Pi;
    Sum_PiLi[Index]     += Pi*Li;
  }

  int sum_modes = 0;
  
  for(j=0; j<modBinNumb; j++){
    mean_modBin[j] /= modes_perbin[j];

    sum_modes      += modes_perbin[j];
    
    // For a matrix A, paramater vector, (P_0, P_2)^T, P and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.

    // A reads as (a b) for a = sum_Modes 1, b = sum_Modes L_i, c = sum_Modes L_i, d = sum_Modes L_i**2.
    //            (c d)

    // and B = (b1, b2)^T for b1 = sum_Modes measured Pi, b2 = sum_Modes (measured Pi)*Li

    detA[j]        = modes_perbin[j]*Sum_Li2[j] - Sum_Li[j]*Sum_Li[j];

    Monopole[j]    = (1./detA[j])*( Sum_Li2[j]*Sum_Pi[j] - Sum_Li[j]*Sum_PiLi[j]);
    Quadrupole[j]  = (1./detA[j])*( -Sum_Li[j]*Sum_Pi[j] + modes_perbin[j]*Sum_PiLi[j]);

    // if(log10(detA[j]) > -6.0)
    printf("\n%le \t %le \t %le \t %d", mean_modBin[j], Monopole[j], Quadrupole[j], modes_perbin[j]);
  }

  printf("\n%d \t %d", sum_modes, n0*n1*n2);
  
  return 0;
}
