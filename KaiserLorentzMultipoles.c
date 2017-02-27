/* Integrals of the form 

    int_0^1 mu^n * (1. + 0.5*pow(k*mu*sigma_pair, 2.))^-1 dmu, where n is the "order".
*/

double muOrderZero(double ks){
    if(ks > 0.04)  return pow(2., 0.5)*atan(ks*pow(2., -0.5))/ks;

    else{
      return 1.0 - 0.166667*pow(ks, 2.) + 0.05*pow(ks, 4.) - 0.0178571*pow(ks, 6.);
    }
}

double muOrderTwo(double ks){
      if(ks > 0.051675)  return 2./pow(ks,2.) - 2.*sqrt(2.)*atan(ks*pow(2., -0.5))/pow(ks, 3.); 

      else{
        return 1./3. - 0.1*pow(ks, 2.) + 0.0357143*pow(ks, 4.) - 0.0138889*pow(ks, 6.);
      }
}

double muOrderFour(double ks){
      if(ks > 0.0512935)  return -4./pow(ks, 4.) +(2./3.)/pow(ks, 2.) + 4.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 5.); 

      else{
        return 1./5. - 0.0714286*pow(ks, 2.) + 0.0277778*pow(ks, 4.) - 0.0113636*pow(ks, 6.);
      }
}

double muOrderSix(double ks){
      if(ks > 0.039329)  return 8./pow(ks,6.) - (4./3.)*pow(ks, -4.) + 0.4*pow(ks, -2.) -8.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 7.);

      else{
         return 1./7. - 0.0555556*pow(ks, 2.) + 0.0227273*pow(ks, 4.) - 0.00961538*pow(ks, 6.);
      }
}

double muOrderEight(double ks){
    if(ks>0.181839)  return -16.*pow(ks,-8.) + (8./3.)*pow(ks, -6.) - 0.8*pow(ks, -4.) + (2./7.)*pow(ks, -2.) + 16.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 9.);

    else{
      return 1./9. - 0.0454545*ks*ks + 0.0192308*pow(ks, 4.) - 0.00833333*pow(ks, 6.);
    }
}

double muOrderTen(double ks){
    if(ks>0.14)  return (2./315.)*(35./pow(ks,2.) -90./pow(ks,4.) +252./pow(ks,6.) -840./pow(ks,8.) +5040./pow(ks,10.) -5040.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks,11.));

    else{
        return 1./11. - ks*ks/26. + pow(ks, 4.)/60. - pow(ks, 6.)/136. + pow(ks, 8.)/304.;
    }
}

// Multipoles for the Kaiser-Lorentz redshift space distortion model.
double kaiserLorentz_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}

double kaiserLorentz_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}

double kaiserLorentz_Hexfactor(double ks, double beta){
    return (9./8.)*(35.*beta*beta*muOrderEight(ks)  + 10.*beta*(7. -3.*beta)*muOrderSix(ks) + (35. - 60.*beta + 3.*beta*beta)*muOrderFour(ks) + 6.*(beta - 5.)*muOrderTwo(ks) + 3.*muOrderZero(ks));
}

double kaiserLorentz_Octfactor(double ks, double beta){
    return (13./16.)*(231.*beta*beta*muOrderTen(ks) + (462.*beta - 315.*beta*beta)*muOrderEight(ks)  + (231.-630.*beta +105.*beta*beta)*muOrderSix(ks) + (-315.+210.*beta-5.*beta*beta)*muOrderFour(ks) + (105.-10.*beta)*muOrderTwo(ks) - 5.*muOrderZero(ks));
}


double kaiserLorentz_multipole(double ks, double beta, int monoQuad){
    // Mono, L_0 corresponds to 0. Quad, L_2 corresponds to 2.  Hex, L_4 corresponds to 4.
    switch(monoQuad){
        case 0:
            return kaiserLorentz_Monofactor(ks, beta);
        case 2:
            return kaiserLorentz_Quadfactor(ks, beta);
        case 4:
            return kaiserLorentz_Hexfactor(ks, beta);
        case 6:
            return kaiserLorentz_Octfactor(ks, beta);
    }
}


int prep_kaiserLorentSpline(){
  sp_kL_N = 1000;
  
  sp_kL_ks       = malloc(sp_kL_N*sizeof(double));    
  sp_kL_mu0      = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu2      = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu4      = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu6      = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu8      = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu10     = malloc(sp_kL_N*sizeof(double));

  sp_kL_mu0_2d   = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu2_2d   = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu4_2d   = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu6_2d   = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu8_2d   = malloc(sp_kL_N*sizeof(double));
  sp_kL_mu10_2d  = malloc(sp_kL_N*sizeof(double));

  for(j=0; j<sp_kL_N; j++){
     sp_kL_ks[j]  = (double) j*(20./sp_kL_N); // sp_kL_N bins between k*s = 0. and 20.;     

    sp_kL_mu0[j]  = muOrderZero(sp_kL_ks[j]);
    sp_kL_mu2[j]  = muOrderTwo(sp_kL_ks[j]);
    sp_kL_mu4[j]  = muOrderFour(sp_kL_ks[j]);
    sp_kL_mu6[j]  = muOrderSix(sp_kL_ks[j]);
    sp_kL_mu8[j]  = muOrderEight(sp_kL_ks[j]);
    sp_kL_mu10[j] = muOrderTen(sp_kL_ks[j]);
  }

  spline(sp_kL_ks,  sp_kL_mu0,  sp_kL_N, 1.0e31, 1.0e31,  sp_kL_mu0_2d);
  spline(sp_kL_ks,  sp_kL_mu2,  sp_kL_N, 1.0e31, 1.0e31,  sp_kL_mu2_2d);
  spline(sp_kL_ks,  sp_kL_mu4,  sp_kL_N, 1.0e31, 1.0e31,  sp_kL_mu4_2d);
  spline(sp_kL_ks,  sp_kL_mu6,  sp_kL_N, 1.0e31, 1.0e31,  sp_kL_mu6_2d);
  spline(sp_kL_ks,  sp_kL_mu8,  sp_kL_N, 1.0e31, 1.0e31,  sp_kL_mu8_2d);
  spline(sp_kL_ks, sp_kL_mu10,  sp_kL_N, 1.0e31, 1.0e31, sp_kL_mu10_2d);
  
  return 0;
}


double seqsp_kLmu(double ks, int order, int* klo){
  if(ks > 19.9)  return 0.0; // Make smoother? 
  
  switch(order){
   case 0:
    return sequential_splint(sp_kL_ks, sp_kL_mu0,   sp_kL_mu0_2d, sp_kL_N, ks, klo);

   case 2:
    return sequential_splint(sp_kL_ks, sp_kL_mu2,   sp_kL_mu2_2d, sp_kL_N, ks, klo);

   case 4:
    return sequential_splint(sp_kL_ks, sp_kL_mu4,   sp_kL_mu4_2d, sp_kL_N, ks, klo);

  case 6:
    return sequential_splint(sp_kL_ks, sp_kL_mu6,   sp_kL_mu6_2d, sp_kL_N, ks, klo);

  case 8:
    return sequential_splint(sp_kL_ks, sp_kL_mu8,   sp_kL_mu8_2d, sp_kL_N, ks, klo);

  case 10:
    return sequential_splint(sp_kL_ks, sp_kL_mu10, sp_kL_mu10_2d, sp_kL_N, ks, klo);
  }
}

 
double noRSD(double ks, double beta, int monoQuad){
    // Mono, L_0 corresponds to 0. Quad, L_2 corresponds to 2.  Hex, L_4 corresponds to 4.
    switch(monoQuad){
        case 0:
            return 1.0;
        case 2:
            return 0.0;
        case 4:
            return 0.0;
        case 6:
            return 0.0;
        case 8:
            return 0.0;
    }
}

int setLorentzianRSD(){
    pt2RSD_k  = &kaiserLorentz_multipole;
    
    sprintf(theoryRSD_flag, "LorentzianRSD");

    return 0;
}

// Print model to file. 
int DispersionModel_Multipoles(){
    double k;
    double ks            = 0.0;
    double beta          = 0.5;
    double velDispersion = 2.0;
    
    // sprintf(filepath, "%s/Data/Multipoles/KaiserLorentzMultipoles_Pk_%s_beta_%.2f_velDispersion_%.2f.dat", root_dir, theoryPk_flag, beta, velDispersion);
    // sprintf(filepath, "%s/Data/maskedRSD_draftwork/KaiserLorentzian_multipoleConvergence_beta_%.2lf_sigma_%.2lf.dat", root_dir, beta, velDispersion);
    sprintf(filepath, "%s/W1_Spectro_V7_1/KaiserLorentzian_multipole_beta_%.2lf_sigma_%.2lf.dat", root_dir, beta, velDispersion);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<100000; j++){
        k   = pow(10., -3.)*j;
    
        ks  = k*velDispersion;

        fprintf(output,"%le \t %le \t %le \t %le \t %le \n", k, (*pt2Pk)(k)*kaiserLorentz_multipole(ks,beta,0), (*pt2Pk)(k)*kaiserLorentz_multipole(ks,beta,2), (*pt2Pk)(k)*kaiserLorentz_multipole(ks,beta,4), (*pt2Pk)(k)*kaiserLorentz_multipole(ks,beta,6));
        
        // fprintf(output,"%le \t %le \n", ks, muOrderTwo(ks));
    }
    
    fclose(output);
        
    return 0;
}
