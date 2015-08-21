/* Integrals of the form 

    int_0^1 mu^n * (1. + 0.5*pow(k*mu*sigma_pair, 2.))^-1 dmu

    where n is the "order".
*/

double muOrderZero(double ks){
    if(ks > 0.04)  return pow(2., 0.5)*atan(ks*pow(2., -0.5))/ks;

    else{
      return 1.0 - 0.166667*pow(ks, 2.) + 0.05*pow(ks, 4.) - 0.0178571*pow(ks, 6.);
    }
}


double muOrderTwo(double ks){
      if(ks > 0.051675)  return pow(ks, -3.)*(2.*ks - 2.*sqrt(2.)*atan(ks*pow(2., -0.5))); 

      else{
        return 1./3. - 0.1*pow(ks, 2.) + 0.0357143*pow(ks, 4.) - 0.0138889*pow(ks, 6.);
      }
}


double muOrderFour(double ks){
      if(ks > 0.0512935)  return pow(ks, -5.)*(-4.*ks +(2./3.)*pow(ks, 3.) + sqrt(32.)*atan(ks/sqrt(2.))); 

      else{
        return 0.2 - 0.0714286*pow(ks, 2.) + 0.0277778*pow(ks, 4.) - 0.0113636*pow(ks, 6.);
      }
}


double muOrderSix(double ks){
      if(ks > 0.039329)  return pow(ks, -7.)*(8.*ks - (4./3.)*pow(ks, 3.) + 0.4*pow(ks, 5.) -8.*sqrt(2.)*atan(ks/sqrt(2.)));

      else{
         return 0.142857 - 0.0555556*pow(ks, 2.) + 0.0227273*pow(ks, 4.) - 0.00961538*pow(ks, 6.);
      }
}


double muOrderEight(double ks){
    if(ks>0.181839)  return pow(ks, -9.)*(-16.*ks + (8./3.)*pow(ks, 3.) - 0.8*pow(ks, 5.) + (2./7.)*pow(ks, 7.) + 16.*sqrt(2.)*atan(ks/sqrt(2.)));

    else{
      return 1./9. - 0.0454545*ks*ks + 0.0192308*pow(ks, 4.) - 0.00833333*pow(ks, 6.);
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
    return (9./8.)*(35.*beta*beta*muOrderEight(ks)  -30.*beta*beta*muOrderSix(ks) + 3.*beta*beta*muOrderFour(ks) + 70.*beta*muOrderSix(ks) - 60.*beta*muOrderFour(ks)  + 6.*beta*muOrderTwo(ks) + 35.*muOrderFour(ks) - 30.*muOrderTwo(ks) + 3.*muOrderZero(ks));
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
    }
}


int setLorentzianRSD(){
    pt2RSD_k  = &kaiserLorentz_multipole;
    
    sprintf(theoryRSD_flag, "LorentzianRSD");

    return 0;
}


// Print model to file. 

int DispersionModel_Multipoles(){
    double waveNumber    = 0.0;
    
    sprintf(filepath, "%s/Data/Multipoles/KaiserLorentzMultipoles_Pk_%s_beta_%.2f_velDispersion_%.2f.dat", root_dir, theoryPk_flag, beta, velDispersion);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<100000; j++){
        waveNumber = pow(10., -5.)*j;

        fprintf(output, "%le \t %le \t %le \t %le \t %le \n", waveNumber, (*pt2Pk)(waveNumber), (*pt2Pk)(waveNumber)*kaiserLorentz_Monofactor(waveNumber*velDispersion, beta), (*pt2Pk)(waveNumber)*kaiserLorentz_Quadfactor(waveNumber*velDispersion, beta), (*pt2Pk)(waveNumber)*kaiserLorentz_Hexfactor(waveNumber*velDispersion, beta));
        
    }
    
    fclose(output);
    
    /*
    sprintf(filepath, "%s/Data/Multipoles/muOrder_nLorentzIntegrals.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<1000000; j++){
        waveNumber = pow(10., -5.)*j;
        
        fprintf(output, "%.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf\n", waveNumber, (*pt2Pk)(waveNumber), muOrderZero(velDispersion*waveNumber), muOrderTwo(velDispersion*waveNumber), muOrderFour(velDispersion*waveNumber), muOrderSix(velDispersion*waveNumber), muOrderEight(velDispersion*waveNumber));
    }
    
    fclose(output);
    */
    
    return 0;
}
