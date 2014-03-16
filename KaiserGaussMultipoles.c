/* Integrals of the form 

    int_0^1 mu^n * exp(-k*k*mu*mu*sigma*sigma) dmu

    where n is the "order".

*/

double muOrderZero(double ks){
    // limits were established with signa = 2.*3./sqrt(2.), limits should scale propto velDispersion/(2.*3./sqrt(2.))
    if(ks > 0.0914289)  return 0.5*pow(pi, 0.5)*gsl_sf_erf(ks)/ks;

    else{
        return 1.0 - 0.333333*pow(ks, 2.) + 0.1*pow(ks, 4.) - 0.0238095*pow(ks, 6.);
    }
}


double muOrderTwo(double ks){
    if(ks > 0.15727469)  return pow(4.*ks*ks*ks, -1.)*(pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-1.*ks*ks));
    
    else{
        return 1./3. - ks*ks/5. + pow(ks, 4.)/14.;
    }
}


double muOrderFour(double ks){
    if(ks>0.2418305)  return pow(8.*pow(ks, 5.), -1.)*(3.*pow(pi, 0.5)*gsl_sf_erf(ks) -6.*ks*exp(-1.*ks*ks) - 4.*pow(ks, 3.)*exp(-ks*ks));

    else{
        return 1./5. - ks*ks/7. + pow(ks, 4.)/18.;
    }
}


double muOrderSix(double ks){
    if(ks>0.335168)  return pow(16.*pow(ks, 7.), -1.)*(15.*pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-ks*ks)*(15. + 10.*ks*ks + 4.*pow(ks, 4.)));

    else{
         return 1./7. - ks*ks/9. + pow(ks, 4.)/22.;
    }
}


double muOrderEight(double ks){
    if(ks>0.4326645)  return (105.*pow(pi, 0.5)*gsl_sf_erf(ks)*pow(32.*pow(ks, 9.), -1.) - 2.*ks*exp(-ks*ks)*pow(32.*pow(ks, 9.), -1.)*(8.*pow(ks, 6.) + 28.*pow(ks, 4.) + 70.*pow(ks, 2.) + 105.));

    else{ 
        // Numerical result diverges for ks < 0.2, replace by Taylor expansion at ksigma=0.5
         return 1./9. - ks*ks/11. + pow(ks, 4.)/26.;
    }
}



// Multipoles for the Kaiser-Gauss redshift space distortion model.

double kaiserGauss_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}


double kaiserGauss_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}


double kaiserGauss_Hexfactor(double ks, double beta){
    return (9./8.)*(35.*beta*beta*muOrderEight(ks)  -30.*beta*beta*muOrderSix(ks) + 3.*beta*beta*muOrderFour(ks) + 70.*beta*muOrderSix(ks) - 60.*beta*muOrderFour(ks)  + 6.*beta*muOrderTwo(ks) + 35.*muOrderFour(ks) - 30.*muOrderTwo(ks) + 3.*muOrderZero(ks));
}

// Print model to file. 

int kaiser_nonlinearSuppression_Multipoles(){
    double waveNumber    = 0.0;
    
    sprintf(filepath, "%s/Data/Multipoles/KaiserGaussMultipoles_Beta_%.2f_velDispersion_%.2f.dat", root_dir, beta, velDispersion);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<100000; j++){
        waveNumber = pow(10., -5.)*j;

        fprintf(output, "%le \t %le \t %le \t %le \t %le \n", waveNumber, (*pt2Pk)(waveNumber), (*pt2Pk)(waveNumber)*kaiserGauss_Monofactor(waveNumber*velDispersion, beta), (*pt2Pk)(waveNumber)*kaiserGauss_Quadfactor(waveNumber*velDispersion, beta), (*pt2Pk)(waveNumber)*kaiserGauss_Hexfactor(waveNumber*velDispersion, beta));
        
    }
    
    fclose(output);
    
    sprintf(filepath, "%s/Data/Multipoles/muOrder_nGaussIntegrals.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<1000000; j++){
        waveNumber = pow(10., -5.)*j;
        
        fprintf(output, "%.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf \t %.12lf\n", waveNumber, (*pt2Pk)(waveNumber), muOrderZero(velDispersion*waveNumber), muOrderTwo(velDispersion*waveNumber), muOrderFour(velDispersion*waveNumber), muOrderSix(velDispersion*waveNumber), muOrderEight(velDispersion*waveNumber));
    }
    
    fclose(output);
    
    return 0;
}
