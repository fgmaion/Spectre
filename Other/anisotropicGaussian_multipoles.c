double anisotropicGaussian_monopoleIntegrand(double mu){
    return 0.5*exp(-0.5*pow(anisoGauss_delta*mu, 2.)*(0.5/pow(anisoGauss_bsigma, 2.) - 0.5/pow(anisoGauss_asigma, 2.)));
}


double anisotropicGaussian_quadrupoleIntegrand(double mu){
    // (2l + 1)
    return 5.*0.5*LegendrePolynomials(mu, 2)*exp(-0.5*pow(anisoGauss_delta*mu, 2.)*(0.5/pow(anisoGauss_bsigma, 2.) - 0.5/pow(anisoGauss_asigma, 2.)));
}


double anisotropicGaussian_hexadecapoleIntegrand(double mu){
    // (2l + 1)
    return 9.*0.5*LegendrePolynomials(mu, 4)*exp(-0.5*pow(anisoGauss_delta*mu, 2.)*(0.5/pow(anisoGauss_bsigma, 2.) - 0.5/pow(anisoGauss_asigma, 2.)));
}


double anisotropicGaussian_octupoleIntegrand(double mu){
    // (2l + 1)
    return 13.*0.5*LegendrePolynomials(mu, 6)*exp(-0.5*pow(anisoGauss_delta*mu, 2.)*(0.5/pow(anisoGauss_bsigma, 2.) - 0.5/pow(anisoGauss_asigma, 2.)));
}


int anisotropicGaussian_multipoles(){
    anisoGauss_asigma = 60.;
    anisoGauss_bsigma = 90.;
    
    double monopole, quadrupole, hexadecapole, octupole;

    sprintf(filepath, "%s/Data/stacpolly/anisotropicGaussian_monopolePrediction.dat", root_dir);
    
    output = fopen(filepath, "w");

    for(j=0; j<1000; j++){
        anisoGauss_delta  = j*1.;
    
        monopole     = exp(-0.25*pow(anisoGauss_delta/anisoGauss_asigma, 2.))*qromb(anisotropicGaussian_monopoleIntegrand,     -1.0, 1.0); 
        
        quadrupole   = exp(-0.25*pow(anisoGauss_delta/anisoGauss_asigma, 2.))*qromb(anisotropicGaussian_quadrupoleIntegrand,   -1.0, 1.0); 

        hexadecapole = exp(-0.25*pow(anisoGauss_delta/anisoGauss_asigma, 2.))*qromb(anisotropicGaussian_hexadecapoleIntegrand, -1.0, 1.0); 
        
        octupole     = exp(-0.25*pow(anisoGauss_delta/anisoGauss_asigma, 2.))*qromb(anisotropicGaussian_octupoleIntegrand,     -1.0, 1.0); 

        fprintf(output, "%e \t %e \t %e \t %e \t %e \n", anisoGauss_delta, monopole, quadrupole, hexadecapole, octupole);
    }

    fclose(output);

    return 0;
}
