double sig2_integrand(double lnk){
    // integrand for the dlnk integral used to calculate sigma^2(R).
    double     k;
    double a_eff;
    
    // halofit units for k. 
    k     = exp(lnk)*2998.;
    // for(k=1.0e2; k<=1.0e7; k*=1.5)  printf("\n%e \t %e \t %e", k/2998.0, k*k*k/2./pi/pi*P_L(0.5, k), k*k*k/2./pi/pi*P_NL(0.5, k));

    a_eff = 1./(1.+z_eff); 

    // halofit units for k. 
    return k*k*k/2./pi/pi*P_L(a_eff, k)*pow(spherical_tophat(k/2998., sig2_R), 2.);
}


double splint_sig2(double R){
    // splint sigma^2(R)
    double Interim;

    splint(radii, sig2, sig2_2D, sigres, R, &Interim);
    
    return Interim;
}


double splint_msig2(double M){
    // splint sigma^2(M)
    double       R;
    
    R = pow(3.*M/(4.*pi*rhobox), 1./3.);
    
    return  splint_sig2(R);
}


double fs_tinker(double sigma, double z){
    // returns Tinker 2008 'universal' f(sigma).
    // possible divergent wrt integrals over total mass. see appendix B of tinker 2010.

    double Delta = 200.;
    double A     = 0.186*pow(1.+ z_eff, -0.14);
    double a     = 1.470*pow(1.+ z_eff, -0.06);
    
    // redshift evolution of b.
    double alpha = exp(-pow(0.75/log(200./75.), 1.2));
    double b     = 2.570*pow(1.+ z_eff, -alpha);
    double c     = 1.190;
    
    double fs;
    
    fs           = A*(1. + pow(sigma/b, -a))*exp(-c*pow(sigma, -2.)); 

    return fs;
}


double fs_Jenkins(double sigma){
    return 0.315*exp(-fabs(pow(log(1./sigma) + 0.61, 3.8)));
}


double nm_tinker(double M, double z){
    // Tinker 2008 'universal' mass function.

    double R;
    double error;
    double d_dMln_invSig;
    
    double sigma;
    
    R             = pow(3.*M/(4.*pi*rhobox*pow(1.+z, 3.)), 1./3.);

    sigma         = sqrt(splint_msig2(M));

    d_dMln_invSig = pow(-6.*splint_msig2(M), -1.)*pow(3./(4.*pi*rhobox*pow(1.+z, 3.)*M*M), 1./3.)*dfridr(&splint_sig2, R, 0.5, &error);

    // redshift evolution of rhobox is (1+z)**3, 
    return fs_tinker(sigma, z)*(rhobox/M)*pow(1.+z, 3.)*d_dMln_invSig;
}


double splint_ln_invSig(double lnM){
    // splint sigma^2(R)
    double Interim;

    splint(lnm_jenkins, ln_invSig_jenkins, ln_invSig_jenkins_2D, sigres, lnM, &Interim);
    
    return Interim;
}


double press_Schechter(double M){
    // assumes an EdS cosmology, \delta_sc = 1.686 at (z=0);
    
    double nu;
    double sig;
    double ln_invsig;
    double delta_sc = 1.686;
    
    double nu_fnu;
    
    ln_invsig = splint_ln_invSig(log(M));

    // sigma(M) 
    sig       = 1./exp(ln_invsig);
    
    nu        = pow(delta_sc/sig, 2.);
    
    nu_fnu    = sqrt(nu/(2.*pi))*exp(-0.5*nu); 
    
    // deal with dnu
    return 0.;
}


int mass_fn(){
    double     k;
    double  mass;
    double error;
    double   sig;
    
    double  jenkins_dndm;
    double  tinker_dndm;
    
    sigres  = 100;
    
    sig2    = malloc(sigres*sizeof(*sig2));
    radii   = malloc(sigres*sizeof(*radii));
    sig2_2D = malloc(sigres*sizeof(*sig2_2D));
    
    lnm_jenkins             = malloc(sigres*sizeof(double));
    ln_invSig_jenkins       = malloc(sigres*sizeof(double));
    ln_invSig_jenkins_2D    = malloc(sigres*sizeof(double));
    
    for(j=0; j<sigres; j++){
       // R global variable for passing under the hood. 
       sig2_R   = 1. + 1.*j;
    
       radii[j] = sig2_R;

       // sigma^2(R) from integral dln k. 
       sig2[j]  = qromb(sig2_integrand, -7., 8.);
    
       // redshift zero.
       lnm_jenkins[j]         = log((4.*pi/3.)*pow(radii[j], 3.)*rhobox);
    
       ln_invSig_jenkins[j]   = log(pow(sig2[j], -0.5)); 
    }
    
    // for calculation of sigma^2(R)
    spline(radii, sig2, sigres, 1.0e31, 1.0e31, sig2_2D);
    
    // Jenkins mass function. 
    spline(lnm_jenkins, ln_invSig_jenkins, sigres, 1.0e31, 1.0e31, ln_invSig_jenkins_2D);
    
    // Jenkins & Tinker mass fn's.
    sprintf(filepath, "%s/Data/500s/Jenkins_nm.dat", root_dir);
    
    output = fopen(filepath, "w");

    for(j=0; j<110; j++){  
        mass = pow(10., 11.4 + 4.*0.01*j);
        
        sig  = sqrt(splint_msig2(mass));
    
                                      // d/d_lnM ln(sig^-1).
        jenkins_dndm = fs_Jenkins(sig)*rhobox*dfridr(&splint_ln_invSig, log(mass), 0.5, &error);
        
        tinker_dndm  = fs_tinker(sig, 0.0)*rhobox*dfridr(&splint_ln_invSig, log(mass), 0.5, &error);
    
        fprintf(output, "%e \t %e \t %e \n", log10(mass), log10(jenkins_dndm*log(10.)/mass), log10(tinker_dndm*log(10.)/mass));
    }
    
    fclose(output);
    
    return 0;
}


double b_Kaiser(double ln_nu){
    // Cole and Kaiser halo bias. 
    double delta_sc = 1.686; // assumes EdS, z=0?
    double nu;
    
    double sig;
    double ln_invsig;
    
    // ln_invsig = splint_ln_invSig(log(M));

    // sigma(M) 
    // sig       = 1./exp(ln_invsig);    
        
    // watch out for nu, nu sq.! (pg 21. Cooray, Tinker bias).
    // nu        = delta_sc/sig;
    
    nu        = pow(10., ln_nu);
    
    return 1. + (nu*nu - 1.)/delta_sc;
}


double b_tinker(double ln_nu){
    // Cole and Kaiser halo bias. 
    double delta_sc = 1.686; // assumes EdS, z=0?
    double nu;

    double Delta = 200.;
    double y     = log10(Delta);

    double A     = 1.0 + 0.24*y*exp(-1.*pow(4./y, 4.));
    double a     = 0.44*y - 0.88;
    double B     = 0.183;
    double b     = 1.5;
    double C     = 0.019 + 0.107*y + 0.19*exp(-1.*pow(4./y, 4.));
    double c     = 2.4;
    
    double delta_ac = 1.0; // ??
    
    double b_nu;
    
    nu       = pow(10., ln_nu);
    
    // ASSUMES MASS FUNCTIONS DEFINED IN APPENDIX C of tinker et al., change. 
    b_nu     = 1. - A*pow(nu, a)/(pow(nu, a) + delta_ac) + B*pow(nu, b) + C*pow(nu, c);
    
    return b_nu;
}


int halobias_fn(){
    double ln_nu;
    
    /*
    sprintf(filepath, "%s/Data/500s/halobias_pk.dat", root_dir);
    
    output = fopen(filepath, "w");

    for(j=1; j<100000; j++){  
        kval = 0.00001*j;
            
        fprintf(output, "%e \t %e \t %e\n", kval, pow(2998., 3.)*P_L(1.0, kval*2998.), pow(2998., 3.)*P_NL(0.9999, kval*2998.));
    }
    
    fclose(output); */
    
    
    sprintf(filepath, "%s/Data/500s/halobias.dat", root_dir);
    
    output = fopen(filepath, "w");

    for(j=1; j<100; j++){  
        // base 10. 
        ln_nu = -0.4 + 0.01*j*1.;
            
        fprintf(output, "%e \t %e \t %e \n", ln_nu, b_Kaiser(ln_nu), b_tinker(ln_nu));
    }
    
    fclose(output);
    
    return 0;
}
