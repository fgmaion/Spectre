// Gaussian model
// eqn. (2.4) of arXiv 9412062.
float A(float kappa){
    if(kappa < 0.01)    return 1. - kappa*kappa/3.;

    else                return sqrt(pi)*gsl_sf_erf(kappa)*pow(2.*kappa, -1.);
}

float B(float kappa){
    if(kappa < 0.01)    return 1. - 3.*kappa*kappa/5.;

    else                return 3.*pow(2.*kappa*kappa, -1.)*(A(kappa) - exp(-kappa*kappa));
}

float C(float kappa){
    if(kappa < 0.01)    return 1. - 5.*kappa*kappa/7.;

    else                return 5.*pow(2.*kappa*kappa, -1.)*(B(kappa) - exp(-kappa*kappa));
}

float D(float kappa){
    // Need taylor expansion. 
    return B(kappa)*pow(kappa, -2.);
}

float E(float kappa){
    // Need taylor expansion. 
    return C(kappa)*pow(kappa, -2.);
}


float i_n(int n){
    switch(n){
        case 0:
            return  1.0;
        case 2:
            return -1.0;
        case 4: 
            return  1.0;
    }
}

// Correlation fn. monopole contribution from A, B, C.

float transformIntegrand(float x){
    return i_n(besseltransform_order)*pow(2.*pi*pi, -1.)*x*x*(*pt2_pkIntegrand)(x)*gsl_sf_bessel_jl(besseltransform_order, x*q0)*(*pt2Pk)(x/velDispersion);   
}


float inverse_transformIntegrand(float x){
    // In this instance (*pt2_pkIntegrand)(x) will point to the correlation function of a given multipole order. 
    return i_n(besseltransform_order)*4.*pi*x*x*(*pt2_pkIntegrand)(x)*gsl_sf_bessel_jl(besseltransform_order, x*q0); 
}


double Integrate(double a, double b, double dx){
    int    bins;
    double Interim = 0.0;
    
    bins = (int) ceil((b-a)/dx);

    for(j=0; j<bins; j++){  
        Interim += dx*(*transformIntegrand)((j+0.5)*dx);
    }

    return Interim;
}


int Mono_xi(){    
    double q0Array[100], mono_xi[100], quad_xi[100], clipped_xi[100], clipped_xi2D[100];
    
    double kmax   =  50., Integral, xInterval;
    
    velDispersion = 0.5;
    
    xInterval     = 0.05;
    
    sprintf(filepath, "%s/Data/SpectralDistortion/Gaussian_zSpaceCorrfn_sigma_%.2f_beta_%.2f.txt", root_dir, velDispersion, beta);
    
    output = fopen(filepath, "w");
    
    for(ii=0; ii<100; ii++){
        q0                    = pow(10., 0. + ii*0.02)/velDispersion;
        
        q0Array[ii]           = q0;
                
        besseltransform_order = 0;
        
        mono_xi[ii]           = 0.0;
        
        pt2_pkIntegrand       = &A;
        
        // Integral           = qromb(transformIntegrand, 0.0, kmax);
    
        Integral              = Integrate(0.0, kmax, xInterval);
    
        mono_xi[ii]          += Integral;
        
        pt2_pkIntegrand       = &B;
    
        Integral              = Integrate(0.0, kmax, xInterval);
    
        mono_xi[ii]          += (2./3.)*beta*Integral;
        
        pt2_pkIntegrand       = &C;

        Integral              = Integrate(0.0, kmax, xInterval);
    
        mono_xi[ii]          += (1./5.)*beta*beta*Integral;
        
        
        // Quadrupole
        besseltransform_order = 2;
        
        quad_xi[ii]           = 0.0; 
    
        pt2_pkIntegrand       = &A;
    
        Integral              = Integrate(0.0, kmax, xInterval);
        
        quad_xi[ii]           += -5./2.*Integral;
        
        pt2_pkIntegrand       = &B;

        Integral              = Integrate(0.0, kmax, xInterval);
    
        quad_xi[ii]          += ((5./2.) + (4./3. - 3.)*beta)*Integral;
    
        pt2_pkIntegrand       = &C;
    
        Integral              = Integrate(0.0, kmax, xInterval);
        
        quad_xi[ii]          += (3.*beta + beta*beta)*Integral;
        
        pt2_pkIntegrand       = &D;
        
        Integral              = Integrate(0.0, kmax, xInterval);
        
        quad_xi[ii]          += (-15./4.)*beta*beta*Integral;
        
        pt2_pkIntegrand       = &E;
        
        Integral              = Integrate(0.0, kmax, xInterval);
        
        quad_xi[ii]          += (15./4.)*beta*beta*Integral;
        
        mono_xi[ii]          *= pow(velDispersion, -3.);
        quad_xi[ii]          *= pow(velDispersion, -3.);
    }

    // Quadrupole.
    
    double var, u0;
    
    var = mono_xi[1];
    
    u0       = inverse_erf(2.*sqrt(0.5) -1.);
    
    for(jj=0; jj<100; jj++){
        q0               = 1. + jj*0.01;
    
        clipped_xi[jj]   = 0.25*pow(1.0 + gsl_sf_erf(u0), 2.)*quad_xi[jj];
        
        // clipped_xi[jj]  += pow(var, -1.)*(2.*mono_xi[jj]*quad_xi[jj] + 0.285714*quad_xi[jj]*quad_xi[jj]);
    
        fprintf(output, "%f \t %f \t %f \t %f\n", q0Array[jj]*velDispersion, mono_xi[jj], quad_xi[jj], clipped_xi[jj]);
    }

    fclose(output);
    
    spline(q0Array, clipped_xi, 99, 1.0e31, 1.0e31, clipped_xi2D);

    return 0;
}
