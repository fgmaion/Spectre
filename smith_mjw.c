#include <malloc.h>
#include <assert.h>


int set_testparams(){
    Om_m    =      0.25;
    Om_v    = 1. - Om_m;
    Om_b    =     0.045;
    h       =      0.73;
    sigma_8 =       0.9; 
    p_index =       1.0;
         
    return 0;
}


int halofit(){
    // The `halofit' code models the nonlinear evolution of cold matter 
    // cosmological power spectra. The full details of the way in which 
    // this is done are presented in Smith et al. (2003), MNRAS, 341, 1311. 

    // The code `halofit' was written by R. E. Smith & J. A. Peacock. 

    // Updated to use full E+H CDM 17.11.2004

    // Further edits 17.11.2008 by JAP

    // set_testparams();

    double z0 = 0.85; // 0.75
    
    double     aexp; // expansion factor for desired redshift. 
    
    double    Om_mz;
    double    Om_vz;
    
    // usual Om_m problem. Calculation of the shape parameter.
    pk_gamma = Om_m*h*exp(-Om_b*(1. + sqrt(2.*h))/Om_m);
    
    double sig8_app, gg, gg0, plin, pnlin; 
    double rknl, rneff, rncur, dummy, rk;
    
    printf("\n\nHalofit.");
    
    printf("\n\nOm_m: %.3e, Om_b: %.3e, Om_v: %.3e, h: %.3e, sigma 8: %.3e, p index: %.3e", Om_m, Om_b, Om_v, h, sigma_8, p_index);
    
    printf("\n\nDesired redshift: %.2f, scale factor: %.2e, Gamma: %.8e, Om_m(z): %.2e, Om_v(z): %.2e, CPT D+: %.2e\n", z0, aexp, pk_gamma, Om_mz, Om_vz, gg);
    
    for(ii=0; ii<1; ii++){
        // possibility of multiple redshifts. 
    
        aexp  = 1./(1. + z0);          
    
        Om_mz = omega_m(aexp, Om_m, Om_v);
        Om_vz = omega_v(aexp, Om_m, Om_v);
    
        // checked against jap.
        // printf("\nz: %d, scale factor: %.6e, Om_m(z): %.6e, Om_v(z): %.6e", ii, aexp, Om_mz, Om_vz);
        
        // checked against jap. 
        gg0 = CarollPressTurner(Om_m,  Om_v);    
        gg  = CarollPressTurner(Om_mz, Om_vz);
        
        // printf("\nz: %d, D_+(z): %.6e, D_+(0): %.6e", i, gg, gg0);
    
        // embedded in BunnWhite(). 
        amp = 1.;
    
        // for sigma_8 apparent calc. amp is always 1. sigma_8 calculated from Del^2 at redshift 0.
        sig8_app = sigint(8.0);
        
        // sigma_8 (defined today.), checked against jap.  
        // printf("\napparent sigma_8: %.16e", sig8_app);
    
        // bit baffled by aexp here. 
        amp = (aexp*gg/gg0)*(sigma_8/sig8_app);
    
        // printf("\namplitude correction: %.6e, apparent sigma_8: %.6e, current: %.6e, scaled to z=0: %.6e", amp, sig8_app, sigint(8.0), sigint(8.0)*(gg0/gg));
    
        rknlCalc(&rknl, &rneff, &rncur);
        
        
        sprintf(filepath, "%s/Data/500s/EisensteinHu_halofit_pk_%.2f.dat", root_dir, z0);
    
        output = fopen(filepath, "w");
        
        for(jj=0; jj<10000; jj++){
            rk    = -4.5 + 6.5*(jj - 1)/9900.;
        
            rk    = pow(10., rk);
            
            plin  = BunnWhite(rk);
        
        
            halofit_fittingfn(rk, rknl, rneff, rncur, Om_mz, Om_vz, plin, &pnlin);
        
        
            plin  /= pow(rk, 3.)/(2.*pi*pi);
    
            pnlin /= pow(rk, 3.)/(2.*pi*pi);
        
            fprintf(output, "%e \t %e \t %e \n", rk, plin, pnlin);
        }
    
        fclose(output);
    }

    return 0;
}


double BondEfstathiou(double rk, double gams, double p_index, double sig8){
    // Bond & Efstathiou (1984) approximation to the linear CDM power spectrum. 
    double p_gam, rkeff, q, q8, tk, tk8;

    rkeff = 0.172 + 0.011*log(gams/0.36)*log(gams/0.36);

    q     = pow(10., -20.) + rk/gams;

    q8    = pow(10., -20.) + rkeff/gams;

    tk    = 1./pow(1. + pow(6.4*q + pow(3.0*q, 1.5) + pow(1.7*q, 2.), 1.13), 1./1.13);

    tk8   = 1./pow(1. + pow(6.4*q8 + pow(3.0*q8, 1.5) + pow(1.7*q8, 2.), 1.13), 1./1.13);

    p_gam = sig8*sig8*pow(q/q8, 3. + p_index)*tk*tk/tk8/tk8;
      
    return p_gam;
}


double omega_m(double aa, double om_m0, double om_v0){
    // evolution of omega matter with expansion factor
    double omega_m, omega_t;
      
    omega_t = 1.0 + (om_m0 + om_v0 - 1.0)/(1. - om_m0 - om_v0 + om_v0*aa*aa + om_m0/aa);
      
    omega_m = omega_t*om_m0/(om_m0+om_v0*aa*aa*aa);
      
    return omega_m;
}


double omega_v(double aa, double om_m0, double om_v0){
    // evolution of omega lambda with expansion factor
    double omega_v, omega_t;
      
    omega_t=1.0+(om_m0+om_v0-1.0)/(1.-om_m0-om_v0+om_v0*aa*aa+om_m0/aa);
      
    omega_v=omega_t*om_v0/(om_v0+om_m0/aa/aa/aa);
      
    return omega_v;
}


double CarollPressTurner(double om_m, double om_v){
    // growth factor for linear fluctuations, given by Caroll, Press & Turner 1992, ARAA, 30, 499.
    double gg;

    gg = 2.5*om_m/(pow(om_m, 4./7.) - om_v + (1.0 + om_m/2.)*(1. + om_v/70.));
    
    return gg;
}


double sigint(double r){
    double t, y, rk, x, w, d2, sigint;
    
    int    nint = 10000;
    
    double sum1 = 0.0;
    
    // \sigma_8^2 = \int Del^2(k) W^2(k) d(lnk)
    // W(k) = 3/y^3 (sin(y) yada yada yada).
    
    // r radius of smoothing, e.g. 8;
    
    for(i=1; i<nint; i++){
        t     = (i - 0.5)/nint;
            
        // range of y, 2*(10^4) to 1.5*10^(-4)
        y     = -1.0 + 1.0/t;

        rk    = y;

        //  COBE normalised Del^2 at z=0, when amp = 1. otherwise scaled to appropriate redshift. 
        d2    = BunnWhite(rk); 

        x     = y*r;
        
        // smoothing kernel.
        w     = (3./x/x/x)*(sin(x) - x*cos(x));
    
        sum1 += w*w*d2/y/t/t;
    }
    
    sum1  /= nint;
    
    sigint = sqrt(sum1);

    return sigint;
}


double rknl_int_power(double rn){
    // nonlinear wavenumber k_{sigma} for power-law scale-free models, Del2 = (k/k_0)^(3+n). eqn. (55) of Smith et al., pg 15.
    // similarly, n_eff, (C7) 
      
    double rknl_int_pow, a, gammln, arg;
    
    int    ifail;
      
    arg = (3. + rn)/2.;
      
    // replaced call to NR gammaln. by gsl. computs natural log of gamma fn. using lanczos formula. This is calculating the required (1+n)/2 ! I'd imagine.   
    a   = exp(gsl_sf_lngamma(arg));
      
    rknl_int_pow = pow(0.5*a, -1./(3. + rn));
      
    return rknl_int_pow;
}


int rknlCalc(double* rknl, double* rneff, double* rncur){
    // bisection calc to find where sigma=1., as defined in eqns. C5 and C6 of Smith et al.
    // printf("\n\nComputing effective spectral quantities.");

    double xlogr1, xlogr2, rmid;
    double sig, d1, d2;
    
    double diff = pow(10., 99.);
    
    // range for bisection. 
    xlogr1 = -3.0;
    xlogr2 =  3.5; 
    
    while(fabs(diff)>0.001){ 
        rmid   = (xlogr2 + xlogr1)/2.0;
      
        rmid   = pow(10., rmid);
  
        // calculate k_nonlinear, n_eff and curvature at \sigma = 1., placing answers in sig, d1, d2.  
        wint(rmid, &sig, &d1, &d2);
  
        diff   = sig - 1.0;
      
        // having found (sigma-1.0) is close enough to 1 ...
        if(fabs(diff) < 0.001){
            *rknl  =  1./rmid;
        
            *rneff = -3. - d1;
            *rncur =     - d2;                  
        
            // we are done. 
            break;
        }
      
        else if(diff>0.001){
            xlogr1 = log10(rmid);
        }
      
        else if(diff<-0.001){
            xlogr2 = log10(rmid);
        }
        
        else{
            printf("\n\nError in halofit calc. of rknl, gradient and curvature.");
        
            return 1;
        }
    }

    printf("\nrknl [h/Mpc] = %.6e, rneff = %.6e, rncur = %.6e", *rknl, *rneff, *rncur);

    return 0;
}


int wint(double r, double* sig, double* d1, double* d2){
    // The subroutine wint finds the effective spectral quantities
    // r_knl, r_neff and r_ncur. rknl is the 'radius' of 
    // the Gaussian filter at which the variance is unity.
    // rneff is defined as the first derivative of the variance, calculated 
    // at the nonlinear wavenumber and similarly the rncur is the second
    // derivative at the nonlinear wavenumber, curvature. 

    double sum1, sum2, sum3, t, y, x, w1, w2, w3;
    double rk, del2;
    
    int    i, nint, ipflag;

    nint = 3000;
      
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
      
    for(i=1; i<nint; i++){
        t    = (i-0.5)/nint;
        y    = -1.0 + 1.0/t;
        rk   = y;
        
        del2 =  BunnWhite(rk); // p_lin(rk);                                                                      
        
        x    = y*r;
        
        w1   =            exp(-x*x);
        w2   =            2.*x*x*w1;
        w3   = 4.*x*x*(1. - x*x)*w1;
        
        sum1 += w1*del2/y/t/t;
        sum2 += w2*del2/y/t/t;
        sum3 += w3*del2/y/t/t;
    }
      
    sum1 = sum1/nint;
    sum2 = sum2/nint;
    sum3 = sum3/nint;
      
    *sig = sqrt(sum1);
      
    *d1  = -sum2/sum1;
    *d2  = -sum2*sum2/sum1/sum1 - sum3/sum1;
    
    // printf("\n\nSigma: %.2e, gradient: %.2e, curvature: %.2e", *sig, *d1, *d2);
    
    return 0;
}


int halofit_fittingfn(double rk, double rknl, double rn, double rncur, double om_mz, double om_vz, double plin, double* pnl){
    // halo model nonlinear fitting formula as described in Appendix C of Smith et al. (2003)
    // rk, rn, rncur, rknl, om_mz, om_vz AND plin are input parameters. output are onehalo term ph, two halo pq and pnl.  
    // rknl: k non-linear. rn first derivatitve, rncur: second derivative (curvature).  
    
    double gam, a, amod, b, c, xmu, xnu, alpha, beta, f1, f2, f3, f4;
    double f1a, f2a, f3a, f4a, f1b, f2b, f3b, f4b, frac;
      
    double pq, ph;
      
    // double y, p_cdm;
    // double om_m, om_v, om_b, h, p_index, gams, sig8, amp, om_mz, om_vz;
      
    // int    ipflag;

    gam = 0.86485 + 0.2989*rn + 0.1631*rncur;
      
    // eqn. (c9))
    a   = 1.4861 + 1.83693*rn + 1.67618*rn*rn + 0.7940*rn*rn*rn + 0.1670756*rn*rn*rn*rn - 0.620695*rncur;
    a   = pow(10., a);      
      
    // eqn. c10
    b   = pow(10., (0.9463  + 0.9466*rn + 0.3084*rn*rn - 0.940*rncur));
      
    c   = pow(10., (-0.2807 + 0.6669*rn + 0.3214*rn*rn - 0.0793*rncur));
      
    xmu = pow(10., (-3.54419 + 0.19086*rn));
      
    xnu = pow(10., (0.95897 + 1.2857*rn));
      
    alpha = 1.38848 + 0.3701*rn - 0.1452*rn*rn;
      
    beta  = 0.8291 + 0.9854*rn + 0.3400*pow(rn, 2.);

    // Here we need omega's at redshift of interest.
    // Omega dependent functions. C17 and C18. 
    
    if(fabs(1. - om_mz) > 0.01){ 
        f1a = pow(om_mz, -0.0732);
        f2a = pow(om_mz, -0.1423);
        f3a = pow(om_mz,  0.0725);
        
        f1b = pow(om_mz, -0.0307);
        f2b = pow(om_mz, -0.0585);
        f3b = pow(om_mz,  0.0743);     
        
        frac = om_vz/(1. - om_mz); 
        
        f1 = frac*f1b + (1.-frac)*f1a;
        f2 = frac*f2b + (1.-frac)*f2a;
        f3 = frac*f3b + (1.-frac)*f3a;
    }
    
    else{    
        f1 = 1.0;
        f2 = 1.;
        f3 = 1.;
    }
    
    y    = (rk/rknl);

    // eqn. C4.
    ph   = a*pow(y, 3.*f1)/(1. + b*pow(y, f2) + pow(f3*c*y, 3.-gam));
    
    // eqn. C3.
    ph   = ph/(1. + xmu*pow(y, -1.) + xnu*pow(y, -2.));
    
    // eqn. C2.
    pq   = plin*pow(1. + plin, beta)/(1. + plin*alpha)*exp(-y/4.0 - y*y/8.0);

    *pnl = pq + ph;

    return 0;
}


double BunnWhite(double rk){
    // COBE normalised Del^2(k), assuming Eisenstein & Hu transfer function. allowing for tilt. dh == \delta_H, eqn. (A1) Eisenstein & Hu.
    double tilt, dh, pk;
    
    tilt = p_index - 1.;

    if(Om_v > 0.){
      // cosmological constant. lambda = 1 - Om_m. eqn. (A3) 
      dh = 1.94*pow(10., -5.)*pow(Om_m, -0.785 - 0.05*log(Om_m))*exp(-0.95*tilt - 0.169*tilt*tilt);
    }  
    
    else{
      // lambda = 0.0, eqn. (A2).
      dh = 1.95*pow(10., -5.)*pow(Om_m, -0.35 - 0.19*log(Om_m) - 0.17*tilt)*exp(-tilt-0.14*tilt*tilt);
    }
                   // c/H_0
    pk = dh*dh*pow(3000.*rk, 3. + p_index)*EisensteinHu(rk, Om_m, Om_b)*EisensteinHu(rk, Om_m, Om_b);

    // when amp = 1., corrsponds to Del^2(k) at redshift zero. 
    return amp*amp*pk;
}


double EisensteinHu(double yy, double om_m, double om_b){
    //    "the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])", JAP, remember uses k/h, whereas they use pure k
    //**  om_m (equiv to Om_0 in E&Hu) is the total matter density parameter - i.e. CDM + baryons. **//
    
    double rk;
    double thet, b1, b2, zd, ze, rd, re, rke, s, rks, q, y, g;
    double ab, a1, a2, ac, bc, f, c1, c2, tc, bb, bn, ss;
    double tb, tk_eh;
    
    // h^-1 Mpc to Mpc I imagine. 
    rk   = yy*h;

    // doesn't seem necessary. 
    // e    = exp(1.);
      
    // T_CMB is 2.7 thet_{2.7} K.
    thet = 2.728/2.7;
    
    b1   = 0.313*pow(om_m*h*h, -0.419)*(1. + 0.607*pow(om_m*h*h, 0.674)); //eqn. (4)
    b2   = 0.238*pow(om_m*h*h, 0.223);
    zd   = 1291.*(1. + b1*pow(om_b*h*h,b2))*pow(om_m*h*h, 0.251)/(1. + 0.659*pow(om_m*h*h, 0.828));
    
    // z_equality. eqn. (2)
    ze  = 2.5*pow(10., 4.)*om_m*h*h/pow(thet, 4.);
    
    // R eqn. (5)
    rd  = 31500.*om_b*h*h/pow(thet, 4.)/zd;  // evaluated at redshift, zd.
    re  = 31500.*om_b*h*h/pow(thet, 4.)/ze;      // ''                   , ze.
    
    // k equality, eqn. (3)
    rke = 7.46*pow(10., -2.)*om_m*h*h/pow(thet, 2.);
    
    s   = (2./3./rke)*sqrt(6./re)*log(    (sqrt(1. + rd) + sqrt(rd + re))/(1. + sqrt(re))    ); // s equality, eqn (6) for R=re
    
    // kSilk eqn. (7)
    rks = 1.6*pow(om_b*h*h, 0.52)*pow(om_m*h*h, 0.73)*(1. + pow(10.4*om_m*h*h, -0.95) );

    // eqn. (10)
    q   = rk/13.41/rke;
      
    y   = (1. + ze)/(1. + zd);
    
    // G(y), eqn. 15. at y = (1+z_eq)/(1+z_d).
    g   =  y*(-6.*sqrt(1.+y) + (2.+3.*y)*log((sqrt(1.+y) + 1.)/(sqrt(1.+y) -1.) ));
    
    // alpha_b
    ab  = g*2.07*rke*s/pow(1.+rd, 0.75);

    a1 = pow(46.9*om_m*h*h, 0.670)*(1. + pow(32.1*om_m*h*h, -0.532));  // eqns. 11-12
    a2 = pow(12.0*om_m*h*h, 0.424)*(1. + pow(45.0*om_m*h*h, -0.582));
    
    // alpha_c
    ac = pow(a1, -om_b/om_m)*pow(a2, -pow(om_b/om_m, 3.));

    b1 = 0.944/(1. + pow(458.*om_m*h*h, -0.708));
    b2 = pow(0.395*om_m*h*h, -0.0266);
    // beta_c,              // O_c = om_m - om_b.
    bc = 1./(1. + b1*(  pow(1. - om_b/om_m, b2) -1.) );

    // eqn. (18)
    f  = 1./(1. + pow(rk*s/5.4, 4.));

    // eqn. 20, a_c=1.
    c1 = 14.2    + 386./(1.+69.9*pow(q,1.08));
    // eqn. (20)
    c2 = 14.2/ac + 386./(1.+69.9*pow(q,1.08));
    
    tc = f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) + (1-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q);

    // beta_b, eqn. (24).
    bb = 0.5 + (om_b/om_m) + (3.-2.*om_b/om_m)*sqrt(pow(17.2*om_m*h*h, 2.) + 1.);
    
    // beta_node eqn. (23)
    bn = 8.41*pow(om_m*h*h, 0.435);
    
    // eqn (22).
    ss = s/pow(1. + pow(bn/rk/s, 3.), 1./3.);
    
    // \tilde T(k; 1,1), eqn. 19
    tb = log(e + 1.8*q)/(log(e + 1.8*q) + c1*q*q)/(1. + pow(rk*s/5.2, 2.));
    
    // eqn. 21 baryon transfer function. 
    tb = (tb + ab*exp(-pow(rk/rks, 1.4))/(1. + pow(bb/rk/s, 3.)))*sin(rk*ss)/rk/ss;

    // eqn. (8)
    tk_eh = (om_b/om_m)*tb + (1.-om_b/om_m)*tc;

    return tk_eh;
}
