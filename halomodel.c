double splint_nfw_inversion(double q){
    double Interim;

    splint(q_nfwinversion, r_nfwinversion, r2D_nfwinversion, 1000, q, &Interim);
    
    return Interim;
}


double ukm_nfw_profile(double k){
    // For u(k|m) = (1/m)\int_Vol \rho(x|m) exp^{-ikx}
    // Spherically symmetric profile, truncated at the virial radius. 
    // u(k|m) = \int_0^{r_{vir}} dr 4pi r^2 [sin(kr)/kr] rho(r|m)/m
    
    double      Interim; 
    double           rs;
    
    rs       = nfw_rvir/nfw_conc; 
    
    Interim  = sin(k*rs)*(gsl_sf_Si((1. + nfw_conc)*k*rs) - gsl_sf_Si(k*rs)) - sin(nfw_conc*k*rs)*pow((1.+ nfw_conc)*k*rs, -1.) + cos(k*rs)*(gsl_sf_Ci((1.+ nfw_conc)*k*rs) - gsl_sf_Ci(k*rs));

    // for the NFW profile, defined by rhos, rs and c, mass is not independent. 
    Interim /= log(1. + nfw_conc) - nfw_conc/(1.+ nfw_conc);

    return Interim;
}


double haloModel_pk(double k, double nbar, double beta, int Order){
    // assumes linear bias of 1, therefor beta = f. Here 2. is the standard deviation of the gaussian variate added to z co-ordinate for non-linear RSD.

    double Interim;
    double Dplus = 0.5;
    
    // sets nbar. 
    Interim = nbar; 
                                                                  
    return Interim*pow(ukm_nfw_profile(k), 2.)*(*pt2RSD_k)(k*2., 0.0, Order) + pow(Dplus, 2.)*(*pt2Pk)(k)*(*pt2RSD_k)(k*2., beta, Order);
}


int prep_NFWhaloCat(int halo_number, int total_gal){
    // halo catalogue construction. 
    
    double rs;
    
    pt2RSD_k = &kaiserLorentz_multipole;

    rhobox   = 2.775*pow(10., 11.)*Om_m;             // units of h^-1 Mstar (h^-1 Mpc)^{-3}
    Mbox     = 2.775*pow(10., 11.)*Om_m*TotalVolume; // units of h^-1 Mstar.
 
    Mhalo    = Mbox/halo_number;
    Mpart    = Mbox/total_gal;
    
    // fixed by volume of the box, cosmology and specified halo comoving density. 
    printf("\n\nResolution in halo     mass: %e [10^11 h^-1 Mstar]", Mhalo*pow(10., -11.));
    printf("\nResolution in particle mass: %e [10^11 h^-1 Mstar]",   Mpart*pow(10., -11.));
    
    nfw_rvir = pow(3.*Mhalo/(4.*pi*Delta_crit*rhobox), 1./3.);

    r_nfwinversion    = malloc(1000*sizeof(*r_nfwinversion));        //** r_nfwinversion[j] is in units of r_s **//
    q_nfwinversion    = malloc(1000*sizeof(*q_nfwinversion));
    r2D_nfwinversion  = malloc(1000*sizeof(*r2D_nfwinversion));
    
    rs       = nfw_rvir/nfw_conc;
    
    for(j=0; j<1000; j++){  
        r_nfwinversion[j] = (j/1000.)*(nfw_rvir/rs);                //** r_nfwinversion[j] is in units of r_s
    
        q_nfwinversion[j] = (log(1.+ r_nfwinversion[j]) - r_nfwinversion[j]*pow(1. + r_nfwinversion[j], -1.))*pow(log(1. + nfw_conc) - nfw_conc/(1. + nfw_conc), -1.);
    }
    
    spline(q_nfwinversion, r_nfwinversion, 1000, 1.0e31, 1.0e31, r2D_nfwinversion);
    
    // rand_number   = total_gal;
    
    // printf("\n\nTotal Gals in catalogue: %d", rand_number);
    
    // assign_randmemory();
    
    // hod_disp = malloc(rand_number*sizeof(double));
    
    // prep_DisplacementCalc();

    return 0;
}

/*
int HaloCatalogue_NFWprofile(int total_gal){
    // halo catalogue construction. 

    int    gals_perhalo;
    
    double x, y, z;
    double r, q, mu, phi, theta, midr;
    
    double xhalo, yhalo, zhalo;
    int    xlabel, ylabel, zlabel;

    // NFW profile    
    double alpha = 1., beta = 2.;
    double rs, rhos;
    double dispersion;
    
    rs       = NFW_rvir/NFW_conc;
    
    // Halo resolution determines, 
    rhos = Mhalo/(4.*pi*rs*rs*rs*(log(1.+NFW_conc) - NFW_conc/(1.+ NFW_conc)));
 
    printf("\n\nr_vir %e [h^-1 Mpc], rs %e [h^-1 Mpc], rho_s %e [10^11 h^-1 Mstar (h^-1 Mpc)^3]\n", NFW_rvir, rs, rhos*pow(10., -11.));
    
    gals_perhalo = (int) ceil(Mhalo/Mpart);
    
    printf("\n\n%d Gals per halo", gals_perhalo);
    
    int    randCount = 0;
    
    DisplacementCalc();
    
    for(k=0; k<haloNumber; k++){    
        xhalo        = 500.*gsl_rng_uniform(gsl_ran_r);
        yhalo        = 500.*gsl_rng_uniform(gsl_ran_r);
        zhalo        = 500.*gsl_rng_uniform(gsl_ran_r);  
                
        xlabel       = (int) floor(xhalo/CellSize);
        ylabel       = (int) floor(yhalo/CellSize);
        zlabel       = (int) floor(zhalo/CellSize);
        
        boxlabel     = xlabel + n2*ylabel + n2*n1*zlabel;
        
        xhalo       +=     Dplus*xDisplacement[boxlabel];
        yhalo       +=     Dplus*yDisplacement[boxlabel];
        zhalo       +=     Dplus*zDisplacement[boxlabel]; 
        
        // linear RSD, (1+f)*\psi_z.
        zhalo       += fGR*Dplus*zDisplacement[boxlabel]; 
        
        for(i=0; i<gals_perhalo; i++){
            hod_disp[randCount] = (1. + fGR)*Dplus*zDisplacement[boxlabel];
                
            q = gsl_rng_uniform(gsl_ran_r);
            
            r = splint_NFWinversion_rq(q);
            
            r *= rs;
            
            phi = 2.0*pi*gsl_rng_uniform(gsl_ran_r);
            
            mu  = 2.*(gsl_rng_uniform(gsl_ran_r) - 0.5);            
            
            theta = acos(mu);
            
            x   =  r*sin(theta)*cos(phi) + xhalo;
            y   =  r*sin(theta)*sin(phi) + yhalo;
            z   =  r*mu                  + zhalo;
            
            dispersion = rand_gaussian(gsl_ran_r, 2.);
            
            z  +=  dispersion;
            
            hod_disp[randCount] += dispersion;

            periodic_bounds(&x, 2);
            periodic_bounds(&y, 1);
            periodic_bounds(&z, 0);
            
            rand_x[randCount] = x;
            rand_y[randCount] = y;
            rand_z[randCount] = z;
                
            randCount        += 1;
        }
    }
    
    
    sprintf(filepath, "%s/Data/stacpolly/NFW_halo_catalogue_%d.dat", root_dir, loopCount);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<rand_number; j++) fprintf(output, "%e \t %e \t %e \t %e\n", rand_x[j], rand_y[j], rand_z[j], hod_disp[j]);
    
    fclose(output);    
    
    return 0;
}*/
