double splint_NFWinversion_rq(double q){
    double Interim;

    splint(q_Nfwinversion, r_Nfwinversion, rq2D_Nfwinversion, 1000, q, &Interim);
    
    return Interim;
}


double NFW_densityprofile_NormedFourier(double k){
    // For u(k|m) = (1/m)\int_Vol \rho(x|m) exp^{-ikx}

    // Spherically symmetric profile, truncated at the Virial radius. 

    // u(k|m) = \int_0^{r_{vir}} dr 4pi r^2 [sin(kr)/kr] rho(r|m)/m
    
    double      Interim; 
    double           rs;
    
    rs       = NFW_rvir/NFW_conc;
    
    printf("\n%e \t %e", NFW_rvir, NFW_conc);
    
    Interim  = sin(k*rs)*(gsl_sf_Si((1. + NFW_conc)*k*rs) - gsl_sf_Si(k*rs)) - sin(NFW_conc*k*rs)*pow((1.+ NFW_conc)*k*rs, -1.) + cos(k*rs)*(gsl_sf_Ci((1.+ NFW_conc)*k*rs) - gsl_sf_Ci(k*rs));

    // for the NFW profile, defined by rhos, rs and c, mass is not independent. 
    Interim /= log(1. + NFW_conc) - NFW_conc/(1.+ NFW_conc);

    return Interim;
}


double haloModel_pk(double k, double beta, int Order){
    // assumes linear bias of 1, therefor beta = f. Here 2. is the standard deviation of the gaussian variate added to z co-ordinate for non-linear RSD.

    double Interim;
    double Dplus = 0.5;
    
    Interim = TotalVolume/haloNumber; 
                                                                  
    return Interim*pow(NFW_densityprofile_NormedFourier(k), 2.)*kaiserGauss_multipole(k*2., 0.0, Order) + pow(Dplus, 2.)*(*pt2Pk)(k)*kaiserGauss_multipole(k*2., beta, Order);
}


int prep_NFWhaloCat(int total_gal){
    double rs;

    rhobox = 2.775*pow(10., 11.)*Om_m;             // units of h^-1 Mstar (h^-1 Mpc)^{-3}
    Mbox   = 2.775*pow(10., 11.)*Om_m*TotalVolume; // units of h^-1 Mstar.
 
    Mhalo  = Mbox/haloNumber;
    Mpart  = Mbox/total_gal;
    
    // fixed by volume of the box, Cosmology and further specified halo comoving density. 
    printf("\n\nResolution in halo     mass: %e [10^11 h^-1 Mstar]", Mhalo*pow(10., -11.));
    printf("\nResolution in particle mass: %e [10^11 h^-1 Mstar]",   Mpart*pow(10., -11.));
    
    NFW_rvir = pow(3.*Mhalo/(4.*pi*Delta_crit*rhobox), 1./3.);

    r_Nfwinversion    = malloc(1000*sizeof(*r_Nfwinversion));        // r_Nfwinversion[j] is in units of r_s
    q_Nfwinversion    = malloc(1000*sizeof(*q_Nfwinversion));
    rq2D_Nfwinversion = malloc(1000*sizeof(*rq2D_Nfwinversion));
    
    rs       = NFW_rvir/NFW_conc;
    
    for(j=0; j<1000; j++){  
        r_Nfwinversion[j] = (j/1000.)*(NFW_rvir/rs); // r_Nfwinversion[j] is in units of r_s
    
        q_Nfwinversion[j] = (log(1.+ r_Nfwinversion[j]) - r_Nfwinversion[j]*pow(1. + r_Nfwinversion[j], -1.))*pow(log(1. + NFW_conc) - NFW_conc/(1. + NFW_conc), -1.);
    }
    
    spline(q_Nfwinversion, r_Nfwinversion, 1000, 1.0e31, 1.0e31, rq2D_Nfwinversion);
    
    rand_number   = total_gal;
    
    printf("\n\nTotal Gals in catalogue: %d", rand_number);
    
    assign_randmemory();
    
    hod_disp = malloc(rand_number*sizeof(double));
    
    prep_DisplacementCalc();

    return 0;
}


int HaloCatalogue_NFWprofile(int total_gal){
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
}


int NFW_NGP(){
    for(j=0; j<n0*n1*n2; j++) densityArray[j] = 0.0;
    
    // printf("\nmin: %e \t %e", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    
    for(j=0; j<rand_number; j++){
        boxlabel = boxCoordinates(rand_x, rand_y, rand_z, j);
    
        // printf("\nbox label: %d", boxlabel);
    
        densityArray[boxlabel] += 1.;
    }
    
    for(j=0; j<n0*n1*n2; j++) densityArray[j] /= rand_number;
        
    return 0;
}


int NFW_NGP_overdensity(){
    for(j=0; j<n0*n1*n2; j++) densityArray[j] = 0.0;
    
    for(j=0; j<rand_number; j++){
        boxlabel = boxCoordinates(rand_x, rand_y, rand_z, j);
    
        densityArray[boxlabel] += 1.;
    }
    
    for(j=0; j<n0*n1*n2; j++) densityArray[j] /= CellVolume;
        
    // nbar
    for(j=0; j<n0*n1*n2; j++) densityArray[j] /= (rand_number/TotalVolume);
    
    for(j=0; j<n0*n1*n2; j++) densityArray[j] -= 1.;
        
    return 0;
}


int ukm_calc(){    
     // NFW_NGP();
    NFW_NGP_overdensity();
    
    // The true density field multiplied by a mask, W(x). 
    for(j=0; j<n0*n1*n2; j++) in[j][0] = densityArray[j];
    
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    fftw_execute(p);
    
    polarPk_modeCount = 0;
        
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;
                
                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(kSq, 0.5);
                
                mu                                 = k_z/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;      
                
                H_kReal                            = out[Index][0];
                H_kImag                            = out[Index][1];
                
                // pk calc. 
                H_kReal                           /= n0*n1*n2;
                H_kImag                           /= n0*n1*n2;
                
                WindowFunc                         = 1.;

                if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/NyquistWaveNumber)/(pi*k_x*0.5/NyquistWaveNumber);
                if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/NyquistWaveNumber)/(pi*k_y*0.5/NyquistWaveNumber);
                if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/NyquistWaveNumber)/(pi*k_z*0.5/NyquistWaveNumber);

                // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2. 
                H_kReal                           /= pow(WindowFunc, 1.);
                H_kImag                           /= pow(WindowFunc, 1.);
                                
                PkArray[Index][0]                  = kmodulus;
                PkArray[Index][1]                  = pow(H_kReal, 2.) + pow(H_kImag, 2.);
                
                // pk calc. 
                PkArray[Index][1]                 *= TotalVolume; 
                
                // PkArray[Index][1]                 -= (TotalVolume/rand_number);

                // Issue with mu for a zeroth length vector being ill defined. 
	            if(kmodulus > 0.000001){
	                // Only half the modes are independent. 
	            	if(k_z>0.){
	            	    // One hemi-sphere is independent, e.g. k_z >= 0.
		                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		                polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		                polar2Dpk[polarPk_modeCount][2]    = PkArray[Index][1];
		            
		                polarPk_modeCount                 += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y > 0.0)){
                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		                polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		                polar2Dpk[polarPk_modeCount][2]    = PkArray[Index][1];
		            
		                polarPk_modeCount                 += 1;
		            }
		            
		            else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
		                // on the line k_z=k_y=0, one half is independent, k_x>=0.
		                                        // in the k_z=0 plane one semi-circle is independent, k_y>0.
		                polar2Dpk[polarPk_modeCount][0]    = kmodulus;
		                polar2Dpk[polarPk_modeCount][1]    = fabs(mu);
		                polar2Dpk[polarPk_modeCount][2]    = PkArray[Index][1];
		            
		                polarPk_modeCount                 += 1;
		            }    
	            }
	        }
        }
    }
    
    sprintf(filepath,"%s/Data/stacpolly/NFW_profile_500_pk_%d.dat", root_dir, loopCount);

    Monopole(filepath);
    
    sprintf(filepath,"%s/Data/stacpolly/NFW_profile_500_quad_%d.dat", root_dir, loopCount);
    MultipoleCalc(kBinNumb, meanKBin, kMonopole, kQuadrupole, polar2Dpk, polarPk_modeCount, filepath, kbinInterval, 0.0, 1.0, 1);
    
    return 0;
}


int NFWprofile_oneHalo_pairCount(){
    // particles tracing halo: 100218
    
    assignMemory_xi();

    grow_randTree();

    CountPairs_rMu(rr_0, rr_2, rr_4, rr_meanr, rr_meanmu, randTree, randTree, 1);

    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_onehalo_DD.dat", root_dir);

    output = fopen(filepath, "w");

    for(j=0; j<nlogbins; j++){  
        for(i=0; i<nlinbins; i++)  fprintf(output, "%e \t", rr[j][i]);
    
        fprintf(output, "\n");
    }

    fclose(output);

    return 0;
}


int NFWprofile_oneHalo_xiCalc(){
    // particles tracing halo: 100218
    sprintf(surveyType, "NFW_profile_onehalo_xi");
    
    assignMemory_xi();

    // load 'DD'
    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_onehalo_DD.dat", root_dir);

    inputfile = fopen(filepath, "r");

    for(j=0; j<nlogbins; j++){
      for(i=0; i<nlinbins; i++)  fscanf(inputfile, "%le \t", &gg[j][i]);

      fscanf(inputfile, "\n");
    }

    fclose(inputfile);

    Vipers_Num = rand_number;

    // load RR.
    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_RR.dat", root_dir);

    inputfile = fopen(filepath, "r");

    for(j=0; j<nlogbins; j++){
      for(i=0; i<nlinbins; i++)  fscanf(inputfile, "%le \t", &rr[j][i]);

      fscanf(inputfile, "\n");
    }

    fclose(inputfile);

    rand_number = 1181096;
    
    
    landy_szalay();

    xiMonopole(landy_xi, gg_meanmu, xi0);
    
    
    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_xi_multipoles_mean.dat", root_dir);

    output = fopen(filepath, "w");

    for(j=0; j<nlogbins; j++){  
        fprintf(output, "%e \t %e \n", logrbins[j], xi0[j]);
    }
  
    fclose(output);

    return 0;
}

/*
int NFW_profileOneHalo_xi(){
    FFTlogRes        = 4096;
    
    ** Change **
    // pt2Pk            = &ukm2;
    
    FFTlog_memory(FFTlogRes, beta, velDispersion);
    
    xi_mu(mono_config);
    
    
    sprintf(filepath, "%s/Data/stacpolly/NFW_profile_fftlog_multipoles_xi.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<mono_config->N; j++){  
        if((mono_config->krvals[j][1] > 0.1) && (mono_config->krvals[j][1] < 100.)){
            fprintf(output, "%e \t %e \n", mono_config->krvals[j][1], mono_config->xi[j][0]);
        }
    }
    
    fclose(output); 
    
    printf_u_ofrm();

    return 0;
}
*/
