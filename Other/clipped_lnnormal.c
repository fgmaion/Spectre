double J0_calc(double u0, double sigma){
    return 0.5*(1. + gsl_sf_erf(u0 - sigma/sqrt(2.)));
}


double recursion(double J_lo, int n, double u0, double sigma){
    return J_lo - HermitePolynomial(u0, n-1)*exp(-pow(u0 - sigma/sqrt(2.), 2.))*pow(2., -n/2.)*pow(sigma, -1.*n)/sqrt(pi);
}


double series_term(double xi_g, double Jn, int n){
    return pow(xi_g, n+1.)*Jn*Jn/gsl_sf_fact(n+1);
}


double dG0_given_d0(double d0, double sigma){
    // fractional Gaussian threshold, u0, given the variance of the underlying Gaussian and the threshold applied to the lognormal field: d0. 
    return sigma*sigma/2. + log(1. + d0);
}


double u_0_given_dG0(double dG0, double sigma){
    // fractional Gaussian threshold, u0, given the variance of the underlying Gaussian and the threshold applied to the lognormal field: d0. 
    return dG0/sqrt(2.)/sigma;
}


double f_given_u0(double u0, double sigma){
    // amplitude suppression of the clipped lognormal model: \xi_C = f^2 \xi_G given the normalised threshold u0. 
    return 0.5*(1. + gsl_sf_erf(u0 - sigma/sqrt(2.)));
}


int print_amplitude_threshold_relation(){
    double d0, dG0, u0, f;
    
    // double sig = underlyingGaussian_sigma;
    
    double sig = 1.30;
    
    sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/amplitude_threshold_sig_%.2e.dat", root_dir, sig);
    
    output = fopen(filepath, "w");

    for(j=0; j<1000; j++){
        d0  = 0.01*j - 1.;
    
        dG0 = dG0_given_d0(d0, sig);
    
        u0  = u_0_given_dG0(dG0, sig);
    
        f   = f_given_u0(u0, sig);
        
        fprintf(output, "%e \t %e \t %e \n", d0, u0, f);
    }

    fclose(output);

    return 0;
}


int kvals_pairup(int mono_order, double kVals[]){
    // Given the k values of the measured P(k) multipoles, find closest fftlog modes. 
    double     diff;
    double min_diff;
    
    fftlog_indices = malloc(mono_order*sizeof(*fftlog_indices));
    
    
    for(i=0; i<mono_order; i++){  
        min_diff = pow(10., 99.);
    
        for(j=0; j<FFTlogRes; j++){  
            diff = fabs(mono_config->krvals[j][0] - kVals[i]);
            
            if(diff<min_diff){
                min_diff = diff;
            
                fftlog_indices[i]                 = j;
            }
        }
        
        if(100.*min_diff/kVals[i] > 1.)  printf("\nError matchup in k, %.2e percentage match in k of %.2e", 100.*min_diff/kVals[i], kVals[i]);
    }
    
    printf("\n\nk vals match up.\n");
    
    printf("\nMatch up between observed k vals and FFTlog vals complete.\n");
    
    // for(j=0; j<mono_order; j++)  printf("%e \t %e \n", kVals[j], mono_config->krvals[fftlog_indices[j]][0]);
    
    return 0;
}


double regression_calc(){
    sprintf(filepath,"%s/W1_Nagoya_v6_mocks_work/clipped_lnnormal/HOD_cube_256_pk_galweights_depletionfactor_%.2lf_foldfactor_%.1lf_d0_%.1lf_R_%.1lf.dat", root_dir, depletion_factor, Jenkins_foldfactor, appliedClippingThreshold, clipping_smoothing_radius);
    
    inputfile = fopen(filepath, "r");
    
    int line_no;
    
    ch         = 0;
    line_no    = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  line_no += 1;
    } while (ch != EOF);

    rewind(inputfile); 
    
    double        kvals[line_no];
    double clipped_mono[line_no];
    
    for(j=0; j<line_no; j++)   fscanf(inputfile, "%lf \t %lf \t %*d \n", &kvals[j], &clipped_mono[j]);
    
    // for(j=0; j<line_no; j+=2)  printf("%lf \t %le \n", kvals[j], clipped_mono[j]);
    
    fclose(inputfile);
    

    // initialise FFTlog.
    FFTlogRes     =               4096;
    
    fsigma8       =                0.0;  
    bsigma8       =        0.557*1.495;
    velDispersion =               0.00;
    
    FFTlog_memory(FFTlogRes, fsigma8/bsigma8, velDispersion);
    
    xi_mu(mono_config);
    
    // non-linear halofit p(k) to p(k) of the underlying Gaussian field in the lognormal model.  
    invert_lnnorm(mono_config, mono_config);

    pk_mu(mono_config);
    // mono_config-> pk now corresponds to the p(k) of the underlying Gaussian field in the 
    // lognormal model.  
    
    // variance of the underlying Gaussian. 
    varCalc(mono_config, &variance);
    
    kvals_pairup(line_no, kvals);
    
    // The regression calc.
    double sum_Di_xiG;
    double sum_xiG_sq;
    
    double f;
    
    for(i=0; i<line_no; i++){
        if((kvals[i] >= 0.01) && (kvals[i] <= 0.06)){
            sum_Di_xiG += clipped_mono[i]*mono_config->pk[fftlog_indices[i]][0];
            
            sum_xiG_sq += mono_config->pk[fftlog_indices[i]][0]*mono_config->pk[fftlog_indices[i]][0];
        }
    }
    
    f = sqrt(sum_Di_xiG/sum_xiG_sq);
    
    printf("\n\nf found to be: %e", f);
    
    return f;
} 


int Gauss_varlognormal(){
    // after smoothing, calculate the variance of the underlying Gaussian field < log(1 + delta_ln) > = -sig^2 /2. 
    // assuming delta field is lognormal.  Smoothed on scale GaussianFilter_radius.
    underlyingGaussian_sigma = 0.0;
 
    double mean_density;
 
    mean_density = Vipers_Num/TotalVolume;
 
    prep_grid();
 
    // clean overdensity before galaxy assignment. 
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;
 
    // assign galaxies using NGP.  overdensity contains the galaxy counts. 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){ 
            boxlabel                          = boxCoordinates(xCoor, yCoor, zCoor, j);
            
            overdensity[boxlabel][0]         += 1.;
        }
    }
  
    for(j=0; j<n0*n1*n2; j++){
        // counts 
        overdensity[j][0] /= (CellVolume*mean_density);
        
        overdensity[j][0] -= 1.;
    } 
 
    // Second flag to enforce zero mean for the field, 1 for true. 
    // smoothing does not affect the mean, homogeneous mode, in general.
    Gaussian_filter(clipping_smoothing_radius, 0);
 
    // <1 + delta> = 1: (1 + \delta) should be rescaled, not subtracted.
    double exp_onedelta = 0.0;
    
    for(j=0; j<n0*n1*n2; j++)  exp_onedelta += (1. + overdensity[j][0]);
    
    exp_onedelta /= n0*n1*n2;
    
    printf("\nExpectation of 1 + delta: %.4lf", exp_onedelta);
    
 
    for(j=0; j<n0*n1*n2; j++)  underlyingGaussian_sigma += log(1. + overdensity[j][0]);
    
    underlyingGaussian_sigma /= n0*n1*n2;
    
    underlyingGaussian_sigma  = sqrt(-2.*underlyingGaussian_sigma);
    
    printf("\n\n\sig = %.4lf", underlyingGaussian_sigma);
    
    
    for(j=0; j<n0*n1*n2; j++){
        if(overdensity[j][0]  > appliedClippingThreshold){
           overdensity[j][0]  = appliedClippingThreshold; 
        }
    }
    
    return 0;
}


int clipped_lnnormal_main(){
    print_amplitude_threshold_relation();
    
    // regression_calc();

    return 0;
}
