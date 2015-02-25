int inputLinearPk(){
    // Interpolated theoretical P(k) on a regular grid in k. 
    lineark  = realloc(lineark,          470*sizeof(*lineark));
    linearPk = realloc(linearPk,         470*sizeof(*linearPk));
              
    // Second derivates of HOD P(k) for cubic spline.           
    linear2d = realloc(linear2d,         470*sizeof(*linear2d));

    sprintf(filepath, "%s/Data/HODTheoryPk/camb_matterPk.dat", root_dir);
    
    inputfile = fopen(filepath, "r");
    
    for(j=0; j<470; j++) fscanf(inputfile, "%le \t %le \n", &lineark[j], &linearPk[j]);
    
    for(j=0; j<470; j++) linearPk[j] *= linearBias*linearBias;
    
    fclose(inputfile); 
    
    
    spline(lineark, linearPk, 470, 1.0e31, 1.0e31, linear2d);
    
    pt2Pk = &linearPk_Gaussian;
    
    // Approx. model to HOD P(k), inflationary power law with exponential truncation. 
    pt2Xi = &Pk_powerlaw_truncated_xi;
    
    sprintf(theoryPk_flag, "linear_-20.0");
    
    return 0;
}


double splintLinearPk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    
    // if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 1.00)*pow(linearBias/1.495903, 2.); 
    if(k<0.0004)  return 4.675*pow(10., 6.)*pow(k, 0.96)*pow(linearBias/1.495903, 2.); 

    else{
        splint(lineark, linearPk, linear2d, 469, k, &Interim);
    
        return Interim*pow(linearBias/1.495903, 2.);
    }
}


double linearPk_Gaussian(double k){
    return splintLinearPk(k)*exp(-0.5*pow(3.*k, 2.));
}


double splintHODpk(double k){
    // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec). 
    if(k>20.)          return pk_hiA*pow(k, 3. + pk_hin);
    
    else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon); 

    else{
        double Interim;
    
        splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);
    
        return Interim;
    }
}


double sigma8_integrand(double k){
    return pow(2.*pi, -3.)*splintHODpk(k)*4.*pi*k*k*spherical_tophat(k, 8.)*spherical_tophat(k, 8.);
}


double sigma8_calc(){
    double Interim;
    
    Interim = qromb(&sigma8_integrand, 0.0000001, 200.);
    
    return sqrt(Interim);
}


int inputHODPk(){
    sprintf(filepath, "%s/Data/500s/EisensteinHu_halofit_pk_%.2f.dat", root_dir, z_eff);
    
    inputfile  = fopen(filepath, "r");
    
    ch         = 0;
    pk_lineNo  = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  pk_lineNo += 1;
    } while (ch != EOF);

    rewind(inputfile);
    
    // Interpolated theoretical P(k) on a regular grid in k. 
    sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
    sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
    
    // Second derivates of HOD P(k) for cubic spline.           
    sdlt2d = realloc(sdlt2d,         pk_lineNo*sizeof(*sdlt2d));

    for(j=0; j<pk_lineNo; j++)       fscanf(inputfile, "%le \t %*le \t %le \n", &sdltk[j], &sdltPk[j]);
    
    // for(j=0; j<pk_lineNo; j++)  printf("\n%e \t %e", sdltk[j], sdltPk[j]);
    
    fclose(inputfile);
    
    
    spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);
        

    powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);
    
    powerlaw_regression(pk_lineNo, 20., 100., 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

    
    pt2Pk      = &splintHODpk;
    
    aexp       = 1./(1. + z_eff);  
    
    app_sigma8 = 0.0;

    app_sigma8 = sigma8_calc()/linearGrowth_factor(log(aexp));
    
    printf("\n\napp_sigma_8(%.2e): %.4e, sigma_8(0.0): %.4e", z_eff, sigma8_calc(), app_sigma8);
    
    return 0;
}


double HODPk_Gaussian(double k){
    return  0.02*splintHODpk(k)*exp(-0.5*pow(3.*k, 2.));
}
