double HubbleCnst(double z){
    // Units of h km s^-1 Mpc^-1   
    return 100.*pow(Om_v + Om_m*pow(1.+z, 3.), 0.5);
}


double z_chi_integrand(double z){
    return pow(Om_v + Om_m*pow(1.+ z, 3.) + Om_r*pow(1.+z, 4.) - (Om_tot - 1.)*pow(1.+z, 2.), -0.5);
}


int chi_zcalc(){
    sprintf(filepath, "%s/Data/zComovingDistance/500s_hodmocks_params_%.2lf_%.2lf_%.2lf.dat", root_dir, Om_m, Om_v, Om_b);
    
    inputfile = fopen(filepath, "rb");
  
    if(inputfile == NULL){
        printf("\nCreating new z - comoving distance interpolation file.");
      
        pt2zChiIntegrand = &zChi_Integrand;
        
        // R_0*r = integral (c/H_0)*pow(Om_Lambda + Om_m*(1+z)^3 + Om_r*(1+z)^4 - (Om_tot -1)*(1+z)^2, -0.5) dz 
        for(i=1000; i>0; i--)  z_Array[1000-i]            =  2.0 + pow(1000.0, -1.0)*i*(0.0 - (2.0));
        
        // redshift 0.
        ComovingDistance_z[0] = 0.0;
        
        for(i=999; i>0; i--)  ComovingDistance_z[1000-i]  =  pow(100.0/lightSpeed_kmpersec, -1.0)*qromb(pt2zChiIntegrand, 0.0, z_Array[1000-i]);            
        
        // Test. 
        // qromb(pt2zChiIntegrand, 0.0, UpperIntegralLimit[i]) in units of [H_0*R_0/c] converted to [h^-1 Mpc] by pow(100.0/lightSpeed_kmpersec, -1.0) factor.

        printf("\nComoving distance to redshift 1.5: %e Mpc\n\n", pow(100.0*h/lightSpeed_kmpersec, -1.0)*qromb(pt2zChiIntegrand, 0.0, 1.5));
        
        output       = fopen(filepath, "wb");

        fwrite(z_Array,            sizeof(z_Array[0]),            sizeof(z_Array)/sizeof(z_Array[0]),                       output);
      
        fwrite(ComovingDistance_z, sizeof(ComovingDistance_z[0]), sizeof(ComovingDistance_z)/sizeof(ComovingDistance_z[0]), output);
      
        fclose(output);
    }

    else{
        // written as a binary file. 
        printf("\n\nReading z - comoving distance interpolation file.");
      
        fread(z_Array, sizeof(z_Array[0]), sizeof(z_Array)/sizeof(z_Array[0]), inputfile);
      
        fread(ComovingDistance_z, sizeof(ComovingDistance_z[0]), sizeof(ComovingDistance_z)/sizeof(ComovingDistance_z[0]), inputfile);
      
        fclose(inputfile);
    }
  
    // Working correctly, tested against Ned Wright Cosmology calculator.

    // First array must be a monotonically increasing function, start from redshift zero rather than redshift 2.0
    spline(z_Array, ComovingDistance_z, nPoints, 1.0e31, 1.0e31, z_ComovingDistance_2derivatives);
    
    spline(ComovingDistance_z, z_Array, nPoints, 1.0e31, 1.0e31, ComovingDistance_z_2derivatives);

    UniverseAge();

    linearGrowthRate();  // Must match declaration in cosmology_planck2015.h (cosmology_valueaddedmocks.h); This is NOT automatically ensured. //  
    
    return 0;
}


double interp_comovingDistance(double z){
    double InterimInterp_yVal;
    
    splint(z_Array, ComovingDistance_z, z_ComovingDistance_2derivatives, nPoints, z, &InterimInterp_yVal);
    
    return InterimInterp_yVal;  // Returns comoving distance at redshift z in h^-1 Mpc. 
}


double interp_inverseComovingDistance(double r){
    double InterimInterp_yVal;
    
    splint(ComovingDistance_z, z_Array, ComovingDistance_z_2derivatives, nPoints, r, &InterimInterp_yVal);
    
    return InterimInterp_yVal;  // Returns z at comoving distance, [h^-1 Mpc].
}
