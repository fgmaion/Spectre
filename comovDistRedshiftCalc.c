int comovDistReshiftCalc(){
    sprintf(filepath, "%s/Data/zComovingDistance/HODmock_params.dat", root_dir);
    inputfile = fopen(filepath, "rb");
  
    if(inputfile == NULL){
      printf("\nCreating new z - Comoving distance interpolation file.");

      // R_0*r = integral (c/H_0)*pow(Om_Lambda + Om_m*(1+z)^3 + Om_r*(1+z)^4 - (Om_tot -1)*(1+z)^2, -0.5) dz 
      for(i=1000; i>0; i--)  z_Array[1000-i]            =  2.0 + pow(1000.0, -1.0)*i*(0.0 - (2.0));
      for(i=1000; i>0; i--)  ComovingDistance_z[1000-i] =  pow(100.0/lightSpeed_kmpersec, -1.0)*qromb(pt2Integrand, 0.0, z_Array[1000-i]);            
    
      //Test. 
      //qromb(pt2Integrand, 0.0, UpperIntegralLimit[i]) in units of [H_0*R_0/c] converted to [h^-1 Mpc] by pow(100.0/lightSpeed_kmpersec, -1.0) factor.

      printf("\nComoving distance to redshift 1.5: %f Mpc\n\n", pow(100.0*h/lightSpeed_kmpersec, -1.0)*qromb(pt2Integrand, 0.0, 1.5));

      output       = fopen(filepath, "wb");

      fwrite(z_Array,            sizeof(z_Array[0]),            sizeof(z_Array)/sizeof(z_Array[0]),                       output);
      fwrite(ComovingDistance_z, sizeof(ComovingDistance_z[0]), sizeof(ComovingDistance_z)/sizeof(ComovingDistance_z[0]), output);
      fclose(output);
    }

    else{
      printf("\n\nReading z - Comoving distance interpolation file.");
      fread(z_Array, sizeof(z_Array[0]), sizeof(z_Array)/sizeof(z_Array[0]), inputfile);
      fread(ComovingDistance_z, sizeof(ComovingDistance_z[0]), sizeof(ComovingDistance_z)/sizeof(ComovingDistance_z[0]), inputfile);
      fclose(inputfile);
    }
  
    // Working correctly, Ned Wright Cosmology calculator.

    // First array must be a monotonically increasing function, start from redshift zero rather than redshift 2.0
    spline(z_Array, ComovingDistance_z, nPoints, 1.0e31, 1.0e31, z_ComovingDistance_2derivatives);
    spline(ComovingDistance_z, z_Array, nPoints, 1.0e31, 1.0e31, ComovingDistance_z_2derivatives);

    return 0;
}


float Integrand(float x){
    return pow(Om_v + Om_m*pow(1.+ x, 3.) + Om_r*pow(1.+ x, 4.) - (Om_tot -1.)*pow(1.+ x, 2.), -0.5);
}


// Returns comoving distance at redshift z in h^-1 Mpc. 
double interp_comovingDistance(double z){
    float InterimInterp_yVal;
    splint(z_Array, ComovingDistance_z, z_ComovingDistance_2derivatives, nPoints, (float) z, &InterimInterp_yVal);
    return InterimInterp_yVal;
}


double interp_inverseComovingDistance(double r){
    float InterimInterp_yVal;
    splint(ComovingDistance_z, z_Array, ComovingDistance_z_2derivatives, nPoints, (float) r, &InterimInterp_yVal);
    return InterimInterp_yVal;
}   // Returns z at comoving distance, [h^-1 Mpc]. 


double HubbleCnst(double z){
    // Units of h km s^-1 Mpc^-1   
    return 100.*pow(Om_v + Om_m*pow(1. + z, 3.), 0.5);
}
