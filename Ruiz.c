// Compute F_AP(z) and f_sig8(z) for Ruiz plot, given e.g. Planck cosmology. 
// F_AP(z)   = (1+z) H(z) D_A(z) / c
// f_sig8(z) = f(z)*sig8(0)*D_+
// D_A(z)    = (1+z)^-1. * R_0 S_k(r), assuming flat S_k(r) = r.

// Useful functions. 
// f_Om_545(double lna)
// HubbleCnst(double z) // Units of h km s^-1 Mpc^-1
// interp_comovingDistance(double z)  // Returns comoving distance at redshift z in h^-1 Mpc.
// linearGrowth_factor(double ln a)

// cosmology_planck2015.h
// double  Om_v      =    0.69;                 
// double  Om_r      =    0.00;                                                                                    
// double  Om_m      =    0.31;        // total matter, cdm + baryons.                                                        // double  Om_b      =   0.048;                                                                         
// double  h         =   0.673;                                                                        
// double  sigma_8   =    0.82;                                                                                     

int Ruiz_locus(){
  double lightspeed = 2.9979*pow(10., 5.);
  double z, ln_a, f, H, chi, Da, F_AP, Dplus, f_sig8;

  // sprintf(filepath, "%s/W1_Spectro_V7_2/Ruiz_locus_GR.dat", root_dir);

  // output = fopen(filepath, "w");

  // for(j=0; j<100; j++){
  // z     = 1.5*j/100.;

    z     = 1.05;

    ln_a  = -log(1. + z);

    f     = f_Om_545(ln_a);

    // Units of km s^-1 Mpc^-1  
    H     = h*HubbleCnst(z); 

    chi   = interp_comovingDistance(z);  // Returns comoving distance at redshift z in h^-1 Mpc
 
    Da    = pow(1.+z, -1.)*chi/h;        // Returns angular diameter distance at redshift z in Mpc  

    F_AP  = (1. + z)*Da*H/lightspeed;

    Dplus = linearGrowth_factor(ln_a);

    f_sig8 = f*sigma_8*Dplus;

    printf("\n\nFSIG8= %.4lf", f_sig8);

    // fprintf(output, "%e \t %e \t %e \n", z, f_sig8, F_AP);
    // }

  // fclose(output);

  return 0;
}
