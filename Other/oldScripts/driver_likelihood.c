#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include "header.h"
#include "cosmology_planck2015.h"
#include "comovDistRedshiftCalc.c"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "assignMemory.c"
#include "KaiserMultipoles.c"    
#include "KaiserLorentzMultipoles.c"
#include "FFTw.c"
#include "FFT_log.h"
#include "FFT_log.c"
#include "toymodel_pk_xi.c"
#include "matter_pk.c"
#include "VIPERS_window.c"
#include "VIPERS_kwindow.c"
#include "VIPERS_window_jointfield.c"
#include "MultipoleCovariance.c"
#include "MultipoleCovariance_eigenvecs.c"
#include "ChiSq_minimisation.c"
#include "ChiSq_eval.c"
#include "ChiSq_input.c"
#include "posteriors_1D.c"
#include "posteriors_2D.c"
#include "Metropolis_algorithm.c"
#include "Alcock_Paczynski.c"
#include "Ruiz.c"
#include "combining_clipped_fsig8.c"
#include "freeMemory.c"


int main(int argc, char **argv){
  sprintf(root_dir,      "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");

  AxisLimsArray[0][0]       =        0.0; 
  AxisLimsArray[1][0]       =      800.0;
 
  AxisLimsArray[0][1]       =        0.0;   
  AxisLimsArray[1][1]       =      800.0;                                                                                                    

  AxisLimsArray[0][2]       =        0.0;            // h^-1 Mpc
  AxisLimsArray[1][2]       =      800.0;
  
  fieldFlag                 =       atoi(argv[1]);
  appliedClippingThreshold  =       atof(argv[2]);
  lo_zlim                   =       atof(argv[3]);   // previously 0.6<z<0.9, 0.7<z<1.1
  hi_zlim                   =       atof(argv[4]);
  Jenkins_foldfactor        =       atof(argv[5]);   // Apply Jenkins folding to increase spatial resolution of mesh.
  ChiSq_kmax                =   0.1*atof(argv[6]);   // Easier to pass int in bash. 

  if(fieldFlag == 1){
    LowerRAlimit            =     30.175;            // Navgoya v6 + Samhain 
    UpperRAlimit            =     38.797;
    CentreRA                =     34.492;

    LowerDecLimit           =     -5.970;     
    UpperDecLimit           =     -4.171;     
    CentreDec               =     -5.091; 
    
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
  }
  
  else if(fieldFlag ==4){
    LowerRAlimit            =    330.046; // Really parent boundary limits. 
    UpperRAlimit            =    335.389;
    CentreRA                =    332.638;
  
    LowerDecLimit           =      0.862;     
    UpperDecLimit           =     2.3696;     
    CentreDec               =      1.583; 
    
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
  }
  
  W1area                    =     10.692;    // Nagoya v7 - overlapping Samhain; v7 and v6 are identical area. 
  W4area                    =      5.155;

  // data
  // W1area                 =      10.763;  // sq. degs. from Ben.
  // W4area                 =       5.155;
  
  TotalW1W4area             =      W1area + W4area; 
  
  fft_size                  =       256;     
  
  z_eff                     =      0.75;                  // set to weighted average of all galaxies in sample (Anderson et. al.)

  linearBias                =      1.53;                  // 1.32, appropriate for 0.6 to 0.9;
  velDispersion             =       3.0;                  // units of h^-1 Mpc rather than 300 km s^-1
  beta                      =     0.541;                  // 2dF measurement, beta = 0.43, beta = 0.542 for the cube.
  
  min_fsigma8               =      0.05;                  //  Priors on the model params.  
  max_fsigma8               =      0.80;

  min_velDisperse           =      0.00;                  // CHANGED FROM 0.00 13/02/2017
  max_velDisperse           =      6.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)

  min_bsigma8               =      0.05;                  // FOR GRANETT 2D POSTERIOR.
  max_bsigma8               =      1.00;                  // Previously 0.2 < b \sig_8 < 1.6
                                                          // 0.05 < b s8 < 1.0 (13/02/17)
  min_alpha_pad             =    0.9999;
  max_alpha_pad             =    1.0001;

  min_epsilon_pad           =   -0.0001;
  max_epsilon_pad           =    0.0001;
 
  min_A11Sq                 =      0.99;                  // distinct from bias due to spectral distortion. 
  max_A11Sq                 =      1.01;

  paramNumber               =       3.0;  // # of fitted params. -> dof in X^2. 

   Res                      =         2;  // Likelihood resolution [voxel number].
  dRes                      =       2.0;  // Previously 16: 13/02/17

   Res_ap                   =         1;  // Resoltuion in AP.
  dRes_ap                   =       1.0;

  ChiSq_kmin                =      0.02;
                 
  jenkins_fold_kjoin        =       0.4;  // k at which P(k) switches from unfolded to folded.     
          
  hiMultipoleOrder          =         2;  // Fit monopole (1) or mono + quad (2).

  chi_interval              =     16.00;  // Change to equal increments in volume?
  
  if(lo_zlim > 0.8){
    nz_smoothRadius         =      50.0;
  }

  else{
    nz_smoothRadius         =      100.;
  }
  
  fkpPk                     =    8000.0;            // FKP P_0 of interest [h^-1 Mpc]^3.

  modkMax                   =      1.00;            // Binned P(k) variables 
  muBinNumb                 =       100;
  kbin_no                   =        40;            

  CatalogNumber             =       306;            // Number of mocks. 
  
  clipping_smoothing_radius =       2.0;
   
  // analysis of VIPERS data or mock catalogues.
  data_mock_flag            =         0;
  
  printf("\n\n\nField flag: %d, d0: %.1lf, lo-z: %.1lf, hi-z: %.1lf, fold factor: %.1lf, chi sq. k_max: %.1lf", fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, Jenkins_foldfactor, ChiSq_kmax);  
  
  gsl_rng_env_setup();  // Random variable generation.

  gsl_ran_T = gsl_rng_default;
  
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);
  
  // Cosmology from cosmology_planck2015.h or cosmology_valueaddedmocks.h// Selection parameters.
  comovDistReshiftCalc();
  
  UniverseAge();
  
  linearGrowthRate();
  
  // inputHODPk();  //  Must match Cosmology in cosmology_planck2015.h, or cosmology_valueaddedmocks.h
  inputLinearPk();  //  This is NOT automatically ensured.               
     
  prep_fftw();
  
  prep_pkRegression(-2., log10(modkMax), kbin_no, 0);
  
  load_CovarianceMatrix(100, 1); // Number of mocks, starting mock. 
  
  prep_VIPERS_maskMultipoles();
  
  prep_VIPERS_jmaskMultipoles();

  W1_Spectro_V7_2_chiSq_minimisation();
    
  // ensemble_fsig8(2, 1, 8); //  ensemble_fsig8(int mockNumber, int start, int totalCats)  
  
  // calc_bestfit_fsig8(1, 0.4);  // Calc. best-fit fsig8 to all d0={1000, 10, 6, 4}.
  // calc_bestfit_fsig8(1, 0.6);
  // calc_bestfit_fsig8(1, 0.8);
  
  // metropolis_mcmc(mock_number); 
    
  printf("\n\n");

  return 0; 
}
