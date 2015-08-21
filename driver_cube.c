// Stacpolly run. 
// #define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

#include <stdbool.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_expint.h>

// #include "omp.h"

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#define   AUXfn_funcs  "/disk1/mjw/Aux_functions/Aux_functions.c"
#include  AUXfn_funcs

#include "Scripts/header.h"

// #include "Scripts/cosmology_planck2015.h"
#include "Scripts/cosmology_valueaddedmocks.h"

#include "Scripts/comovDistRedshiftCalc.c"

#include "Scripts/Jenkins_fold.c"
#include "Scripts/GridParams.c"

#include "Scripts/assignMemory.c"

#include "Scripts/load_mask.c"

#include "Scripts/assignAcceptance.c"

#include "Scripts/CoordinateCalc.c"
#include "Scripts/CoordinateCalcCube.c"

#include "Scripts/overdensity_calc.c"
#include "Scripts/CloudInCell.c"
// #include "Scripts/BasisChange.c"
// #include "Scripts/CalcCellraDec.c"

#include "Scripts/KaiserMultipoles.c"
// #include "Scripts/KaiserGaussMultipoles.c"
#include "Scripts/KaiserLorentzMultipoles.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw.c"

#include "Scripts/nbar.c"
#include "Scripts/nbar_smooth.c"
// #include "Scripts/MockAvgComovingDensity.c"

#include "Scripts/AgeOftheUniverse.c"
#include "Scripts/linearGrowthRate.c"
// #include "Scripts/growthfactor_derivative.c"

#include "Scripts/smith_mjw.h"
#include "Scripts/smith_mjw.c"

#include "Scripts/toymodel_pk_xi.c"
#include "Scripts/matter_pk.c"
#include "Scripts/Clipped_zSpace.c"

// #include "Scripts/ArtificialWf.c"
// #include "Scripts/BootStrap.c"

// #include "Scripts/correlation_fns.c"

#include "Scripts/randGen.c"

#include "Scripts/FFT_log.h"

#include "Scripts/clipped_lnnormal.c"

#include "Scripts/FFT_log.c"
#include "Scripts/FFTw_3Dwf.c"

/*#include "Scripts/cubature/cubature.h"*/
/*#include "Scripts/FFT_log_zeldovich.h"*/
/*#include "Scripts/FFT_log_zeldovich.c"*/

/*#include "Scripts/anisotropicGaussian_multipoles.c"*/

/*#include "Scripts/HOD_mock_theoryExp.c"*/

#include "Scripts/MultipoleCovariance.c"
#include "Scripts/MultipoleCovariance_eigenvecs.c"
#include "Scripts/ChiSq_minimisation.c"
#include "Scripts/posteriors_1D.c"
// #include "Scripts/posteriors_2D.c"

/*#include "Scripts/MonteCarlo_SSPOC.c"*/
/*#include "Scripts/AngularSelectionCats.c"*/
/*#include "Scripts/SaundersDeproject.c"*/

#include "Scripts/libkdtree.h"
#include "Scripts/kdtree_xi_mom.c"
#include "Scripts/buildTree.c"
#include "Scripts/libkdtree.c"

/*#include "Scripts/mockGalaxyCats.c"*/

#include "Scripts/halomodel.c"
#include "Scripts/VIPERS_window.c"

#include "Scripts/libkdpoly.h"
#include "Scripts/libkdpoly.c"
#include "Scripts/spec_weights.c"
#include "Scripts/fkp.c"

// #include "tinker.c"

#include "Scripts/freeMemory.c"

int main(int argc, char **argv){
  // char* s = getenv("ROOTDIR");
  sprintf(root_dir,      "/disk1/mjw/HOD_MockRun");

  sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_HOD_Mocks");
  // sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_ValueAddedHOD");

  // MPI_Init(&argc,&argv);
  // MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &process_number);
 
  // With orientation of the -z Cartesian axis to the line of sight. 
  // lower_xlimit & upper_xlimit
  // AxisLimsArray[0][0]   =    1550.0;                                                  // h^-1 Mpc
  // AxisLimsArray[1][0]   =    2180.0;                                                  // h^-1 Mpc

  // lower_ylimit & upper_ylimit
  // AxisLimsArray[0][1]   =    -170.0;                                                  // h^-1 Mpc
  // AxisLimsArray[1][1]   =     170.0;                                                  // h^-1 Mpc

  // lower_zlimit & upper_zlimit
  // AxisLimsArray[0][2]   =     -75.0;                                                  // h^-1 Mpc
  // AxisLimsArray[1][2]   =     -10.0;                                                  // h^-1 Mpc

  // Embedding volume for P(k) measurement. Stefano basis. 
  AxisLimsArray[0][0]   =       0.0; 
  AxisLimsArray[1][0]   =    1000.0;
 
  AxisLimsArray[0][1]   =       0.0;   
  AxisLimsArray[1][1]   =    1000.0;                                                     // h^-1 Mpc                                                 

  AxisLimsArray[0][2]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]   =    1000.0;

  // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
  // stefano_trans_x       =      +50.; 
  // stefano_trans_y       =     +250.; 
  // stefano_trans_z       =    -1425.;
  
  // stefano_trans_x       =     +250.;
  // stefano_trans_y       =     +250.;
  // stefano_trans_z       =    -1700.;
  
  /*
  // W1 catalogue. new 500s mocks. parent
  LowerRAlimit          =       30.1; 
  UpperRAlimit          =       38.8;
  CentreRA              =      34.45;

  LowerDecLimit         =      -5.95;     
  UpperDecLimit         =      -4.15;     
  CentreDec             =      -5.05;  
  
  W1area                =      15.66;      
  */
  
  // W1 catalogue. new 500s mocks. Nagoya v4 spectroscopic mask. Nagoya v5 is identical except for a single pointing removed, to be reobserved.
  LowerRAlimit          =      30.17; 
  UpperRAlimit          =       38.8;
  CentreRA              =     34.487;
  
  LowerDecLimit         =      -5.38;     
  UpperDecLimit         =      -4.17;     
  CentreDec             =      -4.77;     
  
  // calculated in spec_weights.c
  // W1area             =      7.471;
  
  // official release info, maybe photo mask is included. v4 Nagoya.
  W1area                =      7.017;  // sq. degs.
  W4area                =      5.150;
  
  TotalW1W4area         =      W1area + W4area; 
  
  // Cell size, comoving distance, h^-1 Mpc. 
  fft_size              =       256;     

  // Selection parameters.
  lo_MBlim              =     -90.5;                  // -20.5 < M_B < -19.5
  hi_MBlim              =     -00.0;
  
  lo_zlim               =      0.60;                  // previously 0.6<z<0.9, 0.7<z<1.1 
  hi_zlim               =      0.90;
  z_eff                 =      0.75;                  // set to volume avg. redshift of the survey? prior to FKP weights.

  linearBias            =      1.53;                  // 1.32, appropriate for 0.6 to 0.9;
  velDispersion         =       3.0;                  // units of h^-1 Mpc rather than 300 km s^-1
  beta                  =     0.541;                  // 2dF measurement, beta = 0.43, beta = 0.542 for the cube. 

  // Priors on the model params.
  min_fsigma8           =      0.00;
  max_fsigma8           =      0.80;

  min_velDisperse       =      0.00;
  max_velDisperse       =      6.00;

  // 4 param likelihood required. change this. 
  min_bsigma8           =       0.8;
  max_bsigma8           =       1.6;
 
  // distinguished from linear bias by spectral distortion.  
  min_A11Sq             =      0.99;
  max_A11Sq             =      1.01;

  // Number of fitted parameters, defines degrees of freedom in chi sq. expectation. 
  paramNumber           =       3.0;

  // Resolution of the Likelihood evaluation [voxel number].
   Res                   =       16;
  dRes                  =      16.0;

  ChiSq_kmin            =      0.04;
  ChiSq_kmax            =      0.20; 
     
  // Fit solely the monopole (1) or both mono and Quad (2).
  hiMultipoleOrder      =         2;

  // Comoving number density, n(z), measurement. Change to equal increments in volume?
  chi_interval          =     16.00;

  // Apply Jenkins folding to increase spatial resolution of mesh. 
  Jenkins_foldfactor    =       1.0;

  // FKP P(k) of interest.
  fkpPk                 =    3000.0;            // [h^-1 Mpc]^3.
  meanSampling          =       0.4;

  // Binning interval for P(k).
  modkMax               =       1.0;
  muBinNumb             =       100;
  kbin_no               =        40;

  // Total number of HOD mocks. 
  CatalogNumber         =       306;

  //  Apodise the window fn. to supress ringing of the window ("Gibb's phenomenon").
  GibbsSkinDepth            =         5.0;
    
  // Remember to smooth.
  appliedClippingThreshold  =         3.0;    
  
  clipping_smoothing_radius =         2.0;

  GaussianFilter_radius     =         0.0;

  // subsample HOD catalogue 
  depletion_factor          =        1.25;

  // nfw halo generation
  nfw_conc                  =         1.0; 

  // Correlation fn's, logarithmic binning in r. 
  zerolog  =                 log10(0.001);
  maxlog   =                   log10(2.0);      // hiRes: 20., lowRes: 2000.
  logbinsz =                 log10( 1.01);      // previously 1.4, must be >1.0 otherwise log gives 0. or -ve.
  
  nlogbins =  (int) ceil((maxlog - zerolog)/logbinsz);

  // linear binning in mu. 
  zerolin  =                         0.00;
  maxlin   =                         1.00;
  linbinsz =                         0.05;

  nlinbins =  (int) ceil((maxlin - zerolin)/linbinsz);

  // Random variable generation.    
  gsl_rng_env_setup();

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

  
  comovDistReshiftCalc();
  
  Jenkins_foldEmbeddingVol();
  
  EvaluateGridParameters();
  
  
  cube_rands();
  
  Jenkins_foldRand();
    

  CoordinateCalcCube("/disk1/mjw/HOD_MockRun/Data/HODCube/cube_gal_-20.0.dat");

  assignAcceptanceCube();
  
  Jenkins_foldCat();
  
  // clipped lognormal model investigation.
  // Gauss_varlognormal();
  
  
  calc_clippingweights();
  
  load_clippingweights();
  
  calcCube_overdensity();
  
  gridup_randsmask();
  
  // Calculate the galaxy contribution to the shotnoise, galaxies
  // are weighted by gal_weights.  assumes weights are density independent.  
  // weighted_shotnoise();
  
  cube_PkCalc();
  
  /*
  // AgeOftheUniverse.c
  UniverseAge();
  
  // linearGrowthRate.c
  linearGrowthRate();
  
  // needs checked.
  // halofit(z_eff);
  
  //  WATCH OUT: Loads Value added mocks cosmology P(k), change back to Planck 2015. // 
  inputHODPk();
  
  // clipped_lnnormal_main();
  
  clipped_lnnorm_pkCalc();
  */
  printf("\n\n");
  
  return 0; 
}
