// Stacpolly run. 
#include <stdbool.h>
#include <time.h>

#include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
// #include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_multimin.h>

// #include "omp.h"

//#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
//#include  AUXfn_header

// #define   AUXfn_funcs  "/disk1/mjw/Aux_functions/Aux_functions.c"
// #include  AUXfn_funcs

#define  AUXfn_header "/home/mjw/Aux_functions/header.h"
#include AUXfn_header

#define  AUXfn_funcs  "/home/mjw/Aux_functions/Aux_functions.c"
#include AUXfn_funcs

#include "header.h"

#include "cosmology_planck2015.h"
// #include "Scripts/cosmology_valueaddedmocks.h"

#include "comovDistRedshiftCalc.c"

#include "Jenkins_fold.c"
#include "GridParams.c"

#include "assignMemory.c"

#include "load_maskQmultipoles.c"

#include "assignAcceptance.c"

#include "CoordinateCalc.c"
#include "CoordinateCalcCube.c"

#include "overdensity_calc.c"
#include "CloudInCell.c"
// #include "Scripts/BasisChange.c"
// #include "Scripts/CalcCellraDec.c"

#include "KaiserMultipoles.c"
// #include "Scripts/KaiserGaussMultipoles.c"
#include "KaiserLorentzMultipoles.c"

#include "qSortCompare.c"
#include "FFTw.c"

#include "nbar.c"
#include "nbar_smooth.c"

// #include "Scripts/MockAvgComovingDensity.c"
#include "nbar_fit.c"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
// #include "Scripts/growthfactor_derivative.c"

// #include "Scripts/smith_mjw.h"
// #include "Scripts/smith_mjw.c"

#include "toymodel_pk_xi.c"
#include "matter_pk.c"
#include "Clipped_zSpace.c"

#include "ArtificialWf.c"
// #include "Scripts/BootStrap.c"

// #include "Scripts/correlation_fns.c"

#include "fkp_weights.c"
#include "clipping_weights.c"

#include "randGen.c"

#include "FFT_log.h"

#include "clipped_lnnormal.c"

#include "FFT_log.c"
// #include "Scripts/FFTw_3Dwf.c"

/*#include "Scripts/cubature/cubature.h"*/
/*#include "Scripts/FFT_log_zeldovich.h"*/
/*#include "Scripts/FFT_log_zeldovich.c"*/

/*#include "Scripts/anisotropicGaussian_multipoles.c"*/

/*#include "Scripts/HOD_mock_theoryExp.c"*/

#include "MultipoleCovariance.c"
#include "MultipoleCovariance_eigenvecs.c"
#include "ChiSq_minimisation.c"
#include "posteriors_1D.c"
// #include "Scripts/posteriors_2D.c"

/*#include "Scripts/MonteCarlo_SSPOC.c"*/
/*#include "Scripts/AngularSelectionCats.c"*/
/*#include "Scripts/SaundersDeproject.c"*/

#include "libkdtree.h"
#include "kdtree_xi_mom.c"
#include "buildTree.c"
#include "libkdtree.c"

/*#include "Scripts/mockGalaxyCats.c"*/

// #include "Scripts/halomodel.c"
#include "VIPERS_window.c"
#include "VIPERS_window_jointfield.c"
// #include "Scripts/libkdpoly.h"
// #include "Scripts/libkdpoly.c"
// #include "Scripts/spec_weights.c"
// #include "Scripts/fkp.c"
// #include "Scripts/Bailey.c"
#include "Alcock_Paczynski.c"
#include "super_vipers.c"

#include "paircount_mask.c"
// # include "Scripts/Ruiz.c"

// #include "tinker.c"

#include "freeMemory.c"


int main(int argc, char **argv){
  // char* s = getenv("ROOTDIR");
  sprintf(root_dir,      "/home/mjw/HOD_MockRun");

  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");
  // sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_HOD_Mocks");
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
  AxisLimsArray[0][0]       =        0.0; 
  AxisLimsArray[1][0]       =      800.0;
 
  // ALTERED 800 -> 400 
  AxisLimsArray[0][1]       =        0.0;   
  AxisLimsArray[1][1]       =      800.0;                                                     // h^-1 Mpc                                                 

  // ALTERED 800 -> 400
  AxisLimsArray[0][2]       =        0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]       =      800.0;
  
  fieldFlag                 =     atoi(argv[1]);
  
  if(fieldFlag == 1){
    // W1 catalogue. Nagoya v6 spectroscopic mask (& Samhain).
    LowerRAlimit            =     30.18;   // 30.1893; 
    UpperRAlimit            =     38.90;   // 38.8022;
    CentreRA                =     34.3213; // 34.3213;

    LowerDecLimit           =     -5.99;   // -5.9801;     
    UpperDecLimit           =     -4.16;   // -4.1715;     
    CentreDec               =     -5.1188; 
    
    // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
  }
  
  else if(fieldFlag ==4){
    // W4 catalogue. Nagoya v6 spectroscopic mask (& Samhain). parent boundary limits. 
    LowerRAlimit            =    330.0452; 
    UpperRAlimit            =    335.3890;
    CentreRA                =    332.7049;
  
    LowerDecLimit           =      0.8621;     
    UpperDecLimit           =     2.36950;     
    CentreDec               =      1.5549; 
    
    // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
  }
  
  // official release info, (Nagoya v7 - overlapping Samhain area). With regards angular coverage v7 and v6 are identical. 
  W1area                    =      10.763;  // sq. degs. from Ben. 
  W4area                    =       5.155;

  // Required for <n(z)> calculation.
  TotalW1W4area             =      W1area + W4area; 
  
  fft_size                  =        256;     

  // Selection parameters.
  lo_MBlim                  =     -90.5;                  // -20.5 < M_B < -19.5
  hi_MBlim                  =     -00.0;
  
  lo_zlim                   =      atof(argv[3]);                  // previously 0.6<z<0.9, 0.7<z<1.1 
  hi_zlim                   =      atof(argv[4]);
  z_eff                     =      0.75;                  // set to weighted average of all galaxies in sample (Anderson et. al.)

  // analysis of VIPERS data or mock catalogues.
  data_mock_flag            =         0;
  
  // Random variable generation.    
  gsl_rng_env_setup();

  gsl_ran_T = gsl_rng_default;
  
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);
  
  
  //-- cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h --//
  comovDistReshiftCalc();
  
  // Jenkins_foldEmbeddingVol();
  
  EvaluateGridParameters();
  
  // AgeOftheUniverse.c
  UniverseAge();
  
  // linearGrowthRate.c
  linearGrowthRate();
  
  randoms_maskGen();
  
  // assigns memory for overdensity grid.
  // prep_grid();
  
  // prep_mask();
  
  // choice of mock for loading nbar estimate.                                                                                                                                                                    
    
  printf("\n\n");

  return 0;
}
