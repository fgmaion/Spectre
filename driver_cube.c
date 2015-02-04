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

#include "Scripts/comovDistRedshiftCalc.c"

// #include "Scripts/JenkinsRun.c"
#include "Scripts/GridParams.c"

#include "Scripts/assignMemory.c"

#include "Scripts/load_mask.c"

#include "Scripts/assignAcceptance.c"

#include "Scripts/CoordinateCalc.c"

#include "Scripts/overdensity_calc.c"
#include "Scripts/CloudInCell.c"
// #include "Scripts/BasisChange.c"
// #include "Scripts/CalcCellraDec.c"

// #include "Scripts/KaiserMultipoles.c"
// #include "Scripts/KaiserGaussMultipoles.c"
#include "Scripts/KaiserLorentzMultipoles.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"

#include "Scripts/nbar.c"
// #include "Scripts/MockAvgComovingDensity.c"

// #include "Scripts/AgeOftheUniverse.c"
// #include "Scripts/linearGrowthRate.c"
// #include "Scripts/growthfactor_derivative.c"

#include "Scripts/matter_pk.c"
// #include "Scripts/Clipped_zSpace.c"
#include "Scripts/toymodel_pk_xi.c"

// #include "Scripts/ArtificialWf.c"
// #include "Scripts/BootStrap.c"

// #include "Scripts/correlation_fns.c"

#include "Scripts/randGen.c"

#include "Scripts/FFT_log.h"
#include "Scripts/FFT_log.c"

/*#include "Scripts/cubature/cubature.h"*/
/*#include "Scripts/FFT_log_zeldovich.h"*/
/*#include "Scripts/FFT_log_zeldovich.c"*/

/*#include "Scripts/anisotropicGaussian_multipoles.c"*/

/*#include "Scripts/HOD_mock_theoryExp.c"*/

/*#include "Scripts/MultipoleCovariance.c"*/
/*#include "Scripts/MultipoleCovariance_eigenvecs.c"*/
/*#include "Scripts/ChiSq_minimisation.c"*/
/*#include "Scripts/posteriors_1D.c"*/
/*#include "Scripts/posteriors_2D.c"*/

/*#include "Scripts/MonteCarlo_SSPOC.c"*/
/*#include "Scripts/AngularSelectionCats.c"*/
/*#include "Scripts/SaundersDeproject.c"*/

#include "Scripts/smith_mjw.h"
#include "Scripts/smith_mjw.c"

#include "Scripts/libkdtree.h"
#include "Scripts/kdtree_xi_mom.c"
#include "Scripts/buildTree.c"
#include "Scripts/libkdtree.c"

/*#include "Scripts/mockGalaxyCats.c"*/

#include "Scripts/halomodel_pk.c"
#include "Scripts/VIPERS_window.c"

// #include "tinker.c"

#include "Scripts/hod_cube.c"

#include "Scripts/freeMemory.c"


int main(int argc, char **argv){
  sprintf(root_dir,      "/disk1/mjw/HOD_MockRun");

  sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_HOD_Mocks");

  // Embedding volume for P(k) measurement. Stefano basis. 
  AxisLimsArray[0][0]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][0]   =    1000.0;
 
  AxisLimsArray[0][1]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][1]   =    1000.0;

  AxisLimsArray[0][2]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]   =    1000.0;

  // Cell size, comoving distance, h^-1 Mpc. 
  fft_size              =       512;     

  // Binning interval for P(k).
  modkMax               =      1.00;
  muBinNumb             =       100;
  kbin_no               =        40;

  // Interval in k^2 for perp k binning of 2D P(k).
  perpkInterval         =      0.01;
 
  // Total number of HOD mocks. 
  CatalogNumber         =       306;
    
  // Clipping variables. 
  appliedClippingThreshold  =   1.0;    

  // Random variable generation.    
  gsl_rng_env_setup();

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

  // comovDistReshiftCalc();
  
  // VIPERS_SolidAngle = SolidAngleCalc(LowerDecLimit, UpperDecLimit, UpperRAlimit-LowerRAlimit);

  EvaluateGridParameters();
  
      
  
  sprintf(filepath, "%s/mocks_W1_v9.0_500/cube_-20.0.dat", vipersHOD_dir);
    
  CatalogueInput_Cube(filepath);
    
  calc_overdensityCube();

  cube_PkCalc();
  
  printf("\n\n");

  return 0; 
}
