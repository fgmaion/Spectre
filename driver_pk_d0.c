#define  KBIN_NO 40          // Variables deciding memory allocation. 

#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include <sys/types.h>

#include "header.h"
#include "header_pk.h"
#include "cosmology_planck2015.h"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "comovDistRedshiftCalc.c"
#include "assignAcceptance.c"
#include "struct_regress.c"
#include "volavgs.c"
#include "Initialise.c"
#include "stefanoBasis.c"
#include "max_gal.c"
#include "CoordinateCalc.c"
#include "MultipoleCalc.c"
#include "randGen.c"
#include "nbar.c"
#include "invert_nbar.c"
#include "load_mask.c"
// #include "nbar_smooth.c"
#include "fkp_weights.c"
#include "rand_occupied.c"
#include "clipping_weights.c"
#include "ngp.c"
#include "CloudInCell.c" 
#include "overdensity_calc.c"
#include "FFTw.c"
//#include "old_FFTw.c"
#include "GaussianFilter.c"
#include "assign_pkmemory.c"


int main(int argc, char **argv){  
  thread                    =                   0; 

  data_mock_flag            =                   0;          // analysis of VIPERS data or mock catalogues.       
  
  fieldFlag                 =       atoi(argv[1]);
  d0                        =       atof(argv[2]);
  lo_zlim                   =       atof(argv[3]);          // previously 0.6<z<0.9, 0.9<z<1.2
  hi_zlim                   =       atof(argv[4]);

  smooth_radius             =                 2.0;

  sprintf(root_dir,      "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W%d", fieldFlag);
  
  if(fieldFlag == 1){
    LowerRAlimit            =     30.175;                 // W1 catalogue. Nagoya v6 spectroscopic mask (& Samhain).   
    UpperRAlimit            =     38.797;
    CentreRA                =     34.492;                 // Stefano:  CentreRA                =    34.4519;

    LowerDecLimit           =     -5.970;     
    UpperDecLimit           =     -4.171;     
    CentreDec               =     -5.091;                 // Stefano:  CentreDec               =     -5.07;
  }
    
  else if(fieldFlag == 4){
    LowerRAlimit            =    330.046;                 // W4 catalogue. Nagoya v6 & Samhain mask. parent boundary limits.     
    UpperRAlimit            =    335.389;
    CentreRA                =    332.638;
  
    LowerDecLimit           =      0.862;     
    UpperDecLimit           =     2.3696;     
    CentreDec               =      1.583; 
  }
  
  W1area                    =     10.692;                 // (Nagoya v7 - overlapping Samhain area). Coverage of v7 & v6 identical.       
  W4area                    =      5.155;                 // dec cut at -5.97 in the mocks.

  TotalW1W4area             = W1area + W4area;            // Required for <n(z)> calculation.
  
  stefano_trans_x           =      +100.;
  stefano_trans_y           =      +300.; 
  stefano_trans_z           =     -1720.;                 // Changed from 1500. on 06/03/2017; as dealing with 0.8 < z < 1.0        

  if(1.0 < hi_zlim){             //%% Changed from 0.8 %%
    stefano_trans_z -= 600.;     // Previously 600. Mpc for 0.9 < z < 1.2;
  }
  
  chi_interval              =     16.00;                 // Comoving number density, n(z), measurement.
  
  if(lo_zlim > 0.8){                                    // Change for new 0.6 < z < 0.8 and 0.8 < z < 1.0 limits
    nz_smoothRadius         =      50.0;
  }

  else{
    nz_smoothRadius         =      100.; 
  }

  fkpPk                     =    8000.0;                 // [h^-1 Mpc]^3.
  fft_size                  =       512;                 // Worker 46 works up to 1024. 
  
  logk_min                  =      -2.0;
  logk_max                  =   0.60206;                 // k = 4 hMpc^{-1}.
  
  CatalogNumber             =       306;                 // Total number of HOD mocks.

  start_walltime();
  
  fftw_init_threads();
  
  fftw_plan_with_nthreads(omp_get_max_threads());        // Maximum number of threads to be used; use all openmp threads available.  
  
  Initialise();                                          // Initialise grid, fft params and random generation.

  prep_x2c();                                            // Memory for overdensity, smooth_overdensity and H_k; either double or fftw_complex.
  
  prep_pkRegression();                                   

  prep_CatalogueInput_500s();                            // Requires max. number of gals of ALL mocks analysed simultaneously to be hard coded in.  
  
  prep_nbar();
  
  load_rands_radec(1.0);

  // prep_clipping_calc();
  
  prep_r2c_modes(&unit, 1.0); // unfolded.
  // prep_r2c_modes(&half, 2.0); // one fold.

  struct regress set[2] = {unit, half};
  
  walltime("All prep. done");
  
  for(loopCount=1; loopCount<2; loopCount++){            
    sprintf(filepath, "%s/mock_%03d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, loopCount);

    CatalogueInput_500s(); // mocks 1 to 153 are independent. 
    
    assignAcceptance();  
    
    spline_nbar(0);  // new <n(z)> for each mock. arg 1: bool for smoothed + reflected 2-field avg., arg 2: 'truth' i.e. mock avg.

    // set_clipping_weights(); // basis without rotation. 
    
    StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);  // applied to both gals and rands.  (ra, dec, z) to (x, y, z) in Stefano's basis.
    
    rand_newchi_newbasis();

    // loop over thresholds here. 
    alpha_calc();
    
    calc_fkpweights();  // normalisation of FKP weights set by random catalogue.

    for(fold=0; fold<1; fold++){
      calc_overdensity();
    
      // PkCalc(&set[fold]);

      PkCalc(&unit);
    }
  }
  
  walltime("Wall time at finish");

  printf("\n\n");
  
  return 0; 
}
