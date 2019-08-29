#define   KBIN_NO     40
#define   FOLDFACTOR 2.0       
 
// #include "/home/mjw/Aux_functions/header.h"
// #include "/home/mjw/Aux_functions/Aux_functions.c"

#include <stddef.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include "omp.h"

#include "header.h"
#include "header_pk.h"
// #include "cosmology_planck15.h" // #include "cosmology_rota16.h"
#include "struct_regress.h"
// #include "angular_limits.c"
// #include "AgeOftheUniverse.c"
// #include "linearGrowthRate.c"
// #include "chi_zcalc.c"
// #include "assignAcceptance.c"
// #include "volavgs.c"
#include "Initialise.c"
// #include "stefanoBasis.c"
// #include "max_gal.c"
#include "CoordinateCalc.c"
#include "MultipoleCalc.c"
// #include "randGen.c"
// #include "nbar.c"
// #include "invert_nbar.c"
// #include "load_mask.c"
// #include "get_zeff.c"
// #include "nbar_smooth.c"
// #include "fkp_weights.c"
// #include "rand_occupied.c"
// #include "calc_clippingweights.c"
#include "ngp.c"
#include "CloudInCell.c" 
#include "overdensity_calc.c"
#include "FFTw.c"
// #include "GaussianFilter.c"
#include "assign_pkmemory.c"
#include "assign_binnedpk_memory.c"
// #include "lock_files.c"
// #include "fftw_qlr.c"
// #include "fftw_qlk.c"


int main(int argc, char **argv){  
  (void) argc;                                              // escape compiler unused variable warning. 

  fold                      =                   0;

  thread                    =                   0; 
    
  smooth_radius             =                 2.0;          // smoothing kernel for clipping.

  max_gals                  =             3238855;
  
  //  outputdir                 = getenv("outputdir");  

  sprintf(filepath,      "/global/homes/m/mjwilson/UNIT/sim/HOD_Shadab/HOD_boxes/redshift0.9873/UNIT_DESI_Shadab_HOD_snap97_ELG_v0.txt");

  
  lopad                     =      50.0;                    // stefano_trans_z set such that lo_zlim boundary is lopad [h^-1 Mpc] from embedded vol boundary.
  fkpPk                     =    8000.0;                    // [h^-1 Mpc]^3.  Stefano: 4000 [h^-1 Mpc]^3.

  fft_size                  =       256;                    // Worker 46 works up to 1024. 
  
  logk_min                  =      -2.0;
  logk_max                  =   0.00000;                    // k = 1 hMpc^{-1} :  0.00000;  k = 3 hMpc^{-1} :  0.47712;  k = 4 hMpc^{-1} : 0.60206 
  
  //  start_walltime();
  
  //  printf_branch();
  
  fftw_init_threads();
  
  fftw_plan_with_nthreads(omp_get_max_threads());        // Maximum number of threads to be used; use all openmp threads available.  
  
  //  set_angularlimits(0, fieldFlag);                   // Assumes data cut to mock limits.
  
  Initialise();                                          // Initialise grid, FFT params and random generation.
  
  //  check_radialextent(loChi, hiChi, lopad);
  
  //  stefano_trans_z = -loChi + lopad;                  // Translate lower limit to near the edge of box.
  
  prep_CatalogueInput_500s(max_gals);                    // Max. number of gals of ALL mocks (& data) analysed simultaneously is `hard coded' (with some contingency).  
  
  //  prep_nbar();                                       // assign memory and zero e.g. bins of number of galaxies per chi. 
  
  prep_x2c();                                            // Memory for overdensity, smooth_overdensity and H_k; either double or fftw_complex.                                                                                                                                                                                                                                         
  prep_pkRegression();                                   // set k binning interval arrays.  
  
  // onemock_nbarcalc();                                 // rewrites FieldFlag and loopCount. Code must exit at this point.  
  
  // load_rands_radec(1.0);
  
  // get_zeff();

  // calc_volavg_fkpweights2();

  // delete_lockfile();
  
  // prep_clipping_calc();
  
  prep_r2c_modes(&flat,                    1.0);         // unfolded.
  prep_r2c_modes(&half,             FOLDFACTOR);         // one fold.
  prep_r2c_modes(&quart, FOLDFACTOR*FOLDFACTOR);         // two folds.

  regress set[3] = {flat, half, quart};
  // int     d0s[5] = {2, 4, 6, 10, 1000};

  // int mock_end   = atoi(argv[4]) + atoi(argv[5]) > CatalogNumber ? CatalogNumber : atoi(argv[4]) + atoi(argv[5]);
  
  // walltime("All prep. done");
    
  CatalogueInput_500s(max_gals);  // mocks 1 to 153 are independent. 
  
  // set_zcut();                  // applies z cut      
  // set_deccut();                // applies (problem) dec cut in mocks to data.
      
  // vollim_cutbyMB(-20.5);
    
  // StefanoBasis(Vipers_Num, ra, dec, xCoor, yCoor, zCoor);  // applied to both gals and rands.  (ra, dec, redshift) to (x, y, z) in Stefano's basis.
    
  // print_xy();
      
  // spline_nbar(0);           // new <n(z)> for each mock. arg 1: bool for smoothed + reflected 2-field avg., arg 2: 'truth' i.e. mock avg.
      
  // pt2nz = &interp_nz;       // Gaussian smoothed galaxy counts. 
  // pt2nz = &vollim_nz;q
    
  // accepted_gal();           // count number of accepted galaxies.
    
  // prep_inverseCumulative_nbar();
      
  // calc_clipping_weights(); 
    
  // walltime("Clipping weights done.");
    
  // print_metd0();
  
  // rand_newchi_newbasis();
    
  // assign_randbox();
  
  // fft_qellk(&flat);
        
  // calc_bare_fkpweights();   // fkp_weights in units of alpha. 
    
  // printf("\n\n");
  
  calc_overdensity(max_gals);
  
  PkCalc(&set[fold], 0);
  
  // walltime("Wall time at finish");
  
  // MPI_Finalize();
  
  printf("\n\n");
  
  exit(EXIT_SUCCESS);
}
