#define   KBIN_NO     40  // 12          
#define   FOLDFACTOR 2.0       
 
#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include "header.h"
#include "header_pk.h"
#include "cosmology_planck15.h" // #include "cosmology_rota16.h"
#include "struct_regress.h"
#include "angular_limits.c"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "chi_zcalc.c"
#include "assignAcceptance.c"
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
#include "GaussianFilter.c"
#include "assign_pkmemory.c"
#include "assign_binnedpk_memory.c"
#include "lock_files.c"
// #include "fftw_qlr.c"
// #include "fftw_qlk.c"


int main(int argc, char **argv){  
  (void) argc; // escape compiler unused variable warning. 

  thread                    =                   1; 
  
  fieldFlag                 =       atoi(argv[1]);
  lo_zlim                   =       atof(argv[2]);          // previously 0.6<z<0.9, 0.9<z<1.2
  hi_zlim                   =       atof(argv[3]);

  smooth_radius             =                 2.0;

  outputdir                 = getenv("outputdir");  

  sprintf(root_dir,      "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W%d", fieldFlag);

  // stefano_trans_z set by mean redshift in chi_zcalc.c
  stefano_trans_x           =      +100.;
  stefano_trans_y           =      +300.; 
  
  chi_interval              =     16.00;                    // Comoving number density, n(z), measurement.
  
  if(lo_zlim > 0.8){                                        // Change for new 0.6 < z < 0.8 and 0.8 < z < 1.0 limits
    nz_smoothRadius         =      50.0;
  }

  else{
    nz_smoothRadius         =      100.; 
  }

  lopad                     =      50.0;                    // stefano_trans_z set such that lo_zlim boundary is lopad [h^-1 Mpc] from embeeded volume boundary.
  fkpPk                     =    8000.0;                    // [h^-1 Mpc]^3.  Stefano: 4000 [h^-1 Mpc]^3.

  fft_size                  =       512;                    // Worker 46 works up to 1024. 
  
  logk_min                  =      -2.0;
  logk_max                  =   0.60206;                    // k = 4 hMpc^{-1}.
  
  CatalogNumber             =       153;                    // Total number of (independent) HOD mocks.

  
  
  start_walltime();
  
  printf_branch();
  
  fftw_init_threads();
  
  fftw_plan_with_nthreads(omp_get_max_threads());        // Maximum number of threads to be used; use all openmp threads available.  

  set_angularlimits(0, fieldFlag);                       // Cut data to mock limits.
  
  Initialise();                                          // Initialise grid, fft params and random generation.
  
  check_radialextent(loChi, hiChi, lopad);
  
  stefano_trans_z = -loChi + lopad;                      // Translate lower limit to close to edge of box.
  
  prep_x2c();                                            // Memory for overdensity, smooth_overdensity and H_k; either double or fftw_complex.
  
  prep_pkRegression();                                   // set k binning interval arrays. 
  
  prep_CatalogueInput_500s();                            // Max. number of gals of ALL mocks (& data) analysed simultaneously is 'hard coded'.  
  
  prep_nbar();                                           // assign memory and zero e.g. bins of number of galaxies per chi. 

  load_rands_radec(1.0);
  
  delete_lockfile();
  
  prep_clipping_calc();
  
  prep_r2c_modes(&flat,                    1.0); // unfolded.
  prep_r2c_modes(&half,             FOLDFACTOR); // one fold.
  prep_r2c_modes(&quart, FOLDFACTOR*FOLDFACTOR); // two folds.
  
  regress set[3] = {flat, half, quart};
  int     d0s[5] = {2, 4, 6, 10, 1000};

  int mock_end   = atoi(argv[4]) + atoi(argv[5]) > CatalogNumber ? CatalogNumber : atoi(argv[4]) + atoi(argv[5]);
  
  walltime("All prep. done");

  for(data_mock_flag = 0; data_mock_flag < 2; data_mock_flag++){ // analysis of VIPERS data and mock catalogues.
    trash_nbarshot_file(atoi(argv[4]));

    for(loopCount=atoi(argv[4]); loopCount <= mock_end; loopCount++){            
      sprintf(filepath, "%s/mock_%03d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, loopCount);
      
      if(data_mock_flag == 1){
        DataInput();            // includes ESR input.
      }
      
      else{
        CatalogueInput_500s();  // mocks 1 to 153 are independent. 
      }

      set_zcut();               // applies z cut      
      set_deccut();             // applies (problem) dec cut in mocks to data.
    
      // vollim_cutbyMB(-20.5);
    
      StefanoBasis(Vipers_Num, ra, dec, xCoor, yCoor, zCoor);  // applied to both gals and rands.  (ra, dec, z) to (x, y, z) in Stefano's basis.
    
      // print_xy();

      spline_nbar(0);           // new <n(z)> for each mock. arg 1: bool for smoothed + reflected 2-field avg., arg 2: 'truth' i.e. mock avg.
      
      pt2nz = &interp_nz;       // Gaussian smoothed galaxy counts. 
      // pt2nz = &vollim_nz;
    
      accepted_gal();           // count number of accepted galaxies.
    
      prep_inverseCumulative_nbar();
      
      calc_clipping_weights(); 
    
      walltime("Clipping weights done.");
    
      // print_metd0();
    
      rand_newchi_newbasis();
    
      // assign_randbox();
    
      // fft_qellk(&flat);
        
      calc_bare_fkpweights();   // fkp_weights in units of alpha. 
    
      printf("\n\n");
    
      for(int m=4; m<5; m++){
        d0 = d0s[m];
      
        set_clipping_weights(); // unity weights for d0=1000, else load. 
        
        alpha_calc();
            
        for(fold=0; fold<3; fold++){
          calc_overdensity();
        
          PkCalc(&set[fold], atoi(argv[4]));
        }
      }
    
      if(data_mock_flag == 1) break;
    }
  }
  
  walltime("Wall time at finish");
  
  // MPI_Finalize();
  
  printf("\n\n");

  exit(EXIT_SUCCESS);
}
