#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include "header.h"
#include "cosmology_planck15.h"
#include "header_pk.h"
#include "header_W2.h"
#include "libkdtree.h"
#include "fitsio.h"
#include "angular_limits.c"
#include "chi_zcalc.c"
#include "invert_nbar.c"
#include "load_maskQmultipoles.c"
#include "assignAcceptance.c"
#include "volavgs.c"
#include "max_gal.c"
#include "CoordinateCalc.c"
#include "lock_files.c"
#include "nbar.c"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "kdtree_xi_mom.c"
#include "buildTree.c"
#include "libkdtree.c"
#include "randGen.c"
#include "fitsio.c"
#include "kdtree_qsnapshot.c"


int init_gsl_randgen(){
  gsl_rng_env_setup();  // random generation setup.

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);  // returns a pointer to a newly-created instance of a random number generator of type T

  return 0;
}

int initi_dist_z(){
  chi_zcalc(); // Cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h

  printf("\n\nRedshift limits, lower bound: %.4lf \t %.4lf h^-1 Mpc, \n\t\t upper bound: %.4lf \t %.4lf h^-1 Mpc", lo_zlim, loChi, hi_zlim, hiChi);

  return 0;
}

int set_file_paths(char survey[], int count_res){
  if(count_res < 3){  
    sprintf(surveyType, survey, fieldFlag, fkpPk, lo_zlim, hi_zlim, thread);
  }

  else{
    sprintf(surveyType, survey, fkpPk, lo_zlim, hi_zlim, thread);
  }
  
  return 0;
}

int set_outputfiles(int count_res){
  switch(count_res){
  case 0:
    set_file_paths("Ql_W%d_Nag_v7_specweight_nbar_Pfkp_%.0lf_%.1f_%.1f_thread_%d_hihiRes_hex", 0);
    break;

  case 1:
    set_file_paths("Ql_W%d_Nag_v7_specweight_nbar_Pfkp_%.0lf_%.1f_%.1f_thread_%d_hiRes_hex",   1);
    break;

  case 2:
    set_file_paths("Ql_W%d_Nag_v7_specweight_nbar_Pfkp_%.0lf_%.1f_%.1f_thread_%d_loRes_hex",   2);
    break;

  case 3:
    set_file_paths("Ql_W1W4_Nag_v7_specweight_nbar_Pfkp_%.0f_%.1f_%.1f_thread_%d_hihiRes_hex", 3);
    break;

  case 4:
    set_file_paths("Ql_W1W4_Nag_v7_specweight_nbar_Pfkp_%.0f_%.1f_%.1f_thread_%d_hiRes_hex",   4);
    break;

  case 5:
    set_file_paths("Ql_W1W4_Nag_v7_specweight_nbar_Pfkp_%.0f_%.1f_%.1f_thread_%d_loRes_hex",   5);
    break;
  }

  return 0;
}

int main(int argc, char **argv){
  int              count_res;
  double       sampling_frac;

  double      dilution = 0.0001;

  // Randoms in cats. are W1: 2.5196772 x 10^7; W4: 1.2147352 x 10^7; W1+W4: 6.2245045 x 10^7.
  
  double       max_logs[6]  = {log10(2.0), log10(20.0), log10(2000.0), log10(2.0), log10(20.0), log10(4000.0)}; 
  double sampling_fracs[6]  = {1.00, 0.30, 0.01, 1.00, 0.20, 0.01}; // 1.00, 0.10, 0.01, 1.000, 0.004, 0.005
  
  thread                    =                                   1;
  outputdir                 =                 getenv("outputdir");
  
  sprintf(root_dir,                      "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");
    
  fieldFlag                 =                       atoi(argv[1]);
  lo_zlim                   =                       atof(argv[2]); // previously 0.6<z<0.9, 0.7<z<1.1
  hi_zlim                   =                       atof(argv[3]);
  count_res                 =                       atoi(argv[4]);

  maxlog                    =                 max_logs[count_res];
  sampling_frac             =  dilution*sampling_fracs[count_res];

  // Comoving number density, n(z), measurement. Change to equal increments in (effective) volume?
  chi_interval              =     16.00;
  
  // Change from 150 to match Stefano's 100.
  nz_smoothRadius           =     100.0;

  // FKP p(k) of interest;
  fkpPk                     =    4000.0;            // [h^-1 Mpc]^3.

  CatalogNumber             =       153;
      
  // analysis of VIPERS data or mock catalogues.
  data_mock_flag            =         0;
  
  // Correlation fn's, logarithmic binning in r. 
  zerolog  =               log10(0.001);
  logbinsz =               log10(1.050);    // previously 1.01 and befor that 1.4. must be > 1.0 otherwise log gives 0. or -ve.
  
  // linear binning in mu. 
  zerolin  =                      0.000;
  maxlin   =                      1.000;
  linbinsz =                      0.025;    // previously 0.05


  start_walltime();
  
  printf_branch();
  
  init_gsl_randgen();

  initi_dist_z();
  
  UniverseAge();

  set_angularlimits(0, fieldFlag); // assumes data/randoms cut to mock limits.
  
  prep_nbar();
  
  spline_nbar(1);                  // 306 mock avg; load nbar jointly estimated from W1 and W4, with 150 h^-1 Mpc smoothing and reflection.  

  pt2nz = &interp_nz;              // Gaussian smoothed galaxy counts.

  prep_inverseCumulative_nbar();
  
  set_outputfiles(count_res);
  
  load_maskfits(sampling_frac, count_res);
  // load_homogeneous_rands_window(sampling_frac, count_res); // load randoms with sampling sampling_frac.
  
  delete_lockfile();
  
  assignMemory_xi();

  randWindow_pairCount();  
  
  walltime("Wall time at finish");
  
  printf("\n\n");

  exit(EXIT_SUCCESS);
  
  return 0;
}
