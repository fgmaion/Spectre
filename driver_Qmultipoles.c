
#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include "header.h"
#include "cosmology_planck15.h"
#include "header_pk.h"
#include "header_W2.h"
//#include "header_chi2.h"
#include "libkdtree.h"
#include "comovDistRedshiftCalc.c"
#include "invert_nbar.c"
#include "load_maskQmultipoles.c"
#include "assignAcceptance.c"
#include "volavgs.c"
#include "max_gal.c"
#include "CoordinateCalc.c"
#include "nbar.c"
#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "kdtree_xi_mom.c"
#include "buildTree.c"
#include "libkdtree.c"
#include "randGen.c"


int init_gsl_randgen(){
  gsl_rng_env_setup();  // random generation setup.

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);  // returns a pointer to a newly-created instance of a random number generator of type T

  return 0;
}


int set_file_paths(char survey[]){
  sprintf(surveyType, survey, fieldFlag, fkpPk, lo_zlim, hi_zlim);

  sprintf(logfilepath, "%s/W1_Spectro_V7_4/Qmultipoles/logfiles/%s.dat", root_dir, surveyType);

  return 0;
}


int set_files(count_res){
  switch(count_res){
  case 0:
    set_file_paths("maskmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_hihiRes_hex");
    break;
    
  case 1:
    set_file_paths("maskmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_hiRes_hex");
    break;
    
  case 2:
    set_file_paths("maskmultipoles_W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_loRes_hex");
    break;

  case 3:
    set_file_paths("maskmultipoles_W1W4_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_hihiRes_hex");
    break;
    
  case 4:
    set_file_paths("maskmultipoles_W1W4_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_hiRes_hex");
    break;
    
  case 5:
    set_file_paths("maskmultipoles_W1W4_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_%.2f_xi_%.1f_%.1f_loRes_hex");
    break;
  }

  return 0;
}


int initi_dist_z(){
  comovDistReshiftCalc(); // Cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h

  loChi                 = interp_comovingDistance(lo_zlim);
  hiChi                 = interp_comovingDistance(hi_zlim);

  printf("\n\nRedshift limits, lower bound: %.4lf \t %.4lf h^-1 Mpc, \n\t\t upper bound: %.4lf \t %.4lf h^-1 Mpc", lo_zlim, loChi, hi_zlim, hiChi);

  return 0;
}


int new_logfile(char logfilepath[]){
  int  status;
  char command[200];

  sprintf(command, "rm -f %s", logfilepath);

  status = system(command);

  return status;
}


int fprintf_time(char logfilepath[]){
  time_t        timer;  
  char    buffer[200];
  struct tm*  tm_info;
  
  time(&timer);

  tm_info = localtime(&timer);

  strftime(buffer, 26, "%Y:%m:%d %H:%M:%S", tm_info);

  output = fopen(logfilepath, "a");

  fprintf(output, "Time: \t");

  fprintf(output, "%s \n", buffer);

  fclose(output);

  return 0;
}


int main(int argc, char **argv){
  sprintf(root_dir,      "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");

  int        count_res;
  double sampling_frac;
  
  fieldFlag                 =        atoi(argv[1]);
  lo_zlim                   =        atof(argv[2]);      // previously 0.6<z<0.9, 0.7<z<1.1
  hi_zlim                   =        atof(argv[3]);
  maxlog                    =  log10(atof(argv[4]));    // hiRes: 20., lowRes: 2000.
  sampling_frac             =        atof(argv[5]);
  count_res                 =        atoi(argv[6]);

  
  AxisLimsArray[0][0]       =        0.0; 
  AxisLimsArray[1][0]       =      800.0;
 
  AxisLimsArray[0][1]       =        0.0;   
  AxisLimsArray[1][1]       =      800.0;                                                     // h^-1 Mpc                                                 

  AxisLimsArray[0][2]       =        0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]       =      800.0;
  
  if(fieldFlag == 1){
    // W1 catalogue. Nagoya v6 spectroscopic mask (& Samhain).
    // Neglects boundary issue in mocks. 
    LowerRAlimit            =     30.175; 
    UpperRAlimit            =     38.797;
    CentreRA                =     34.492;

    LowerDecLimit           =     -5.970;     
    UpperDecLimit           =     -4.171;     
    CentreDec               =     -5.091; 
    
    // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
  }
  
  else if(fieldFlag ==4){
    // W4 catalogue. Nagoya v6 spectroscopic mask (& Samhain). parent boundary limits. 
    LowerRAlimit            =    330.046; 
    UpperRAlimit            =    335.389;
    CentreRA                =    332.638;
  
    LowerDecLimit           =      0.862;     
    UpperDecLimit           =     2.3696;     
    CentreDec               =      1.583; 
    
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

  // Comoving number density, n(z), measurement. Change to equal increments in volume?
  chi_interval              =     16.00;
  
  // Change from 150 to match Stefano's 100.
  nz_smoothRadius           =     100.0;

  // FKP p(k) of interest;
  fkpPk                     =    8000.0;            // [h^-1 Mpc]^3.

  CatalogNumber             =       306;
      
  // analysis of VIPERS data or mock catalogues.
  data_mock_flag            =         0;
  
  // Correlation fn's, logarithmic binning in r. 
  zerolog  =               log10(0.001);
  logbinsz =               log10( 1.01);    // previously 1.4, must be >1.0 otherwise log gives 0. or -ve.
  
  nlogbins =  (int) ceil((maxlog - zerolog)/logbinsz);

  // linear binning in mu. 
  zerolin  =                   0.00;
  maxlin   =                   1.00;
  linbinsz =                   0.05;

  nlinbins =  (int) ceil((maxlin - zerolin)/linbinsz);


  init_gsl_randgen();

  initi_dist_z();
  
  UniverseAge();
  
  linearGrowthRate();
      
  loopCount = 10;  // choice of mock for loading nbar estimate.

  prep_nbar();
  
  spline_nbar(0);  // 306 mock avg; load nbar jointly estimated from W1 and W4, with 150 h^-1 Mpc smoothing and reflection.  

  pt2nz = &interp_nz; // Gaussian smoothed galaxy counts.

  prep_inverseCumulative_nbar();
  
  printf("\n\nRun details: field flag: %d, lo z: %.1lf, hi z: %.1lf", fieldFlag, lo_zlim, hi_zlim);
  
  printf("\n\nMax log: %.4lf, Sampling frac.: %.4lf, count res: %d", maxlog, sampling_frac, count_res);
  
  set_files(count_res);

  new_logfile(logfilepath);
  
  fprintf_time(logfilepath);
    
  load_homogeneous_rands_window(sampling_frac); // load randoms with sampling sampling_frac.
  
  assignMemory_xi();

  randWindow_pairCount();  
  
  fprintf_time(logfilepath);
  
  printf("\n\n");

  return 0;
}
