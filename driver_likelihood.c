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
#include "MultipoleCalc.c"
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
#include "ChiSq_input.c"
#include "calc_model.c"
#include "A11_vals.c"
#include "Alcock_Paczynski.c"
#include "Ruiz.c"
#include "combining_clipped_fsig8.c"
#include "freeMemory.c"


int main(int argc, char **argv){
  sprintf(root_dir,              "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir,         "/home/mjw/HOD_MockRun/W1_Spectro_V7_2"); 
  sprintf(covariance_mocks_path, "/home/mjw/HOD_MockRun/W1_Spectro_V7_3");
  sprintf(maskmultipoles_path,   "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");
  
  fieldFlag                 =       atoi(argv[1]);
  appliedClippingThreshold  =       atof(argv[2]);
  lo_zlim                   =       atof(argv[3]);   // previously 0.6<z<0.9, 0.7<z<1.1
  hi_zlim                   =       atof(argv[4]);
  Jenkins_foldfactor        =       atof(argv[5]);   // Apply Jenkins folding to increase spatial resolution of mesh.
  ChiSq_kmax                =   0.1*atof(argv[6]);   // Easier to pass int in bash. 

  W1area                    =              10.692;   // Nagoya v7 - overlapping Samhain v7; v6 and v7 has identical area. 
  W4area                    =               5.155;   // Don't move!
  TotalW1W4area             =      W1area + W4area;
  
  if(fieldFlag == 1){
    LowerRAlimit            =     30.175;            // Navgoya v6 + Samhain 
    UpperRAlimit            =     38.797;            // data:
    CentreRA                =     34.492;            // W1area = 10.763; W4area = 5.155;

    LowerDecLimit           =     -5.970;     
    UpperDecLimit           =     -4.171;     
    CentreDec               =     -5.091; 
    
    stefano_trans_x         =       +50.; 
    stefano_trans_y         =      +250.; 
    stefano_trans_z         =     -1425.;
    
    fracArea                = W1area/TotalW1W4area;
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

    fracArea                = W4area/TotalW1W4area;
  }

  data_mock_flag            =         0;                  // analysis of VIPERS data or mock catalogues. 
  
  min_bsigma8               =      0.05;                  // FOR GRANETT 2D POSTERIOR.
  max_bsigma8               =      1.00;                  // Previously 0.2 < b \sig_8 < 1.6
                                                          // 0.05 < b s8 < 1.0 (13/02/17)   
  min_fsigma8               =      0.05;                  //  Priors on the model params.  
  max_fsigma8               =      0.80;

  min_velDisperse           =      0.00;                  // CHANGED FROM 0.00 13/02/2017
  max_velDisperse           =      6.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)

  min_alpha_pad             =    0.9999;
  max_alpha_pad             =    1.0001;

  min_epsilon_pad           =   -0.0001;
  max_epsilon_pad           =    0.0001;
 
  min_A11Sq                 =      0.99;                  // distinct from bias due to spectral distortion. 
  max_A11Sq                 =      1.01;

  paramNumber               =       3.0;  // # of fitted params. -> dof in X^2. 

   Res                      =        16;  // Likelihood resolution [voxel number].
  dRes                      =      16.0;  // Previously 16: 13/02/17

   Res_ap                   =         1;  // Resoltuion in AP.
  dRes_ap                   =       1.0;

  FFTlogRes                 =       768; // FFTlogRes = 4096;

  ChiSq_kmin                =      0.02;
  hiMultipoleOrder          =         2;  // Fit monopole (1) or mono + quad (2).
  jenkins_fold_kjoin        =       0.4;  // k at which P(k) switches from unfolded to folded.     

  kbin_no                   =        40;
  modkMax                   =      1.00;            // Binned P(k) variables 
  muBinNumb                 =       100;

  clipping_smoothing_radius =       2.0;
   

  double  begin = getRealTime();

  set_chiSq_intervals(); // set e.g. fsig8 interval = (max - min)/interval.
  
  fftw_init_threads();

  fftw_plan_with_nthreads(omp_get_max_threads()); // Maximum number of threads to be used; use all openmp threads available. 

  // Cosmology from cosmology_planck2015.h or cosmology_valueaddedmocks.h// Selection parameters.
  comovDistReshiftCalc();

  // %% Need to Taylor to 0.9<z<1.2 %%
  // inputHODPk();      //  Must match Cosmology in cosmology_planck2015.h, or cosmology_valueaddedmocks.h
  inputLinearPk();      //  This is NOT automatically ensured.               
  
  prep_pkRegression(-2., log10(modkMax), kbin_no, 0);
  
  load_CovarianceMatrix(305, 1); // LOADING FROM W1_SPECTRO_V7_3.  Number of mocks, starting mock. 
  
  prep_VIPERS_maskMultipoles();
  
  prep_VIPERS_jmaskMultipoles();
  
  assign_LikelihoodMemory();  // Assigns memory for xdata, ydata, xtheory, ytheory, ChiSqGrid.

  assign_prep_FFTlog_memory(); // Memory for arrays speeding up FFTlog calc.

  FFTlog_memory(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);

  kvals_matchup();  // Given the k values of the measured P(k) multipoles, find nearest fftlog modes. 

  // default_params();
  // model_compute(0, 0, 0, 0, 0);
  
  calc_models();
  // read_models();
  
  calc_ChiSqs(20);
  
  // test_chiSq();
  
  double  end = getRealTime();

  printf("\n\nWall time: %.6lf", end - begin);
  
  printf("\n\n");

  return 0; 
}
