#define  KBIN_NO 40

#include "/home/mjw/Aux_functions/header.h"
#include "/home/mjw/Aux_functions/Aux_functions.c"

#include "header.h"
#include "FFT_log.h"
#include "header_W2.h"
#include "header_chi2.h"
#include "struct_regress.h"
// #include "struct_multipoles.h"
#include "cosmology_planck15.h"

#include "AgeOftheUniverse.c"
#include "linearGrowthRate.c"
#include "angular_limits.c"
#include "chi_zcalc.c"
#include "assign_binnedpk_memory.c"
#include "likelihood_memory.c"
#include "KaiserMultipoles.c"
#include "KaiserLorentzMultipoles.c"
#include "FFT_log.c"
#include "viperskwindow.c"
#include "viperswindow.c"
#include "vipersjwindow.c"
#include "FFTlog_vipers_clipping.c"
#include "toymodel_pk_xi.c"
#include "linear_pk.c"
#include "nonlinear_pk.c"
#include "get_pkderivedprops.c"
#include "get_allkvals.c"
#include "super_vipers.c"
#include "MultipoleCovariance_eigenvecs.c"
#include "MultipoleCovariance.c"
#include "ChiSq_input.c"
#include "calc_model.c"
#include "ChiSq_minimisation.c"
#include "onedposteriors.c"
#include "clip_icc.c"
#include "lock_files.c"
#include "AlcockPaczynski.c"
// #include "FisherForecast.c"
#include "FisherForecastFAP.c"
#include "Ruiz.c"
#include "camb_call.c"
// #include "combining_clipped_fsig8.c"


int get_main(int set_d0, int set_fieldFlag, double set_lo_zlim, double set_hi_zlim, double set_ChiSq_kmax){  
  thread                    =    1;
  data_mock_flag            =    1;  // analysis of VIPERS data or mock catalogues. 

  z_eff                     =        atof(getenv("ZEFF"));
  outputdir                 =         getenv("outputdir");
  maskmultipoles_path       =        getenv("mask_Qldir");
  
  sprintf(root_dir,                              "/home/mjw/HOD_MockRun");
  sprintf(vipersHOD_dir,         "/home/mjw/HOD_MockRun/W1_Spectro_V7_2"); 
  sprintf(models_path,                                         outputdir);

  // sprintf(covariance_mocks_path,                               outputdir);
  sprintf(covariance_mocks_path, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2");
  
  d0                        =               set_d0;
  fieldFlag                 =        set_fieldFlag;
  lo_zlim                   =          set_lo_zlim;   
  hi_zlim                   =          set_hi_zlim;
  ChiSq_kmax                =       set_ChiSq_kmax;

  
  min_bsigma8               =      0.05;                  // FOR GRANETT 2D POSTERIOR.
  max_bsigma8               =      1.00;                  // Previously 0.2 < b \sig_8 < 1.6
                                                          // 0.05 < b s8 < 1.0 (13/02/17)   
  min_fsigma8               =      0.00;                  // Priors on the model params.  
  max_fsigma8               =      0.80;

  min_velDisperse           =      0.00;                  // CHANGED FROM 0.00 13/02/2017
  max_velDisperse           =      7.00;                  // CHANGED FROM 6.00, 19 JAN. DIFFERS FROM MUNICH CLIPPED RESULTS (I GUESS)
  
  /*
  min_bsigma8               =      0.05;                
  max_bsigma8               =      3.50;                   
                                                          
  min_fsigma8               =      0.05;                  
  max_fsigma8               =      1.80;

  min_velDisperse           =      0.05;                   
  max_velDisperse           =     15.00;                     
  */
  min_alpha_pad             =    0.9999;
  max_alpha_pad             =    1.0001;

  min_epsilon_pad           =   -0.0001;
  max_epsilon_pad           =    0.0001;
 
  min_A11Sq                 =      0.99;                  // distinct from bias due to spectral distortion. 
  max_A11Sq                 =      1.01;

  paramNumber               =       3.0;  // # of fitted params. -> dof in X^2. 

   Res                      =         1;  // Likelihood resolution [voxel number].
  dRes                      =       1.0;  // Previously 16: 13/02/17

   Res_ap                   =         1;  // Resoltuion in AP.
  dRes_ap                   =       1.0;

  FFTlogRes                 =       768;  // FFTlogRes = 4096;
  // FFTlogRes              =      4096;
  
  logk_min                  =      -2.0;
  logk_max                  =   0.60206;  // k = 4 hMpc^{-1}.
  
  ChiSq_kmin                =      0.02;

  hiMultipoleOrder          =         2;  // Fit monopole (1) or mono + quad (2).
  jenkins_fold_kjoin        =       0.4;  // k at which P(k) switches from unfolded to folded.     
  modkMax                   =      1.00;
  
  smooth_radius             =       2.0;

  // Regression to ~ May 2016 -> Catalog number to 305; change file paths of mocks (covariance and chi sq input) and Qmultipoles. Change init_covariance.
  CatalogNumber             =       305; // 153 

  
  start_walltime();

  printf_branch();
  
  // fftw_init_threads();

  // fftw_plan_with_nthreads(omp_get_max_threads()); // Maximum number of threads to be used; use all openmp threads available. 

  set_angularlimits(0, fieldFlag);                   // Cut data to mock limits.
  
  chi_zcalc();            
  
  nonlinear_pk();
  // linear_pk(); 
  
  // get_mocksshotnoise();
  
  prep_FFTlog_memory();    // assign memory for arrays speeding up FFTlog calc; e.g. xi -> pre/post factors. 
  
  set_FFTlog(FFTlogRes, pow(10., -10.), pow(10., 14.), 1., velDispersion);  // assigns values to mono_config etc. 
  
  prep_VIPERS_maskMultipoles();
  
  prep_VIPERS_jmaskMultipoles();
  
  precompute_vipers_clipping_model(FFTlogRes);  // Computes P_R(k), W_0(r), ..., \tilde W_0(k), ..., and \tilde W_0(k) for the joint field.   
  
  get_allkvals(1);    // all kVals, ignoring ChiSq_kmin and ChiSq_kmax, but including folding. 
  
  allkvals_matchup(); // Match all available kVals to FFTlog modes; not just those up to ChiSq_kmax. 
  
  set_chiSq_intervals(); // set e.g. fsig8 interval = (max - min)/interval.
  
  assign_LikelihoodMemory();  // Assigns memory for xdata, ydata, xtheory, ytheory, ChiSqGrid.
  
  // get_mocksclippedamplitudes();  
  
  load_CovarianceMatrix(CatalogNumber, 1);
  
  delete_lockfile();
    
  prep_dlnPR_dlnk();

  // set_clippingvars(1.0);
  
  calc_models();

  kvals_matchup();  // Now match only available modes between ChiSq_kmin and ChiSq_kmax.
  
  set_models();
    
  prewhitenCov();  // Pre-whiten data and covariance.

  // scale_Cov(CatalogNumber);
  
  Covariance_eigenVecs(CatalogNumber);

  delete_lockfile();
  
  // default_params();
  
  // model_compute(0, 0, 0, 0, 0); 
  
  walltime("Wall time at finish");

  // MPI_Finalize();
  
  printf("\n\n");

  return 0;
}


int main(int argc, char** argv){
  (void) argc;

  d0           =       atoi(argv[1]);
  fieldFlag    =       atoi(argv[2]);
  lo_zlim      =       atof(argv[3]);
  hi_zlim      =       atof(argv[4]);
  ChiSq_kmax   =       atof(argv[5]);

  get_main(d0, fieldFlag, lo_zlim, hi_zlim, ChiSq_kmax);  
  
  return 0;
}
