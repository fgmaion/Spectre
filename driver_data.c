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
#include <gsl/gsl_deriv.h>

// #include "omp.h"

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#define   AUXfn_funcs  "/disk1/mjw/Aux_functions/Aux_functions.c"
#include  AUXfn_funcs

#include "Scripts/header.h"

#include "Scripts/cosmology_planck2015.h"
// #include "Scripts/cosmology_valueaddedmocks.h"

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

#include "Scripts/fkp_weights.c"
#include "Scripts/clipping_weights.c"

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
// #include "Scripts/fkp.c"
#include "Scripts/Bailey.c"
#include "Scripts/Alcock_Paczynski.c"

#include "Scripts/paircount_mask.c"

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
  AxisLimsArray[0][0]   =        0.0; 
  AxisLimsArray[1][0]   =      800.0;
 
  AxisLimsArray[0][1]   =        0.0;   
  AxisLimsArray[1][1]   =      800.0;                                                     // h^-1 Mpc                                                 

  AxisLimsArray[0][2]   =        0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]   =      800.0;
  
  fieldFlag             =          1;
  
  if(fieldFlag == 1){
    // W1 catalogue. Nagoya v6 spectroscopic mask (& Samhain).
    LowerRAlimit        =      30.17; 
    UpperRAlimit        =      38.80;
    CentreRA            =     34.487;

    LowerDecLimit       =      -5.95;     
    UpperDecLimit       =      -4.17;     
    CentreDec           =      -5.06; 
    
    // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
    stefano_trans_x     =       +50.; 
    stefano_trans_y     =      +250.; 
    stefano_trans_z     =     -1425.;
  }
  
  else if(fieldFlag ==4){
    // W4 catalogue. Nagoya v6 spectroscopic mask (& Samhain). parent boundary limits. 
    LowerRAlimit        =      330.0; 
    UpperRAlimit        =      335.5;
    CentreRA            =     332.75;
  
    LowerDecLimit       =       0.80;     
    UpperDecLimit       =       2.40;     
    CentreDec           =       1.60; 
    
    // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
    stefano_trans_x     =       +50.; 
    stefano_trans_y     =      +250.; 
    stefano_trans_z     =     -1425.;
  }
  
  // official release info, (Nagoya v7 - overlapping Samhain area). With regards angular coverage v7 and v6 are identical. 
  W1area                =      10.763;  // sq. degs.
  W4area                =       5.155;

  // Required for <n(z)> calculation.
  TotalW1W4area         =      W1area + W4area; 
  
  fft_size              =        256;     

  // Selection parameters.
  lo_MBlim              =     -90.5;                  // -20.5 < M_B < -19.5
  hi_MBlim              =     -00.0;
  
  lo_zlim               =      0.60;                  // previously 0.6<z<0.9, 0.7<z<1.1 
  hi_zlim               =      0.90;
  z_eff                 =      0.75;                  // set to weighted average of all galaxies in sample (Anderson et. al.)

  linearBias            =      1.53;                  // 1.32, appropriate for 0.6 to 0.9;
  velDispersion         =       3.0;                  // units of h^-1 Mpc rather than 300 km s^-1
  beta                  =     0.541;                  // 2dF measurement, beta = 0.43, beta = 0.542 for the cube. 

  // Priors on the model params.
  min_fsigma8           =      0.00;
  max_fsigma8           =      0.80;

  min_velDisperse       =      0.00;
  max_velDisperse       =      6.00;

  min_bsigma8           =       0.2;
  max_bsigma8           =       1.6;
 
  // distinguished from linear bias by spectral distortion.  
  min_A11Sq             =      0.99;
  max_A11Sq             =      1.01;

  // Number of fitted parameters, defines degrees of freedom in chi sq. expectation. 
  paramNumber           =       3.0;

  // Resolution of the Likelihood evaluation [voxel number].
   Res                  =        16;
  dRes                  =      16.0;

  ChiSq_kmin            =      0.09;
  ChiSq_kmax            =      0.30; 
     
  // Fit solely the monopole (1) or both mono and Quad (2).
  hiMultipoleOrder      =         2;

  // Comoving number density, n(z), measurement. Change to equal increments in volume?
  chi_interval          =     16.00;
  nz_smoothRadius       =     150.0;

  // Apply Jenkins folding to increase spatial resolution of mesh. 
  Jenkins_foldfactor        =       1.0;

  // fkp p(k) of interest;
  fkpPk                     =    8000.0;            // [h^-1 Mpc]^3.

  // binned p(k) variables.
  modkMax                   =      1.00;
  muBinNumb                 =       100;
  kbin_no                   =        40;            

  // Total number of HOD mocks. 
  CatalogNumber             =       306;

  //  Apodise the window fn. to supress ringing of the window ("Gibb's phenomenon").
  GibbsSkinDepth            =       5.0;
    
  // Clipping variables. 
  appliedClippingThreshold  =       5.0;
  
  clipping_smoothing_radius =       2.0;
  
  // nfw halo generation
  nfw_conc                  =       1.0; 

  // analysis of VIPERS data or mock catalogues.
  data_mock_flag            =         1;

  // Random variable generation.    
  gsl_rng_env_setup();

  gsl_ran_T = gsl_rng_default;
  
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

  //** cosmology determined by inclusion of cosmology_planck2015.h or cosmology_valueaddedmocks.h **//
  comovDistReshiftCalc();
  
  Jenkins_foldEmbeddingVol();
  
  EvaluateGridParameters();
  
  // AgeOftheUniverse.c
  UniverseAge();
  
  // linearGrowthRate.c
  linearGrowthRate();
  
  // Must match Cosmology declared in cosmology_planck2015.h, cosmology_valueaddedmocks.h //
  // This is NOT automatically ensured. //
  
  inputHODPk();
  // inputLinearPk();
  
  // Must be uncommented for fftlog calcs. 
  // prep_VIPERS_maskMultipoles(data_mock_flag);     
  
  // prep_NFWhaloCat(10000, 100000);
  
  // assign2DPkMemory(muBinNumb, kBinNumb);   

  // knownGRF_mask();
  
  // calculate nbar, in a given cosmology: nbar(chi)
  // nbar_calc(306);
  
  // recall if nbar changes, redshift limits change etc. 
  // z_eff = zeff_calc();
  
  // ^^Turn off Jenkins folding. ^^//
  // randoms_maskGen();
  
  /*
  // number of mocks, starting index. 
  load_CovarianceMatrix(293, 1);
  
  W1_Spectro_V7_ChiSq_minimisation();
  */
  /*
  // number of mocks to be analysed, starting index. mocks are selected at random.
  ensemble_fsig8(30, 1);
  */
  // localTSR_calc_Data();
  
  // mask_area();
  
  // quadrant_sampling();
  
  // assign2DPkMemory();
  
  // toySpoc();

  // AP_correction();
  
  // randoms_slowDFTcalc();
  
  // modeSampling(4000.);
  
  // vipers_fkpCalc(4000.);
  
  // wfPkCalc();
  
  // printf("\n %d", atoi(argv[1]));
  
  // int start = 298;
  
  // clipped_lnnorm_pkCalc();
  
  // 'truth' from averaging over mocks, no smoothing or reflection.
  // nbar_calc(293); 
  
  // Smoothing width, reflection flag. 
  // if(data_mock_flag==0)  for(loopCount=1; loopCount<295; loopCount++)  smoothed_nbar_calc(nz_smoothRadius, 1);
  // if(data_mock_flag==1)  smoothed_nbar_calc(nz_smoothRadius, 1);
  
  // calculate multipole moments of the window. 
  
  // assigns memory for overdensity grid.
  prep_grid();
  
  prep_mask(); 
  
  prep_fftw();
  
  prep_pkRegression(-2., log10(modkMax), kbin_no);
  
  int start = atoi(argv[1]);
  
  paircount_mask(start);
  
  for(loopCount=1; loopCount<1; loopCount++){
    double fkp_norm;
  
    printf("\n\n%d", loopCount);
  
    // data_mock_flag: Nagoya v6 mock catalogues: 0, VIPERS data: 1. Second argument: Gaussian smoothed, Reflected 2-field average for mocks, third argument: 'truth' obtained from averaging mocks.
    spline_nbar(data_mock_flag, 1, 0);
    
    // sampling: 0.01 for lo res, 0.05 for hi res measurement, 0.7 for hihRes measurement. 1.00 for creating mask for P(k), remember to free rands after grid assignment for P(k) measurement.
    load_homogeneous_rands_window(1,  0.1, data_mock_flag); 
    
    Jenkins_foldRand();
    
    if(data_mock_flag==0){
        // mock catalogues.
        if(loopCount<10)        sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_00%d_spec_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, loopCount);
	    else if(loopCount<100)  sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_0%d_spec_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, loopCount);
	    else                    sprintf(filepath, "%s/mocks_W%d_v8.0_500/mock_%d_spec_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, loopCount);
    
	    CatalogueInput_500s(filepath);
    }
    
    if(data_mock_flag==1){
        // %% W1_SPECTRO_V7_0.txt %%//
        DataInput(filepath);
    }
    
    // Choice of redshift from zcos, zpec, zphot, zobs.
    gal_z = &zobs[0];
    
    // load sampling according to local TSR. 
    spec_weights();
    
    if(data_mock_flag==0)  assignAcceptance();
    if(data_mock_flag==1)  assignAcceptance_WX_SPECTRO_V7();
    
    // Convert from (ra, dec, redshift) to (x, y, z) in Stefano's basis. basis choice must be consistent with that used for the mask defined by the randoms. 
    StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
    /* 
    Jenkins_foldCat();
    
    // Normalisation of FKP weights set by random catalogue. 
    fkp_norm = calc_fkpweights();
    
    // Apply this normalisation to gal. weights. 
    set_gal_fkpweights(fkp_norm);
    
    // calc_clippingweights();
    for(j=0; j<Vipers_Num; j++)  clip_galweight[j] = 1.;
        
    // fiberCollision_cat(rand_number, rand_redshift, rand_x, rand_y, rand_z);
    
    // fiberCollision_cat(Vipers_Num, zUtilized, xCoor, yCoor, zCoor);

    // correlationfn();

    // angular_correlationfn();

    // randomGeneration();
    
    // loadRand();
    
    // projectVIPERSsystem();

    // VIPERSbasis(CentreRA, CentreDec, xCoor, yCoor, zCoor, Vipers_Num);
    
    // deplete_homogeneous();
    
    // underclipDensity(appliedClippingThreshold);
    
    // poissonSample_lnNorm();
    
    // Do not forget call to prep Cat. 
    // haloCatalogue_nfwprofile(3000000);
    
    // Note: Dangerous non-commutation of sampling and randomise in load_clustered. may generate biased population if not careful. 
    // load_clustered(1, 0.01);

    // NFW_profile_pairCount();

    // Calculation of multipole moments for the halo model catalogue.
    // ukm_calc();
    
    // CleanNGP();
    
    // NGPCalcCube(xCoor, yCoor, zCoor, Vipers_Num);
    // NGPCalcCube(rand_x, rand_y, rand_z, rand_number);    
    
    // redshift range selection handled by randoms mask. 
    calc_overdensity();
    
    // Gaussianfield();
       
    // lnNormfield();
    
    PkCalc();
    */
    // slowDFTcalc();
  }
  
  // spherical_randDistribution();

  // print_pairwisepdf();

  // hodmodel_xi();

  // calc_zeldovichxi();

  // fog_xiCalc();

  // clippedPkCalc();

  // lnnormPkCalc();

  // spline_realSpaceCorrfn(55);

  // fprintf_fiberCollision();

  // SaundersDeproject(40, 0.1);

  // correlationfn(filepath, rand_x, rand_y, rand_z, rand_number);

  // MockAvg2Dpk(14);

  // MockAverageComovingdensity();

  // fitnz();

  // freeHOD();

  // freeNGP();

  // printHODMatterPk();

  // wfPkCalc();

  // analyticConvTest();

  // AnisoConvolution();
 
  // formPkCube();

  // theoryHexadecapole();

  // multipolesRealisation();

  // Stacpolly run. 
  // MPI_Finalize();
  
  printf("\n\n");

  return 0; 
}
