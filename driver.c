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

#include "Scripts/KaiserMultipoles.c"
// #include "Scripts/KaiserGaussMultipoles.c"
#include "Scripts/KaiserLorentzMultipoles.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"

#include "Scripts/nbar.c"
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

#include "Scripts/randGen.c"

#include "Scripts/FFT_log.h"
#include "Scripts/FFT_log.c"

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

#include "Scripts/halomodel_pk.c"
#include "Scripts/VIPERS_window.c"

// #include "Scripts/libkdpoly.h"
// #include "Scripts/libkdpoly.c"
// #include "Scripts/spec_weights.c"
// #include "Scripts/slowDFT.c"

// #include "tinker.c"

#include "Scripts/freeMemory.c"


int main(int argc, char **argv){
  sprintf(root_dir,      "/disk1/mjw/HOD_MockRun");

  sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_HOD_Mocks");
  // sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_ValueAddedHOD");
 
  // Stacpolly run. 
  // sprintf(root_dir, "/home/mjw/HOD_MockRun");

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
  AxisLimsArray[0][0]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][0]   =     800.0;
 
  AxisLimsArray[0][1]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][1]   =     800.0;

  AxisLimsArray[0][2]   =       0.0;                                                     // h^-1 Mpc
  AxisLimsArray[1][2]   =     800.0;

  // degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
  stefano_trans_x       =     +250.;
  stefano_trans_y       =     +250.;
  stefano_trans_z       =    -1700.;
  
  /*
  // W1 catalogue. new 500s mocks. parent
  LowerRAlimit          =       30.1; 
  UpperRAlimit          =       38.8;
  CentreRA              =      34.45;

  LowerDecLimit         =      -5.95;     
  UpperDecLimit         =      -4.15;     
  CentreDec             =      -5.05;  
  
  W1area                =      15.66;      
  */
  
  // W1 catalogue. new 500s mocks. Nagoya v4 spectroscopic mask. 
  LowerRAlimit          =      30.17; 
  UpperRAlimit          =       38.8;
  CentreRA              =     34.487;
  
  LowerDecLimit         =      -5.38;     
  UpperDecLimit         =      -4.17;     
  CentreDec             =      -4.77;     
  
  W1area                =      7.074;  // sq. degs.
  
  // Cell size, comoving distance, h^-1 Mpc. 
  fft_size              =        512;     

  // Selection parameters.
  lo_MBlim              =     -90.5;                  // -20.5 < M_B < -19.5
  hi_MBlim              =     -00.0;
  
  lo_zlim               =      0.70;                  // previously 0.6<z<0.9, 0.7<z<1.1 
  hi_zlim               =      1.10;
  z_eff                 =      0.85;                  // set to volume avg. redshift of the survey? prior to FKP weights.

  linearBias            =      1.53;                  // 1.32, appropriate for 0.6 to 0.9;
  velDispersion         =       3.0;                  // units of h^-1 Mpc rather than 300 km s^-1
  beta                  =     0.541;                  // 2dF measurement, beta = 0.43, beta = 0.542 for the cube. 

  // Priors on the model params.
  min_fsigma8           =      0.20;
  max_fsigma8           =      0.80;

  min_velDisperse       =      0.00;
  max_velDisperse       =      6.00;

  // 4 param likelihood required. change this. 
  min_bsigma8           =       0.8;
  max_bsigma8           =       1.6;
 
  // distinguished from linear bias by spectral distortion.  
  min_A11Sq             =      0.99;
  max_A11Sq             =      1.01;

  // Number of fitted parameters, defines degrees of freedom in chi sq. expectation. 
  paramNumber           =       3.0;

  // Resolution of the Likelihood evaluation [voxel number].
  Res                   =        16;
  dRes                  =      16.0;

  ChiSq_kmin            =     0.025;
  ChiSq_kmax            =      0.40; 
     
  // Fit solely the monopole (1) or both mono and Quad (2).
  hiMultipoleOrder      =         2;

  // Comoving number density, n(z), measurement. Change to equal increments in volume?
  chi_interval          =     16.00;

  // Apply Jenkins contraction to beat aliasing. 
  JenkinsScalefactor    =       1.0;

  // FKP P(k) of interest;
  fkpPk                 =    1000.0;            // [h^-1 Mpc]^3.
  meanSampling          =       0.4;

  // Binning interval for P(k).
  modkMax               =      1.00;
  muBinNumb             =       100;
  kbin_no               =        40;

  // Interval in k^2 for perp k binning of 2D P(k).
  perpkInterval         =      0.01;
 
  // Total number of HOD mocks. 
  CatalogNumber         =       306;

  //  Apodise the window fn. to supress the Gibb's phenomenon.
  GibbsSkinDepth        =       5.0;
    
  // Clipping variables. 
  appliedClippingThreshold  =  40.0;    

  // nfw halo generation
  nfw_conc              =       1.0; 

  // Correlation fn's, logarithmic binning in r. 
  zerolog  =             log10(0.001);
  maxlog   =             log10(20.0);       // hiRes: 20., lowRes: 2000.
  logbinsz =             log10( 1.01);      // previously 1.4, must be >1.0 otherwise log gives 0. or -ve.
  
  nlogbins =  (int) ceil((maxlog - zerolog)/logbinsz);

  // linear binning in mu. 
  zerolin  =                   0.00;
  maxlin   =                   1.00;
  linbinsz =                   0.05;

  nlinbins =  (int) ceil((maxlin - zerolin)/linbinsz);

  // Random variable generation.    
  gsl_rng_env_setup();

  gsl_ran_T = gsl_rng_default;
  gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

  
  comovDistReshiftCalc();

  // JenkinsCoordinates();

  EvaluateGridParameters();
  
  // cosmology has changed, recalculate linear power spectrum using CAMB. 
  // inputLinearPk();     
  
  // prep_NFWhaloCat(10000, 100000);

  // printf("\n%e", haloModel_pk(0.1, 0.1, 0.5, 0));
  
  // assign2DPkMemory(muBinNumb, kBinNumb);   

  // knownGRF_mask();

  // generate randoms for the mask in a given cosmology. celestial coordinates (ra and dec, spherical)
  // randoms_maskGen();
  
  // sampling: 0.09 for hi res measurement, ? for lo res, 1.00 for creating mask for P(k), free rands for P(k).
  load_homogeneous_rands_window(1, 0.02);
  
  // calculate multipole moments of the window. 
  randWindow_pairCount();
  
  // calculate nbar, in a given cosmology: nbar(chi)
  // nbar_calc(306);
  /*
  spline_nbar();
  
  // calculate volume avg. shot noise for shot noise correction. 
  volavg_invnbar  = calc_volavg_invnbar();
  
  volavg_redshift = interp_inverseComovingDistance(calc_volavg_chi());
  
  printf("\n\nvol. avg. redshift: %.3lf (%.3f h^-1 Mpc)", volavg_redshift, calc_volavg_chi());
  
  // add_fkp();
  
  // correct for the amplitude factor caused by the mask. 
  // calc_DigitalAmplitudeCorrection(); 
  calc_AnalogAmplitudeCorrection();
  */
  // apodiseWindowfn();

  // apodisedVolume        = CellVolume*SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);
  
  // prepBootStrap(n0*n1*n2, Cell_rotatedXvals, Cell_rotatedYvals, Cell_rotatedZvals, 1000.);
  
  UniverseAge();
  
  linearGrowthRate();
  
  inputHODPk();
  
  // Must be uncommented for fftlog calcs. 
  prep_VIPERS_maskMultipoles();         
  
  // VIPERS_mask_cnvldpk();
  
  // VIPERS_mask_intCnsrt();
  
  // halofit();

  // Number of randoms, load positions, load RR.
  // prep_randpairs(20.*30000., 1, 0);

  // initialise_angularCorrelation();

  // randoms_angular_correlationfn();

  // computeCorrelation_fns(2);

  // wfPkCalc();

  // prep_pairwisepdf();

  // print_windowCorrfn();

  // CovarianceMatrix(306, 1);
  
  // ChiSq_minimisation();
  
  // mass_fn();
  
  // halobias_fn();
  
  // polymaskCalc();
  
  for(loopCount=1; loopCount<1; loopCount++){
    printf("\n\n%d", loopCount);
    
    // sprintf(surveyType, "250_fullCube_Gaussian_fogRSD_CellSize_%.2f", CellSize);

    // if(loopCount<10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    // else              sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);
    
    // new 500s. limit to 0.7<z<0.8, limit to linear bias of 1.495903 corresponding to ~ -20.0 mag galaxies. at this mag. volume limited to z=0.85
    if(loopCount<10)        sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_00%d_spec.dat", vipersHOD_dir, loopCount);
    else if(loopCount<100)  sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_0%d_spec.dat",  vipersHOD_dir, loopCount);
    else                    sprintf(filepath, "%s/mocks_W1_v8.0_500/mock_%d_spec.dat",   vipersHOD_dir, loopCount);
    
    CatalogueInput_500s(filepath);
    
    // Choice of redshift from zcos, zpec, zphot, zobs.
    gal_z = &zobs[0];
    
    // Set redshift and absolute mag. cuts. 
    assignAcceptance();
    
    // normalised with respect to accepted galaxies.
    // spec_weights();
    
    // Convert from (ra, dec, redshift) to (x, y, z) in Stefano's basis. basis choice must be consistent with that used for the mask defined by the randoms. 
    StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
    
    // fiberCollision_cat(rand_number, rand_redshift, rand_x, rand_y, rand_z);
    
    // fiberCollision_cat(Vipers_Num, zUtilized, xCoor, yCoor, zCoor);

    // correlationfn();

    // angular_correlationfn();

    // randomGeneration();
    
    // loadRand();
    
    // projectVIPERSsystem();

    // VIPERSbasis(CentreRA, CentreDec, xCoor, yCoor, zCoor, Vipers_Num);
    
    // Gaussianfield();
    
    // deplete_homogeneous();
    
    // underclipDensity(appliedClippingThreshold);
    
    // poissonSample_lnNorm();
    
    // Do not forget call to prep Cat. 
    // haloCatalogue_nfwprofile(3000000);
    
    // Note: Dangerous non-commutation of sampling and randomise in load_clustered. may generate biased population if not careful. //
    // load_clustered(1, 0.01);

    // NFW_profile_pairCount();

    // Calculation of multipole moments for the halo model catalogue.
    // ukm_calc();
    
    // CleanNGP();
    
    // Checked.
    // NGPCalcCube(xCoor, yCoor, zCoor, Vipers_Num);
    // NGPCalcCube(rand_x, rand_y, rand_z, rand_number);
    
    // redshift range selection handled by randoms mask. 
    calc_overdensity();

    // clipDensity(appliedClippingThreshold);
    
    // Apply FKPweights(meanSampling);
  
    PkCalc();
    
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

  // sprintf(filepath, "%s/Data/SpectralDistortion/measuredCorrelationfn_mono_randoms.dat", root_dir);

  // correlationfn(filepath, rand_x, rand_y, rand_z, rand_number);

  // fclose(output);

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
