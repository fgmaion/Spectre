// Stacpolly run. 
// #define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

#include <stdbool.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
// #include "omp.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_expint.h>

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#define   AUXfn_funcs  "/disk1/mjw/Aux_functions/Aux_functions.c"
#include  AUXfn_funcs

#include "Scripts/header.h"

#include "/disk1/mjw/Aux_functions/SVD.c"

#include "Scripts/comovDistRedshiftCalc.c"

#include "Scripts/JenkinsRun.c"
#include "Scripts/GridParams.c"

#include "Scripts/assignMemory.c"

#include "Scripts/CoordinateCalc.c"
#include "Scripts/assignAcceptance.c"

#include "Scripts/Padding.c"
#include "Scripts/NGPCalc.c"
#include "Scripts/CloudInCell.c"
#include "Scripts/BasisChange.c"
#include "Scripts/CalcCellraDec.c"
#include "Scripts/ApodiseWindowFunc.c"

#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

#include "Scripts/KaiserMultipoles.c"
#include "Scripts/KaiserGaussMultipoles.c"
// #include "Scripts/KaiserLorentzMultipoles.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"
#include "Scripts/FFTw_3Dwf.c"
#include "Scripts/axesWfSlices.c"
#include "Scripts/mixingMatrix.c"

#include "Scripts/MeasureWfKernel.c"
// #include "Scripts/sphericalConvolvePk.c"
#include "Scripts/setConvolutionKernels.c"
// #include "Scripts/AnalyticTestConvolution.c"
#include "Scripts/ConvolvePkAnisoWf.c"

#include "Scripts/IntegralConstraintCorrection.c"

#include "Scripts/ComovingNumberDensityCalc.c"
#include "Scripts/MockAvgComovingDensity.c"

#include "Scripts/AgeOftheUniverse.c"
#include "Scripts/linearGrowthRate.c"
#include "Scripts/growthfactor_derivative.c"

#include "Scripts/InvErrorfn.c"
#include "Scripts/MatrixInverse.c"

#include "Scripts/Clipped_zSpace.c"
#include "Scripts/toymodel_pk_xi.c"

#include "Scripts/PowellsRoutine.c"
#include "Scripts/Powells_mnbrak.c"
#include "Scripts/Powells_linmin.c"
#include "Scripts/Powells_Brent.c"

#include "Scripts/fitting_nz.c"
#include "Scripts/randGen.c"

#include "Scripts/ArtificialWf.c"
#include "Scripts/BootStrap.c"

#include "Scripts/correlation_fns.c"

#include "Scripts/redshiftDistribution_NagoyaRandoms.c"

#include "Scripts/FFT_log.h"
#include "Scripts/FFT_log.c"

#include "Scripts/cubature/cubature.h"
#include "Scripts/FFT_log_zeldovich.h"
#include "Scripts/FFT_log_zeldovich.c"

#include "Scripts/anisotropicGaussian_multipoles.c"

/*
#include "Scripts/slowDFT.c"

// #include "Scripts/MockAvgMultipole.c"
#include "Scripts/MultipoleCovariance.c"
#include "Scripts/MultipoleCovariance_Inverse.c"
#include "Scripts/MultipolesRealisation_MultiVariateGauss.c"
#include "Scripts/MarginalisedPosteriors.c"
#include "Scripts/2DPosteriors.c"
#include "Scripts/LikelihoodEval.c"
*/

#include "Scripts/freeMemory.c"

#include "Scripts/MonteCarlo_SSPOC.c"
#include "Scripts/AngularSelectionCats.c"
#include "Scripts/SaundersDeproject.c"

#include "Scripts/libkdtree.h"
#include "Scripts/kdtree_xi_mom.c"
#include "Scripts/buildTree.c"
#include "Scripts/libkdtree.c"

#include "Scripts/mockGalaxyCats.c"

#include "Scripts/NFW_profile.c"
#include "Scripts/VIPERS_window.c"


int main(int argc, char **argv){
sprintf(root_dir,      "/disk1/mjw/HOD_MockRun");
sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_ValueAddedHOD");

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


// Stefano basis. 
AxisLimsArray[0][0]   =       0.0;                                                     // h^-1 Mpc
AxisLimsArray[1][0]   =     500.0;

AxisLimsArray[0][1]   =       0.0;                                                     // h^-1 Mpc
AxisLimsArray[1][1]   =     500.0;

AxisLimsArray[0][2]   =       0.0;                                                     // h^-1 Mpc
AxisLimsArray[1][2]   =     500.0;

// degree of translation for Stefano's co-ordinates, fit survey into surrounding volume. 
stefano_trans_x       =     +250.;
stefano_trans_y       =     +250.;
stefano_trans_z       =    -1650.;

// stefano_trans_x       =     +100.;
// stefano_trans_y       =     +300.;
// stefano_trans_z       =    -1500.;

// limits in right asecension.
// W1 catalogue.
LowerRAlimit          =      30.1;
UpperRAlimit          =      38.8;
CentreRA              =     34.45;

// W4 catalogue.
// LowerRAlimit       =      330.0;
// UpperRAlimit       =      335.5;

// limits in declination. 
// W1 catalogue
LowerDecLimit         =      -5.4;
UpperDecLimit         =     -4.15;
CentreDec             =    -4.775;

// W4 catalogue.
// LowerDecLimit      =      0.82;
// UpperDecLimit      =      2.42;

// Total angular areas, 

// W1: 10.875 sq. degs.
W1area                =    10.875;

// W4:  8.8   sq. degs.
W4area                =       8.8;

// Combined. 
// 19.675 sq degs. Previously 21.47, 20% difference in P(k).
TotalW1W4area         =    19.675; 

// Cell size, comoving distance, h^-1 Mpc.
CellSize              =       2.0;     

// Selection parameters.
absMagCut             =     22.00;
redshiftLowLimit      =      0.60;
redshiftHiLimit       =      0.90;

// Non-linear RSD
velDispersion         =      12.00;                  // units of h^-1 Mpc rather than 300 km s^-1
beta                  =      0.542;                  // 2dF measurement, beta = 0.43, beta = 0.542 for the cube. 

// Priors on the model params.
min_beta              =      0.30;
max_beta              =      0.50;

min_velDisperse       =      1.60;
max_velDisperse       =      1.80;
 
min_A11Sq             =       2.2;
max_A11Sq             =       3.0;

// Resolution of the Likelihood evaluation [voxel number].
 Res                  =        40;
dRes                  =      40.0;

/* linearBias.txt in dir. /HODTheoryPk/
-22.0  2.924277
-21.5  2.199471
-21.0  1.824770
-20.5  1.618792
-20.0  1.495903
-19.5  1.415186
-20.75 1.707502
-20.25 1.550205
-20.15 1.527030
*/

// -20.0 mags lim. sample. See linearBias.txt in HODTheoryPk dir.
linearBias            =  1.495903;

// Comoving number density, n(z), measurement. 
zBinWidth             =      0.03; 
chiBinWidth           =     16.00;
nzSigma               =      50.0;

// Apply Jenkins contraction to beat aliasing. 
JenkinsScalefactor    =       1.0;

// FKP P(k) of interest;
fkpPk                 =    1000.0;                                                  // [h^-1 Mpc]^3, Peeble's convention.
meanSampling          =       0.4;

// Binning interval for P(k).
kbinInterval          =      0.01;
modkMax               =      1.00;
muBinNumb             =       100;

// Interval in k^2 for perp k binning of 2D P(k).
perpkInterval         =      0.01;

wfKernel_minAmp       = pow(10., -2.);
convolution_modkmax   =           0.6;

// Total number of mocks. 
CatalogNumber         =        26;

//  Apodise the window fn. to supress the Gibb's phenomenon.
GibbsSkinDepth        =       5.0;

// Random variable generation.    
gsl_rng_env_setup();

gsl_ran_T = gsl_rng_default;
gsl_ran_r = gsl_rng_alloc(gsl_ran_T);
    
// Clipping variables. Currently underclip. 
appliedClippingThreshold  =       -1.0;    

A11Sq                     =        1.0;


// NFW halo generation
haloNumber = 10000;
NFW_conc   =   1.0; 

// Correlation fn's, logarithmic binning in r. 
zerolog  =             log10(0.001);
maxlog   =             log10(2000.0); // hiRes: 10., lowRes: 2000.
logbinsz =             log10( 1.01); // previously 1.4, must be >1.0 otherwise log gives 0. or -ve.
  
nlogbins =  (int) ceil((maxlog - zerolog)/logbinsz);

// linear binning in mu. 
zerolin  =                   0.00;
maxlin   =                   1.00;
linbinsz =                   0.05;

nlinbins =  (int) ceil((maxlin - zerolin)/linbinsz);

printf("\n\nPair counting binning: %d \t %d", nlogbins, nlinbins);

// padVolume(0.0, 0.0, 0.0);

comovDistReshiftCalc();

// VIPERS_SolidAngle = SolidAngleCalc(LowerDecLimit, UpperDecLimit, UpperRAlimit-LowerRAlimit);

// JenkinsCoordinates();

// Checked.
EvaluateGridParameters();

// assign binning interval in k, and calcualate number of bins required. Checked. 
assignbinninginterval(kbinInterval);

// Checked.
prepNGP();

// VIPERS_Binarymask();
// knownGRF_mask();
// randoms_redshiftDistribution();

// sampling: 0.032 for w^2 hi res measurement, 0.008 for lo res, 1.00 for creating mask for P(k)
// homogeneous_rands_window_volCalc(0.95);
load_homogeneous_rands_window(50000000., 1, 0.0001);
// write_homogeneous_rands_window_gridded(50000000., 1, 1.0);

// FullCube();

// CalcCellraDec();

// CalcCellChi();

// apodiseWindowfn();

// apodisedVolume        = CellVolume*SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);

prepBootStrap(n0*n1*n2, Cell_rotatedXvals, Cell_rotatedYvals, Cell_rotatedZvals, 1000.);

// Applied window fn.
Cell_AppliedWindowFn  = &Cell_SurveyLimitsMask[0];

// Checked.
CalcWfCorrections();

prepFFTw(n0, n1, n2);

prepFFTbinning(kbinInterval);

assign2DPkMemory(muBinNumb, kBinNumb);       

// inputHODPk();

inputLinearPk();

// pt2shot = &lightconeShot;

// Number of randoms, load positions, load RR.
// prep_randpairs(20.*30000., 1, 0);

// loadNagoya_rands();

// initialise_angularCorrelation();

// randoms_angular_correlationfn();
// sprintf(filepath, "%s/redshiftDistribution_NagoyaRandoms.dat", root_dir);

// output = fopen(filepath, "w");

// NFWprofile_oneHalo();

// NFWprofile_oneHalo_pairCount(); 

// NFWprofile_oneHalo_xiCalc();

// NFW_profileOneHalo_xi();

// computeCorrelation_fns(2);

// Change seed before calling. 
// randoms_Sphere(30000., 200.);
// randoms_anisoGauss(50000000., 40., 40., 60.);

// wfPkCalc();

// randWindow_pairCount();

// prep_NFWhaloCat(3000000);

// prep_pairwisepdf();

// && anisotropicGaussian_multipoles();

prep_VIPERS_maskMonopole();

print_W2_2D();

for(loopCount=1; loopCount<1; loopCount++){
    printf("\n\n%d", loopCount);
    
    // sprintf(surveyType, "250_fullCube_Gaussian_fogRSD_CellSize_%.2f", CellSize);

    // if(loopCount<10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    // else              sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);
      
    // Choice of redshift from zcos, zpec, zphot, zobs.
    
    // CatalogueInput(filepath);
    
    // zUtilized = &zcos[0];
    
    // for(jj=0; jj<Vipers_Num; jj++){      //    if((zUtilized[jj] > 0.65) && (zUtilized[jj]<0.75)){
	//        fprintf(output, "%e \n", zUtilized[jj]);
    // 
    // }
    
    // My basis, otherwise use Stefano basis. Never use both in conjuction. 
    // CoordinateCalc();
    
    // Convert from ra, dec, z to x, y, z in Stefano's basis. 
    // StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
   
    // fiberCollision_cat(rand_number, rand_redshift, rand_x, rand_y, rand_z);
    
    // fiberCollision_cat(Vipers_Num, zUtilized, xCoor, yCoor, zCoor);

    // correlationfn();

    // angular_correlationfn();

    // randomGeneration();
    
    // loadRand();
    
    // projectVIPERSsystem();

    // VIPERSbasis(CentreRA, CentreDec, xCoor, yCoor, zCoor, Vipers_Num);
    
    // assignAcceptance();

    // Must be run for all 27 mocks in preparation for P(k) calc.
    // ComovingNumberDensityCalc();
    // pt2nz = &interp_nz;

    // splineGaussfilteredW1_W4_nz();
    
    // splineMockAvg_nz();
    
    // pt2nz = &MockAvg_nz;
    
    Gaussianfield();
    
    // deplete_homogeneous();
    
    // underclipDensity(appliedClippingThreshold);
    
    // poissonSample_lnNorm();
    
    // Do not forget call to prep Cat. 
    // HaloCatalogue_NFWprofile(3000000);
    
    // Note: Dangerous non-commutation of sampling and randomise in load_clustered. may generate biased population if not careful. //
    // load_clustered(1, 0.01);

    // NFW_profile_pairCount();

    // Calculation of multipole moments for the halo model catalogue.
    // ukm_calc();
    
    // CleanNGP();
    
    // Checked.
    // NGPCalcCube(xCoor, yCoor, zCoor, Vipers_Num);
    // NGPCalcCube(rand_x, rand_y, rand_z, rand_number);
    
    // ApplyFKPweights(meanSampling);

    cleanFFTbinning();
  
    PkCalc();
    
    // slowDFTcalc();
}

// VIPERS_mask_intCnsrt();

// VIPERS_mask_cnvldpk();

/*
sprintf(filepath, "%s/Data/stacpolly/spherical_tophat.dat", root_dir);

output = fopen(filepath, "w");

for(j=1; j<1000000; j++) fprintf(output, "%e \t %e \n", 0.0001*j, spherical_tophat_pk(0.0001*j));
 
fclose(output);
 
spherical_randDistribution();
*/

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

// free(rand_x);
// free(rand_y);
// free(rand_z);
// free(rand_chi);

// AnisoConvolution();

// Window func. convolution assuming spherical symmetry/averaging.                                                  
// ConvolveSphericalSymm();

// formPkCube();

// theoryHexadecapole();

// Columns embedded in cube. 
// CovarianceMatrix(64, 20, 3);

// CovarianceInverse();

// multipolesRealisation();

// LikelihoodEval();

// minimiseChiSq();

// Calc_betaPosterior();

// Calc_sigmaPosterior();

// Calc_betaSigmaPosterior();

// kaiser_Multipoles();

// kaiser_nonlinearSuppression_Multipoles();

// clipCorrfn();

// TwoDbinnedClippedPk();

// TwoDbinnedRedshiftSpacePk();

// TwoColumnCompareTest();

// Stacpolly run. 
// MPI_Finalize();
  
    printf("\n\n");

    return 0; 
}
