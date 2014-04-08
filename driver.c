// Stacpolly run. 
// #define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#define   AUXfn_funcs  "/disk1/mjw/Aux_functions/Aux_functions.c"
#include  AUXfn_funcs

#include <stdbool.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_legendre.h>

#include "Scripts/header.h"

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
// #include "Scripts/KaiserGaussMultipoles.c"
#include "Scripts/KaiserLorentzMultipoles.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"
#include "Scripts/FFTw_3Dwf.c"
#include "Scripts/axesWfSlices.c"

#include "Scripts/MeasureWfKernel.c"
// #include "Scripts/sphericalConvolvePk.c"
#include "Scripts/ConvolvePkAnisoWf.c"
#include "Scripts/setConvolutionKernels.c"
// #include "Scripts/AnalyticTestConvolution.c"

#include "Scripts/IntegralConstraintCorrection.c"

#include "Scripts/ComovingNumberDensityCalc.c"
#include "Scripts/MockAvgComovingDensity.c"

#include "Scripts/InvErrorfn.c"

#include "Scripts/Clipped_zSpace.c"

#include "Scripts/PowellsRoutine.c"
#include "Scripts/Powells_mnbrak.c"
#include "Scripts/Powells_linmin.c"
#include "Scripts/Powells_Brent.c"

#include "Scripts/fitting_nz.c"
#include "Scripts/randGen.c"


#include "Scripts/slowDFT.c"
/*
// #include "Scripts/MockAvgMultipole.c"
#include "Scripts/MultipoleCovariance.c"
#include "Scripts/MultipoleCovariance_Inverse.c"
#include "Scripts/MultipolesRealisation_MultiVariateGauss.c"
#include "Scripts/MarginalisedPosteriors.c"
#include "Scripts/2DPosteriors.c"
#include "Scripts/LikelihoodEval.c"
*/
#include "Scripts/freeMemory.c"


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
AxisLimsArray[1][0]   =     800.0;                                                     // h^-1 Mpc

AxisLimsArray[0][1]   =       0.0;                                                     // h^-1 Mpc
AxisLimsArray[1][1]   =     800.0;                                                     // h^-1 Mpc

AxisLimsArray[0][2]   =       0.0;                                                     // h^-1 Mpc
AxisLimsArray[1][2]   =     800.0;                                                     // h^-1 Mpc

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
// CellSize           =    7.8125;     
// CellSize           =     5.208;

// Selection parameters.
absMagCut             =     22.00;
redshiftLowLimit      =      0.60;
redshiftHiLimit       =      0.90;

// Non-linear RSD
velDispersion         =      3.00;                  // units of h^-1 Mpc rather than 300 km s^-1
beta                  =      0.35;                  // 2dF measurement, beta = 0.43
A11Sq                 =      2.58;

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
chiBinWidth           =     15.00;
nzSigma               =      50.0;

// Apply Jenkins contraction to beat aliasing. 
JenkinsScalefactor    =       1.0;

// FKP P(k) of interest;
fkpPk                 =     5000.;                                                  // [h^-1 Mpc]^3, Peeble's convention.

// Binning interval for P(k).
kbinInterval          =      0.01;
modkMax               =      0.80;
muBinNumb             =       100;

// Interval in k^2 for perp k binning of 2D P(k).
perpkInterval         =      0.02;

// Must be odd. 2n+1.
// xwfKernelsize      =        15;
// ywfKernelsize      =         9;
// zwfKernelsize      =         9;

// xwfKernelsize         =       165;
// ywfKernelsize         =        95;
// zwfKernelsize         =        50;

// Total number of mocks. 
CatalogNumber         =        26;

//  Apodise the window fn. to supress the Gibb's phenomenon.
GibbsSkinDepth        =       5.0;

// Random variable generation.    
gsl_rng_env_setup();

gsl_ran_T = gsl_rng_default;
gsl_ran_r = gsl_rng_alloc(gsl_ran_T);

padVolume(0.0, 0.0, 0.0);

// Checked.
comovDistReshiftCalc();

// Checked.
VIPERS_SolidAngle     = SolidAngleCalc(LowerDecLimit, UpperDecLimit, UpperRAlimit-LowerRAlimit);

sprintf(surveyType, "VIPERSparent_StefanoBasis_Cic_GaussSmoothNz_%.1f_SkinDepth_%.1f_VolLim_%.2f_Mesh_%.2f", nzSigma, GibbsSkinDepth, absMagCut, CellSize);

// JenkinsCoordinates();

// Checked.
EvaluateGridParameters();

// Checked.
prepNGP();

// Checked.
// CalcCellraDec();
// StefanoMask();

// apodiseWindowfn();

// apodisedVolume        = CellVolume*SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);

// assign binning interval in k, and calcualate number of bins required. Checked. 
assignbinninginterval(kbinInterval);

prepFFTw(n0, n1, n2);

prepFFTbinning();

// assign2DPkMemory();

sprintf(theoryPk_flag, "HOD_-20.0");
pt2Pk = &splintHODpk;
//pt2Pk = &splintLinearPk;

inputHODPk();
// inputLinearPk();
/*
for(loopCount=1; loopCount<27; loopCount++){
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);
      
    // Choice of redshift from zcos, zpec, zphot, zobs.
    
    CatalogueInput(filepath);
    
    zUtilized = &zcos[0];

    // CoordinateCalc();
    
    StefanoBasis(Vipers_Num, ra, dec, rDist, xCoor, yCoor, zCoor);
    
    StefanoRotated(Vipers_Num, CentreRA, CentreDec, xCoor, yCoor, zCoor);
    
    // randomGeneration();
    
    loadRand();
    
    // Applied window fn.
    Cell_AppliedWindowFn  = &Cell_SurveyLimitsMask[0];

    // Initialsed to zero in header.h
    TotalSurveyedVolume   = SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2)*CellVolume;
    
    printf("\n\nSurveyed volume:  %e", TotalSurveyedVolume/TotalVolume);
    
    // projectVIPERSsystem();

    // VIPERSbasis(CentreRA, CentreDec, xCoor, yCoor, zCoor, Vipers_Num);
    
    assignAcceptance();

    // Must be run for all 27 mocks in preparation for P(k) calc.
    // ComovingNumberDensityCalc();

    // splineGaussfilteredW1_W4_nz();
    splineStefano_nz();
    
    // splineMockAvg_nz();
    pt2nz    =  &interp_Stefano_nz;
    
    // pt2nz = &interp_nz;

    // pt2nz = &MockAvg_nz;
    
    // ApplyFKPweights();

    CleanNGP();

    // Checked.
    NGPCalc();
   
    // Checked.
    CalcWfCorrections();

    cleanFFTbinning();
  
    PkCalc();
    
    // slowDFTcalc();
}
*/
// MockAvg2Dpk(14);

// MockAverageComovingdensity();

// fitnz();

// freeHOD();

// freeNGP();

// printHODMatterPk();

// wfPkCalc();

// analyticConvTest();

// AnisoConvolution();

// Window func. convolution assuming spherical symmetry/averaging.                                                  
// ConvolveSphericalSymm();

// formPkCube();

// theoryHexadecapole();

// CovarianceMatrix();

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

// sprintf(filepath, "%s/Scripts/Normal_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, PkCube);

// sprintf(filepath, "%s/Scripts/zSuppressed_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, clippedPk);

// sprintf(filepath, "%s/Scripts/zDistorted_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, clippedPk);

// TwoColumnCompareTest();

// Stacpolly run. 
// MPI_Finalize();
    
    
    // Inverse of M using LAPACK.
/*
    double** M;
    
    M = (double **) malloc(3*sizeof(*M));
    
    for(j=0; j<3; j++)  M[j] = malloc(3*sizeof(**M));
    
    for(k=0; k<3; k++){
        for(j=0; j<3; j++){
            M[j][k] = 0.0;
        }
    }
    
    M[0][0] = 4.0;
    M[0][1] = 12.0;
    M[0][2] = -16.0;

    M[1][0] = 12.0;
    M[1][1] = 37.0;
    M[1][2] = -43.0;
    
    M[2][0] = -16.0;
    M[2][1] = -43.0;
    M[2][2] =  98.0;

    int    order  = 3;
    int    order2 = 9;

    // length of N+1
    double pivotArray[4];
    
    int    errorHandler; 

    double lapackWorkspace[9];

    double flattenedM[9];

    int count = 0;

    // Note: FORTRAN Ordering. M = {{a, b, c},
    //                              {d, e, f},
    //                              {g, h, i}};

    // Assigned via: M[0][0] = a, M[0][1] = b, M[0][2] = c, M[1][0] = d, etc.

    for(j=0; j<3; j++){
      for(k=0; k<3; k++){
	flattenedM[count] = M[j][k];
	count            += 1;
      }
    }
    */
    // dgetrf_(&order, &order, flattenedM, &order, pivotArray, &errorHandler);
    
    // dgetri_(&order, flattenedM, &order, pivotArray, lapackWorkspace, &order2, &errorHandler);
    /*
    printf("\n %e \t %e \t %e", flattenedM[0], flattenedM[1], flattenedM[2]);
    printf("\n %e \t %e \t %e", flattenedM[3], flattenedM[4], flattenedM[5]);
    printf("\n %e \t %e \t %e", flattenedM[6], flattenedM[7], flattenedM[8]);
    */  
    
    // Cholesky decomposition of a Matrix, with LAPACK. 
    /*
    int INFO, size;
    char UPLO;
    
    UPLO ='L';

    size = 3;

    // Cholesky decomposition. Upper triangular matrix, lower triangular
    // left as original matrix. 

    dpotrf_(&UPLO, &size, flattenedM, &size, &INFO);
    
    count = 0;

    for(j=0; j<3; j++){
      printf("\n");
      for(k=0; k<3; k++){
	printf("%e \t", flattenedM[count]);
	count += 1;
      }
    }
    */
    printf("\n\n");

    return 0; 
}
