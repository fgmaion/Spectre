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

#include "Scripts/NGPCalc.c"
#include "Scripts/BasisChange.c"
#include "Scripts/CalcCellraDec.c"
#include "Scripts/ApodiseWindowFunc.c"

#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

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

#include "Scripts/MockAvg2D_Pk.c"
#include "Scripts/MultipoleCovariance.c"

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
AxisLimsArray[0][0]   =    1673.0;                                                  // h^-1 Mpc
AxisLimsArray[1][0]   =    2200.0;                                                  // h^-1 Mpc

// lower_ylimit & upper_ylimit
AxisLimsArray[0][1]   =    -185.0;                                                  // h^-1 Mpc
AxisLimsArray[1][1]   =     185.0;                                                  // h^-1 Mpc

// lower_zlimit & upper_zlimit
AxisLimsArray[0][2]   =     -25.0;                                                  // h^-1 Mpc
AxisLimsArray[1][2]   =      55.0;                                                  // h^-1 Mpc

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
CentreDec             =    -3.525;

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
CellSize              =       1.0;                                                  

// Selection parameters.
absMagCut             =      20.0;
redshiftLowLimit      =      0.70;
redshiftHiLimit       =      0.90;

// Comoving number density, n(z), measurement. 
zBinWidth             =      0.05; 
chiBinWidth           =     10.00;
nzSigma               =     100.0;

// Apply Jenkins contraction to beat aliasing. 
JenkinsScalefactor    =       1.0;

// FKP P(k) of interest;
fkpPk                 =     5000.;                                                  // [h^-1 Mpc]^3, Peeble's convention.

// Binning interval for P(k).
kbinInterval          =      0.01;
modkMax               =       0.6;
muBinNumb             =       100;

// Interval in k^2 for perp k binning of 2D P(k).
perpkInterval         =      0.02;

// Must be odd. 2n+1.
xwfKernelsize         =        15;
ywfKernelsize         =         9;
zwfKernelsize         =         9;

// Total number of mocks. 
CatalogNumber         =        26;

//  Apodise the window fn. to supress the Gibb's phenomenon.
GibbsSkinDepth        =       5.0;

// Checked.
comovDistReshiftCalc();

// Checked.
VIPERS_SolidAngle     = SolidAngleCalc(LowerDecLimit, UpperDecLimit, UpperRAlimit-LowerRAlimit);

sprintf(surveyType, "VIPERSparent_GaussSmoothMockAvgNz_%.1f_SkinDepth_%.1f", nzSigma, GibbsSkinDepth);

// JenkinsCoordinates();

// Checked.
EvaluateGridParameters();

// Checked.
prepNGP();

// Checked.
CalcCellraDec();

// Applied window fn.
Cell_AppliedWindowFn  = &Cell_SurveyLimitsMask[0];

// Checked.
SumOfVIPERSbools      = SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);

// Initialsed to zero in header.h
TotalSurveyedVolume   = SumOfVIPERSbools*CellVolume;

// apodiseWindowfn();

// apodisedVolume        = CellVolume*SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);

// assign binning interval in k, and calcualate number of bins required. Checked. 
assignbinninginterval(kbinInterval);

prepFFTw(n0, n1, n2);

prepFFTbinning();

// assign2DPkMemory();

for(loopCount=1; loopCount<2; loopCount++){
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CatalogueInput(filepath);
      
    // Choice of redshift from zcos, zpec, zphot, zobs.
    zUtilized = &zcos[0];

    CoordinateCalc();

    VIPERSbasis(CentreRA, CentreDec, xCoor, yCoor, zCoor, Vipers_Num);
    
    assignAcceptance();

    // Must be run for all 27 mocks in preparation for P(k) calc.
    // ComovingNumberDensityCalc();

    // splineGaussfilteredW1_W4_nz();

    // splineMockAvg_nz();
    
    // pt2nz = &interp_nz;

    // pt2nz = &MockAvg_nz;

    // ApplyFKPweights();

    // CleanNGP();

    // Checked.
    // NGPCalc();
  
    // Checked.
    // CalcWfCorrections();

    // cleanFFTbinning();
  
    // PkCalc();
}

// MockAvg2Dpk(14);

// MockAverageComovingdensity();

// fitnz();

// freeHOD();

// freeNGP();

// wfPkCalc();

// printWindowfuncSlices();

// analyticConvTest();

// AnisoConvolution();

// Window func. convolution assuming spherical symmetry/averaging.                                                  
// ConvolveSphericalSymm();

// inputPK();

// formPkCube();

// theoryQuadrupole();

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

    printf("\n\n");

    return 0; 
}
