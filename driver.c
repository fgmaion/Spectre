//Stacpolly run. 
//#define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#include "Scripts/header.h"

#include "Scripts/comovDistRedshiftCalc.c"

#include "Scripts/JenkinsRun.c"
#include "Scripts/GridParams.c"

#include "Scripts/HODrandGen.c"
#include "Scripts/RandCoorCalc.c"
#include "Scripts/randNGP.c"

#include "Scripts/assignMemory.c"

#include "Scripts/CoordinateCalc.c"
#include "Scripts/NGPCalc.c"
#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"
#include "Scripts/FFTw_3Dwf.c"

#include "Scripts/ConvolvePk.c"

#include "Scripts/ComovingNumberDensityCalc.c"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>

#include "Scripts/Clipped_zSpace.c"

#include "Scripts/freeMemory.c"


int main(int argc, char **argv){

sprintf(root_dir, "/disk1/mjw/HOD_MockRun");
sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_ValueAddedHOD");

// Stacpolly run. 
// sprintf(root_dir, "/home/mjw/HOD_MockRun");

// MPI_Init(&argc,&argv);
// MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
// MPI_Comm_size(MPI_COMM_WORLD, &process_number);

// With orientation of the -z Cartesian axis to the line of sight. 
// lower_xlimit & upper_xlimit
AxisLimsArray[0][0]   =       0.0;                                                  // h^-1 Mpc
AxisLimsArray[1][0]   =     600.0;                                                  // h^-1 Mpc

// lower_ylimit & upper_ylimit
AxisLimsArray[0][1]   =       0.0;                                                  // h^-1 Mpc
AxisLimsArray[1][1]   =     600.0;                                                  // h^-1 Mpc

// lower_zlimit & upper_zlimit
AxisLimsArray[0][2]   =       0.0;                                                  // h^-1 Mpc
AxisLimsArray[1][2]   =    3000.0;                                                  // h^-1 Mpc

CellSize              =       3.0;                                                  // Cell size, comoving distance, h^-1 Mpc

// Selection parameters. Volume limited sample between redshift 0.7 and 0.9
redshiftLowLimit      =       0.7;
redshiftHiLimit       =       0.9;
absMagCut             =     -20.0;

// Apply Jenkins contraction to beat aliasing. 
JenkinsScalefactor    =       1.0;

// FKP P(k) of interest;
fkpPk                 =     5000.;                                                  // [h^-1 Mpc]^3, Peeble's convention.

// Binning interval for P(k).
kbinInterval          =     0.005;
modkMax               =       0.5;

// Convolution of P(k) with window fn.
// InterpK_binNumber  =       300;
// MuIntegralPrecision=      9000;

// Must be odd. 2n+1
wfKernelsize          =        9;

// Begin fn. calling. 

VIPERS_SolidAngle     = SolidAngleCalc(-5.4, -4.2, 9.);

comovDistReshiftCalc();

JenkinsCoordinates();

EvaluateGridParameters();

// assign binning interval in k, and calcualte number of bins required. 
assignbinninginterval(kbinInterval);

prepNGP();

// Generate randoms for given grid, will read file if exists. End product: BoolDensity array. 
randNGP(); 

SumOfBoolDensity    = SumDoubleArray(booldensity);

//Initialsed to zero in header.h
TotalSurveyedVolume = SumOfBoolDensity*CellVolume;

prepFFTw(n0, n1, n2);
prepFFTbinning();


for(loopCount=1; loopCount<15; loopCount++){
    if(loopCount < 10)  sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else                sprintf(filepath, "%s/mocks_W1_v1.2/mock_W1_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CoordinateCalc(filepath);
  
    // Orientation of the -z Cartesian axis to the line of sight direction.
    // centerRA in deg, centerDec in deg, xCoors[], yCoors[], zCoors[], len, rotation angle in deg.
    // Note: p hat must lie in the xy plane by assumption. -z hat assumed to be rotated to VIPERS los direction. 
   
    redshiftSpaceRotation(34.5, -4.79, xCoor, yCoor, zCoor, Vipers_Num, 180.0 + 94.79);
    
    printf("\n\nAfter orientation of -z Cartesian axis to LOS ...");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));

    // Check inverse rotation.
    redshiftSpaceRotation(34.5, -4.79, xCoor, yCoor, zCoor, Vipers_Num, -(180.0 + 94.79));
    
    printf("\n\nAfter inverse rotation ...");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
  
    // ComovingNumberDensityCalc();

    // ApplyFKPweights();

    // printf("\n\nHighest weighted overdensity estimate: %e",  DoubleArrayMax(FKPweights, n0*n1*n2));
    // printf("\nLowest weighted overdensity estimate:  %e",    DoubleArrayMin(FKPweights, n0*n1*n2));

    NGPCalc();
  
    // CalcCorrections();
  
    // if(loopCount < 10)  sprintf(filepath, "%s/midK_Del2k_HODMocks_00%d.txt", root_dir, loopCount);
    // else                sprintf(filepath, "%s/midK_Del2k_HODMocks_0%d.txt",  root_dir, loopCount);
  
    // PkCalc(filepath);
}

// sprintf(filepath, "%s/midK_WindowFuncPk_HODMocks.dat", root_dir);
// wfPkCalc(filepath);

// printf("\nFreeing allocated memory.");
// freeHOD();
// freeNGP();

// ConvolveTheory();

// assign2DPkMemory();

// inputLinearPk();

// formPkCube();

// clipCorrfn();

// TwoDbinnedClippedPk();

// TwoDbinnedRedshiftSpacePk();

// sprintf(filepath, "%s/Scripts/Normal_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, PkCube);

// sprintf(filepath, "%s/Scripts/zSuppressed_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, clippedPk);

// sprintf(filepath, "%s/Scripts/zDistorted_PkTheory_KaiserLorentzian_0.6mu0.8.dat", root_dir);
// BinnedPkForMuInterval(0.6, 0.8, filepath, clippedPk);

// freeFFTw();
// freePk();

// free_sdltInterp();

// freeConvolutionMemory();

// freeClipped();  
// free_linear();
// free2dPk();

// TwoColumnCompareTest();

// Stacpolly run. 
// MPI_Finalize();

return 0; 
}
