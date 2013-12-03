// Stacpolly run. 
// #define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

#define   AUXfn_header "/disk1/mjw/Aux_functions/header.h"
#include  AUXfn_header

#include "Scripts/header.h"

#include "Scripts/comovDistRedshiftCalc.c"

#include "Scripts/JenkinsRun.c"
#include "Scripts/GridParams.c"

#include "Scripts/assignMemory.c"

#include "Scripts/CoordinateCalcCube.c"

#include "Scripts/ArtificialWf.c"

#include "Scripts/NGPCalc.c"
#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"
#include "Scripts/FFTw_3Dwf.c"
#include "Scripts/FFTw_3Dwf_pad.c"

#include "Scripts/sphericalConvolvePk.c"
#include "Scripts/ConvolvePkAnisoWf.c"

#include "Scripts/ComovingNumberDensityCalc.c"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_legendre.h>

#include "Scripts/InvErrorfn.c"

#include "Scripts/Clipped_zSpace.c"

#include "Scripts/zCubeCreate.c"
#include "Scripts/rollCube.c"

#include "Scripts/freeMemory.c"


int main(int argc, char **argv){

    sprintf(root_dir, "/disk1/mjw/HOD_MockRun");
    sprintf(vipersHOD_dir, "/disk1/mjw/VIPERS_HOD_Mocks");

    // Stacpolly run. 
    // sprintf(root_dir, "/home/mjw/HOD_MockRun");

    // MPI_Init(&argc,&argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &process_number);

    // Co-ordinate limits (Volume) defined by the random catalogue. 

    // lower_xlimit & upper_xlimit
    AxisLimsArray[0][0]       =        0.0;                                               // h^-1 Mpc
    AxisLimsArray[1][0]       =     1000.0;                                               // h^-1 Mpc

    // lower_ylimit & upper_ylimit
    AxisLimsArray[0][1]       =        0.0;                                               // h^-1 Mpc
    AxisLimsArray[1][1]       =     1000.0;                                               // h^-1 Mpc

    // lower_zlimit & upper_zlimit
    AxisLimsArray[0][2]       =        0.0;                                               // h^-1 Mpc
    AxisLimsArray[1][2]       =     1000.0;                                               // h^-1 Mpc
         
    CellSize                  =        4.0;                                               // Cell size, comoving distance, h^-1 Mpc

    // Selection parameters. Volume limited sample between redshift 0.7 and 0.9
    redshiftLowLimit          =       0.78;
    redshiftHiLimit           =       0.82;
    absMagCut                 =      -20.0;

    // Apply Jenkin's scaling to beat aliasing.
    JenkinsScalefactor        =        1.0;

    // FKP P(k) of interest;
    fkpPk                     =      5000.;                                               // [h^-1 Mpc]^3, Peeble's convention.

    // Binning interval for P(k).
    kbinInterval              =      0.007;
    modkMax                   =        0.5;

    // Convolution of P(k) with window fn.
    // InterpK_binNumber      =        400;
    // MuIntegralPrecision    =       9000;
    wfKernelsize              =         10;

    // padded window fn. calculation.
    sidepad                   =          0; 

    appliedClippingThreshold  =        5.0;
    
    linearBias                = sqrt(2.90);
    
    A11                       =        1.0;
    
    // zCubeCreate();

    comovDistReshiftCalc();

    JenkinsCoordinates();

    EvaluateGridParameters();
    
    // assign binning interval in k, and calculate number of bins required. 
    assignbinninginterval(kbinInterval);

    prepNGP();
    
    sprintf(surveyType, "FullCube_Jenkins%.1f", JenkinsScalefactor);
    
    FullCube();
    // EmbeddedCube(50);
    // Gaussian(250.);
    // PencilBeamSurvey(40, 60, 40, 60);
    // Spherical(250.);

    SumOfBoolDensity    = SumDoubleArray(booldensity);
  
    CoordinateCalcCube();  
    
    // rollcube(xCoor, yCoor, zCoor, Vipers_Num);
    
    // TotalSurveyedVolume initialised to zero in header.h
    TotalSurveyedVolume  = SumOfBoolDensity*CellVolume*pow(JenkinsScalefactor, 3.0);

    prepFFTw(n0, n1, n2);
    prepFFTbinning();

    // NGPCalcCube();
    
    // A11                      =   (1./2.6); // Empirical estimate from suppressed Del2k. 
    // clipDensity(appliedClippingThreshold); // Clipping at 5.0
    
    // CalcCorrections();
    
    // PkCalc();

    // wfPkCalc();
    
    // ConvolveTheory();
    // printWindowfuncSlices();
    
    AnisoConvolution();
    
    // FFTw arrays in and out and binning arrays must be freed and reassigned to the padded size before padded window fn. calc.
    // freeFFTw();
    // freeBinning();
    
    // padwfPkCalc(sidepad);
    
    // assign2DPkMemory();
    
    // inputLinearPk();
    
    // formPkCube();

    // clipCorrfn();

    // InvErrorfnTest();

    // Theory2Dpk();
    
    // Observed2Dpk();

    // Stacpolly run. 
    // MPI_Finalize();

return 0; 
}
