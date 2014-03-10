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
#include <gsl/gsl_sf_legendre.h>

#include "Scripts/header.h"

#include "Scripts/comovDistRedshiftCalc.c"

#include "Scripts/JenkinsRun.c"
#include "Scripts/GridParams.c"

#include "Scripts/assignMemory.c"

#include "Scripts/CoordinateCalc.c"
#include "Scripts/CoordinateCalcCube.c"
#include "Scripts/assignAcceptance.c"

#include "Scripts/ArtificialWf.c"

#include "Scripts/NGPCalc.c"
#include "Scripts/BasisChange.c"
#include "Scripts/CalcCellraDec.c"

#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

#include "Scripts/qSortCompare.c"
#include "Scripts/FFTw_3D.c"
#include "Scripts/FFTw_3Dwf.c"
#include "Scripts/FFTw_3Dwf_pad.c"
#include "Scripts/axesWfSlices.c"

#include "Scripts/MeasureWfKernel.c"
// #include "Scripts/sphericalConvolvePk.c"
#include "Scripts/setConvolutionKernels.c"
// #include "Scripts/AnalyticTestConvolution.c"
#include "Scripts/ConvolvePkAnisoWf.c"

#include "Scripts/IntegralConstraintCorrection.c"

#include "Scripts/ComovingNumberDensityCalc.c"

#include "Scripts/InvErrorfn.c"

#include "Scripts/Clipped_zSpace.c"

#include "Scripts/zCubeCreate.c"
#include "Scripts/rollCube.c"

#include "Scripts/freeMemory.c"

#include "/disk1/mjw/EisensteinHu/power.c"


int main(int argc, char **argv){

    sprintf(root_dir,      "/disk1/mjw/HOD_MockRun");
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
    absMagCut                 =     -20.00;

    // Apply Jenkin's scaling to beat aliasing.
    JenkinsScalefactor        =        1.0;

    // FKP P(k) of interest;
    fkpPk                     =      5000.;                                               // [P(k)] = [h^-1 Mpc]^3, Peeble's convention.

    // Binning interval for P(k).
    kbinInterval              =       0.01;
    modkMax                   =        0.8;
    muBinNumb                 =         50;

    // Must be odd. 2n+1.
    xwfKernelsize             =         35;
    ywfKernelsize             =         35;
    zwfKernelsize             =          3;

    gsl_rng_env_setup();

    gsl_ran_T = gsl_rng_default;
    gsl_ran_r = gsl_rng_alloc(gsl_ran_T);
    
    // Translate cube to boot strap slice samples.
    xtranslateDist            =      120.0;
    ytranslateDist            =      120.0;
    
    // Clipping variables. 
    // appliedClippingThreshold  =        5.0;    
    // linearBias                = sqrt(2.90);
    // A11                       =        1.0;
    
    // zCubeCreate();

    comovDistReshiftCalc();

    // JenkinsCoordinates();

    EvaluateGridParameters();
   
    // assign binning interval in k, and calculate number of bins required. 
    assignbinninginterval(kbinInterval);

    prepNGP();
  
    sprintf(filepath, "%s/Data/HODCube/cube_gal_-20.0.dat", root_dir);
    
    FullCube();
    // EmbeddedCube(50);
    // Gaussian(250.);
    // PencilBeamSurvey(25, 55, 25, 55);
    // Spherical(250.);
    // AnisoGauss(20., 30., 40.);
    
    Cell_AppliedWindowFn  = &Cell_SurveyLimitsMask[0];

    SumOfVIPERSbools      = SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2);

    // TotalSurveyedVolume initialised to zero in header.h
    TotalSurveyedVolume   = SumOfVIPERSbools*CellVolume;

    prepFFTw(n0, n1, n2);
    
    prepFFTbinning();

    assign2DPkMemory();                
                                                                                             
    CoordinateCalcCube(filepath);  
    
    pt2nz = &CubeMeanNumberDensity;
    
    inputPK();

    pt2Pk = &splintMatterPk;
    /*
    sprintf(surveyType, "Cube_Jenkins%.1f", JenkinsScalefactor);
    
    // for(ii=0; ii<8; ii++){ 
       // for(jj=0; jj<8; jj++){
         //   printf("\n%d \t %d", ii, jj);
        
            sprintf(surveyType, "Clipped_zPencilBeamCube_Jenkins%.1f_xtrans_%.2f_ytrans_%.2f", JenkinsScalefactor, ii*xtranslateDist, jj*ytranslateDist);
        
            // rollxy(ii*xtranslateDist, jj*ytranslateDist, Vipers_Num);
    
            // rollcube(xCoor, yCoor, zCoor, Vipers_Num);
    
            // Currently all galaxies accepted. 
            assignAcceptanceCube();

            CleanNGP();

            NGPCalcCube();
	    
	        // clipDensity(5.0);

            // A11                      =   (1./2.6); // Empirical estimate from suppressed Del2k. 
                
            CalcWfCorrections();
            
            cleanFFTbinning();
    
            PkCalc();
     //   }
    //}
    */
    // wfPkCalc();
    
    // Window func. convolution assuming spherical symmetry/averaging.
    // ConvolveSphericalSymm();
    
    // Anisotropic window func. and/or anisotropic P(k) convolution calc. 
    // MeasureAnisoWfKernel();
      
    // printWindowfuncSlices();
    
    // AnisoConvolution();

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
    
    
    sprintf(surveyType, "zPencilBeamCube_Jenkins%.1f_xtrans_%.2f_ytrans_%.2f", JenkinsScalefactor, ii*xtranslateDist, jj*ytranslateDist);
    
    MockAvgMultipole(26);
    
    printf("\n\n");
    
return 0; 
}
