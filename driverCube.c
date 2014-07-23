#include <stdbool.h>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_bessel.h>

// Stacpolly run. 
// #define AUXfn_DIR "/home/mjw/Aux_functions/header.h"

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
#include "Scripts/CoordinateCalcCube.c"
#include "Scripts/assignAcceptance.c"

#include "Scripts/ArtificialWf.c"

#include "Scripts/NGPCalc.c"
#include "Scripts/CloudInCell.c"
#include "Scripts/BasisChange.c"
#include "Scripts/CalcCellraDec.c"

#include "Scripts/FKPweights.c"
#include "Scripts/Windowfn_PkCorrections.c"

#include "Scripts/KaiserMultipoles.c"
// #include "Scripts/KaiserGaussMultipoles.c"
#include "Scripts/KaiserLorentzMultipoles.c"

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

#include "Scripts/AgeOftheUniverse.c"
#include "Scripts/linearGrowthRate.c"
#include "Scripts/growthfactor_derivative.c"

#include "Scripts/InvErrorfn.c"
#include "Scripts/MatrixInverse.c"

#include "Scripts/Clipped_zSpace.c"
#include "Scripts/SpectralDistortion.c"

#include "Scripts/zCubeCreate.c"
#include "Scripts/wedgeMockCreate.c"
#include "Scripts/rollCube.c"

#include "Scripts/MultipoleCovariance.c"
#include "Scripts/Multipoles_EigenVecsCovariance.c"
#include "Scripts/MultipoleCovariance_Inverse.c"
#include "Scripts/MultipolesRealisation_MultiVariateGauss.c"
#include "Scripts/Multipole_minimiseChiSq.c"
#include "Scripts/Multipole_MarginalisedPosteriors.c"
#include "Scripts/Multipole_2DPosteriors.c"

#include "Scripts/BootStrap.c"
#include "Scripts/freeMemory.c"

#include "/disk1/mjw/EisensteinHu/power.c"

#include "Scripts/FFT_log.h"
#include "Scripts/FFT_log.c"


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
    AxisLimsArray[0][0]       =        0.0;                                 // h^-1 Mpc
    AxisLimsArray[1][0]       =      1000.;                                 // h^-1 Mpc

    // lower_ylimit & upper_ylimit
    AxisLimsArray[0][1]       =        0.0;                                 // h^-1 Mpc
    AxisLimsArray[1][1]       =      1000.;                                 // h^-1 Mpc

    // lower_zlimit & upper_zlimit
    AxisLimsArray[0][2]       =        0.0;                                 // h^-1 Mpc
    AxisLimsArray[1][2]       =      1000.;                                 // h^-1 Mpc
                
    CellSize                  =        4.0;                                 // Cell size, comoving distance, h^-1 Mpc

    // Selection parameters. Mag 20.0 galaxies at redshift 0.8;
    redshiftLowLimit          =       0.795;
    redshiftHiLimit           =       0.805;
    absMagCut                 =      -20.00;

    // Apply Jenkin's scaling to beat aliasing.
    JenkinsScalefactor        =        1.0;

    // FKP P(k) of interest.
    fkpPk                     =      5000.;                                               // [P(k)] = [h^-1 Mpc]^3, Peeble's convention.

    // Binning interval for P(k).
    kbinInterval              =       0.01;
    modkMax                   =        1.2;
    muBinNumb                 =         50;

    gsl_rng_env_setup();

    gsl_ran_T                 = gsl_rng_default;
    gsl_ran_r                 = gsl_rng_alloc(gsl_ran_T);
    
    // Non-linear RSD
    beta                      =       0.542;                  // beta = 0.542
    velDispersion             =        5.00;                   // units of h^-1 Mpc rather than 300 km s^-1
    A11Sq                     =       0.567;
    
    // Priors on the model params.
    min_beta                  =      0.45;
    max_beta                  =      0.55;

    min_velDisperse           =       2.0;
    max_velDisperse           =       4.0;
 
    // down weight theory. 
    min_A11Sq                 =       0.9999999;
    max_A11Sq                 =       1.0000001;

    // Resolution of the Likelihood evaluation [voxel number].
     Res                      =        30;
    dRes                      =      30.0;
    
    ChiSqEval_kmin            =      0.02;
    ChiSqEval_kmax            =      0.80; 
    
    // Fit solely the monopole (1) or both mono and Quad (2).
    hiMultipoleOrder          =         2;
    
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
    
    linearBias                =   1.495903;
    
    // Clipping variables. 
    appliedClippingThreshold  =        1.5;    
    // linearBias             = sqrt(2.90);
    
    // zCubeCreate();
    
    // wedgeMockCreate(400., 150., 600., 150., 700., 850., 300., 850., 475., 525.);
    /*
    comovDistReshiftCalc();

    // JenkinsCoordinates();
    
    EvaluateGridParameters();
   
    // assign binning interval in k, and calculate number of bins required. 
    assignbinninginterval(kbinInterval);

    prepNGP();
    
    // No artificially applied window fn. 
    
    FullCube();
    // EmbeddedCube(50);
    // Gaussian(250.);
    // PencilBeamSurvey(25, 55, 25, 55);
    // Spherical(250.);
    // AnisoGauss(20., 30., 40.);
    
    Cell_AppliedWindowFn  = &Cell_SurveyLimitsMask[0];
    
    CalcCellraDec();

    // TotalSurveyedVolume initialised to zero in header.h
    TotalSurveyedVolume   = SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2)*CellVolume;
    
    prepFFTw(n0, n1, n2);
    
    prepFFTbinning();

    assign2DPkMemory();                
    
    // Choice of real or redshift space, zcube or cube. 
    sprintf(filepath, "%s/Data/HODCube/zcube_zvel_gal_-20.0.dat", root_dir);       
                                                                                 
    CoordinateCalcCube(filepath);  
    */
    inputHODPk();
    
    // inputLinearPk();
    
    // setKaiserRSD();
    
    // setLorentzianRSD();
    
    // setGaussianRSD();

    pt2shot = &CubeShot;
    
    // Mean number density is the same when boot strapping. 
    pt2nz   = &CubeMeanNumberDensity;
    
    // Currently all galaxies accepted. 
    /*
    assignAcceptanceCube();
    
    CleanNGP();
    
    // sprintf(surveyType, "zCube_xvel_clipThreshold_%.1e_fullCube", appliedClippingThreshold);
    
    NGPCalcCube();
    
    // Overwrites measured density field with a Gaussian random field, with given P(k).
    Gaussianfield();
    
    clipDensity(appliedClippingThreshold);
       
    prepBootStrap(n0*n1*n2, Cell_rotatedXvals, Cell_rotatedYvals, Cell_rotatedZvals, 1000.);
       
    // for(loopCount=0; loopCount<100; loopCount++){
    
    sprintf(surveyType, "GaussianCube_zvel_MonoQuadPk_powLawXi_KaiserRSD_clipThreshold_%.1e_fullCube", appliedClippingThreshold);
              
    // BootStrapGen(n0*n1*n2, Cell_rotatedXvals, Cell_rotatedYvals, Cell_rotatedZvals, 1000.);
                
    CalcWfCorrections();
            
    cleanFFTbinning();
    
    // ** Shot noise correction is switched off **
    PkCalc();
    // }
    */
    // wfPkCalc();
    
    // kaiser_nonlinearSuppression_Multipoles();
    
    // LikelihoodMemory();
    
    // CovarianceMatrix(100);
        
    // CovarianceEigenVecs();
    
    // ClippingModelling();                                                            // Assigns memory for Clipping prediction. 
    
    // Mono_xi();
    
    // formPkCube();

    // clipCorrfn();
    
    // multipolesRealisation();
    
    // minimiseChiSq();
    
    // Calc_betaPosterior();

    // Calc_sigmaPosterior();

    // Calc_betaSigmaPosterior();
    
    // ConfidenceLimits_2D();
    
    // assign2DPkMemory();

    // InvErrorfnTest();

    // Theory2Dpk();
    
    // Observed2Dpk();

    // Stacpolly run. 
    // MPI_Finalize();
    
    // sprintf(surveyType, "zPencilBeamCube_Jenkins%.1f_xtrans_%.2f_ytrans_%.2f", JenkinsScalefactor, ii*xtranslateDist, jj*ytranslateDist);
    
    // MockAvgMultipole(26);
    
    // growthfactor_derivative();
    
    // Eisenstein & Hu P(k)
    // TFmdm_set_cosm(0.3, 0.047337, 0.0, 0, 0.7, 0.65, 0.0);
    
    FFTlogRes   = 4096;
    
    // Correlation fns. 
    mono   = malloc(FFTlogRes*sizeof(double));
    quad   = malloc(FFTlogRes*sizeof(double));
     hex   = malloc(FFTlogRes*sizeof(double));
    
    // Power spectra
    monop  = malloc(FFTlogRes*sizeof(double));
    quadp  = malloc(FFTlogRes*sizeof(double));
     hexp  = malloc(FFTlogRes*sizeof(double));
    
    double*     r   = malloc(FFTlogRes*sizeof(double));
    double* kVals   = malloc(FFTlogRes*sizeof(double));
    
    cmono   = malloc(FFTlogRes*sizeof(double)); 
    cmonop  = malloc(FFTlogRes*sizeof(double)); 
    
    double* power   = malloc(FFTlogRes*sizeof(double)); 
    double*    Xi   = malloc(FFTlogRes*sizeof(double)); 
    
    /*
    // Correlation functions given input P(k).
    xi_mu(0, r, mono, beta, velDispersion, kVals, FFTlogRes);
    
    xi_mu(2, r, quad, beta, velDispersion, kVals, FFTlogRes);
    
    xi_mu(4, r,  hex, beta, velDispersion, kVals, FFTlogRes);
    
    sprintf(filepath, "%s/Data/SpectralDistortion/TestMonopole_PowerLaw_Multipoles.dat", root_dir);

    output = fopen(filepath, "w");

    for(i=0; i<FFTlogRes; i++) fprintf(output, "%e \t %e \t %e \t %e \n", r[i], mono[i]/Pk_powerlaw_xi(r[i], 5., 1.8), quad[i]/Pk_powerlaw_xi(r[i], 5., 1.8), hex[i]/Pk_powerlaw_xi(r[i], 5., 1.8));

    fclose(output);
    
    // Power spectra given input correlation functions. 
    pk_mu(0, r, mono, beta, velDispersion, kVals, FFTlogRes);
    
    pk_mu(2, r, quad, beta, velDispersion, kVals, FFTlogRes);
    
    pk_mu(4, r,  hex, beta, velDispersion, kVals, FFTlogRes);
    
    sprintf(filepath, "%s/Data/SpectralDistortion/TestMonopole_PowerLaw_pkMultipoles.dat", root_dir);

    output = fopen(filepath, "w");

    for(i=0; i<FFTlogRes; i++) fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \n", kVals[i], mono[i], quad[i], hex[i], Pk_powerlaw(kVals[i], 5., 1.8)*kaiser_multipole(kVals[i], beta, 0), Pk_powerlaw(kVals[i], 5., 1.8)*kaiser_multipole(kVals[i], beta, 2), Pk_powerlaw(kVals[i], 5., 1.8)*kaiser_multipole(kVals[i], beta, 4));

    fclose(output);
    */
    /*
    double u0        = 0.5175669;
    
    double variance  =  4.199712*pow(10., 9.);
    
    // Spectral distortion calculation. 
    // Correlation functions given input P(k).
    xi_mu(0, r, mono, beta, velDispersion, kVals, FFTlogRes, power, Xi);
    
    xi_mu(2, r, quad, beta, velDispersion, kVals, FFTlogRes, power, Xi);
    
    xi_mu(4, r,  hex, beta, velDispersion, kVals, FFTlogRes, power, Xi);
    
    // Clipped monopole calculation to 2nd order. 
    
    for(i=0; i<FFTlogRes; i++){
        cmono[i]  = mono[i];
    
        // cmono[i] *= 0.25*pow(1. + gsl_sf_erf(u0), 2.);
            
        // cmono[i] += (C_n(u0, 1)/variance)*((1./64.)*pow(8.*mono[i] - 4.*quad[i] + 3.*hex[i], 2.) - pow(quad[i], 2.)/20. + mono[i]*(quad[i] - 3.*hex[i]/4.) + 3.*quad[i]*hex[i]/8. - (17./576.)*hex[i]*hex[i]);
    }
    
    rVals    = malloc(FFTlogRes*sizeof(float));
    fcmono   = malloc(FFTlogRes*sizeof(float));
    fcmono2D = malloc(FFTlogRes*sizeof(float));
    
    for(i=0; i<FFTlogRes; i++){
         rVals[i]   = (float)     r[i];
        fcmono[i]   = (float)     Pk_powerlaw_xi(r[i], 5., 1.8); //cmono[i];
    } 
    
    spline(rVals, fcmono, FFTlogRes-1, 1.0e31, 1.0e31, fcmono2D);
        
    sprintf(filepath, "%s/Data/SpectralDistortion/ClippedMonopole_Order2_xi.dat", root_dir);

    output = fopen(filepath, "w");

    for(i=0; i<FFTlogRes; i++) fprintf(output, "%e \t %e \t %e \t %e \t %e\n", r[i], mono[i]/Pk_powerlaw_xi(r[i], 5., 1.8), quad[i]/Pk_powerlaw_xi(r[i], 5., 1.8), hex[i]/Pk_powerlaw_xi(r[i], 5., 1.8), cmono[i]/Pk_powerlaw_xi(r[i], 5., 1.8));

    fclose(output);
    */
    // Power spectra
    monop           = malloc(FFTlogRes*sizeof(double));
    quadp           = malloc(FFTlogRes*sizeof(double));
     hexp           = malloc(FFTlogRes*sizeof(double));
    
    cmonop          = malloc(FFTlogRes*sizeof(double)); 
    
    pk_mu(0, r,  monop, beta, velDispersion, kVals, FFTlogRes);
    
    // sprintf(filepath, "%s/Data/SpectralDistortion/ClippedMonopole_Order2_pk.dat", root_dir);

    // output = fopen(filepath, "w");

    // for(i=0; i<FFTlogRes; i++) fprintf(output, "%e \t %e \t %e \n", kVals[i], monop[i], Pk_powerlaw(kVals[i], 5., 1.8));

    // fclose(output);
    
    printf("\n\n");
    
    return 0;
}


int pk_mu(int mu, double* r, double* multi, double beta, double velDispersion, double* kVals, int FFTlogRes){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    fc = FFTLog_init(FFTlogRes, pow(10., -2.), pow(10., 6.), 0.0, mu + 0.5);
        
    FFTLog_setInput(fc, kVals, r, beta, velDispersion);
    
    // for(i=0; i<fc->N; i++)  fc->xi[i][0]  = pow(2.*pi*r[i], 3./2.)*Pk_powerlaw_xi(r[i], 5., 1.8);
    
    print_fcprop(fc);
    
    FFTLog(fc, fc->xi_forwardplan, fc->xi_backwardplan); 
    
    // reverse array
    for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->pk[i][0],fc->pk[fc->N-i-1][0]);
  
    for(i=0; i<fc->N;   i++) fc->pk[i][0] = (fc->pk[i][0]/(double)fc->N);

    for(i=0; i<fc->N;   i++) multi[i] = pow(-1., mu/2)*fc->pk[i][0]*pow(kVals[i], -1.5);

    // FFTLog_free(fc);    

    return 0;
}


double splint_cmono(double r){
        float Interim;
    
        splint(rVals, fcmono, fcmono2D, fc->N, (float) r, &Interim);
    
        return (double) Interim;
    }


int xi_mu(int mu, double* r, double* multi, double beta, double velDispersion, double* kVals, int FFTlogRes, double* Pk, double* Xi){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    
    fc = FFTLog_init(FFTlogRes, pow(10., -8.), pow(10., 8.), 0.0, mu + 0.5);
    
    FFTLog_setInput(fc, kVals, r, beta, velDispersion);
    
    FFTLog(fc, fc->pk_forwardplan, fc->pk_backwardplan); 
    
    // reverse array
    for(i=0; i<fc->N/2; i++) FFTLog_SWAP(fc->xi[i][0],fc->xi[fc->N-i-1][0]);
  
    for(i=0; i<fc->N;   i++) fc->xi[i][0] = (fc->xi[i][0]/(double)fc->N);

    for(i=0; i<fc->N;   i++) multi[i] = pow(-1., mu/2)*fc->xi[i][0]*pow(r[i], -1.5);
    
    // for(i=0; i<fc->N;   i++)  printf("%e \t %e \t %e \t %e \t %e \t %e \n", kVals[i], (*pt2Pk)(kVals[i]), kaiserLorentz_multipole(kVals[i]*velDispersion, beta, 0), sqrt(pow(kVals[i], 3.)/(8.*pow(pi, 3.)))*(*pt2Pk)(kVals[i])*kaiserLorentz_multipole(kVals[i]*velDispersion, beta, 0), fc->pk[i][0], multi[i]);

    // FFTLog_free(fc);    
    
    return 0;
}
