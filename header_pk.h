// Path to relevant directories. 
char         vipersHOD_dir[200];
char   maskmultipoles_path[200];
char covariance_mocks_path[200];

// Misc. //
int    data_mock_flag;

// r2c or c2c arrays/ 
int     fft_size;

double* overdensity;
double* smooth_overdensity;

fftw_complex* H_k;

//-- embedding volume --// 
double    TotalVolume         = 0.0;
double    TotalSurveyedVolume = 0.0;




double    stefano_trans_x;  // translation parameters for Stefano's co-ordinates.
double    stefano_trans_y;
double    stefano_trans_z;



//-- VIPERS --//
int               max_gals;  // max. number of galaxies (lines) in any mock in covariance calc, i.e. in files in mocks directory.
int           accepted = 0;
int             Vipers_Num;

double        loChi;
double        hiChi;

double        UpperRAlimit;
double        LowerRAlimit;

double        UpperDecLimit;
double        LowerDecLimit;


double*             ra   = NULL;
double*            dec   = NULL;
double*           zobs   = NULL;
double*           zcos   = NULL;
double*            M_B   = NULL;
double*          zflag   = NULL;
int*              type   = NULL;
int*          photoMask  = NULL;

double*           zpec   = NULL;  // Value added catalogue parameters.
double*          zphot   = NULL;
double*          gal_z   = NULL;
double*         sampling = NULL;
double*    fkp_galweight = NULL;
double*   clip_galweight = NULL;

// derived parameters
bool*    Acceptanceflag = NULL;
double*          rDist  = NULL;
double*          xCoor  = NULL;
double*          yCoor  = NULL;
double*          zCoor  = NULL;


// -- Randoms --//
double                   alpha; // ratio of N_rand to N_gal pretty much. 

int          rand_number   = 0;
int          accepted_rand = 0;

double*      rand_ra     = NULL;
double*      rand_dec    = NULL;
double*      rand_chi    = NULL;
double*      rand_x      = NULL;  // Really need x,y,z?
double*      rand_y      = NULL;
double*      rand_z      = NULL;
double*      rand_weight = NULL;

// -- FKP weighting/normalisation --//
int    accepted_gals;
double       fkpPk;
double fkp_norm, daccepted_gals;


// Embedding volume for mock.
int          n0, n1, n2; // (z == 0), (x == 2).

double  AxisLimsArray[2][3];  // Array to hold the coordinate limits of the VIPERS survey.

double        xCellSize;            // Cell size, comoving distance, h^-1 Mpc
double        yCellSize;            // Cell size, comoving distance, h^-1 Mpc
double        zCellSize;            // Cell size, comoving distance, h^-1 Mpc 

double       CellVolume;

double        lo_zlim;              // Selection parameters. Volume limited sample between redshift 0.7 and 0.9
double        hi_zlim;
double          z_eff;
double           aexp;

int       boxlabel;
int          xlabel, ylabel, zlabel;


// -- n(z) calc. -- //
int          chibin_no;
double       chi_interval;
double       nz_smoothRadius; // smoothing length. 

double*      zbins    = NULL;
double*      chibins  = NULL;
double*      Nchi     = NULL;
double*      nbar     = NULL;
double*      nbar_2d  = NULL;
double*      comovVol = NULL;

double        cumulative_nbar[400];
double     cumulative_nbar_2d[400];
double    chi_cumulative_nbar[400];

// -- FFT units --//
double        kIntervalx; // fund_kx.
double        kIntervaly;
double        kIntervalz;

double        xNyquistWaveNumber; // Ny_kx.
double        yNyquistWaveNumber;
double        zNyquistWaveNumber;


// -- Binned modes -- //
int     kbin_no;

double  logk_min;
double  logk_max;
double  logk_interval;

int*    modes_perbin       = NULL;

double* mean_modk          = NULL;
double* binnedPk           = NULL;
double* logk_limits        = NULL;

// -- Multipole decomposition -- //
int     hiMultipoleOrder; // 0: use monopole only, 2: use quadrupole.

double* kLi; // L_2 evaluated for each individual mode. 
double* kM2; // Square NGP/CIC correction for each available mode.  
int*    mode_pkbin_index;  // index in which mode falls in binned p(k) array.

double* detA; //
double* Sum_Pi;
double* Sum_Li;
double* Sum_Li2;
double* Sum_PiLi;
int*    modes_perbin;

double* Hexadecapole;
double* Quadrupole;
double* Monopole;

// -- Clipping -- //
double  fraction_clipped;
double  appliedClippingThreshold;
double  clipping_smoothing_radius;

double* gal_clippingweights;


// -- Functions --
int          comovDistReshiftCalc();
double       SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval);

double       invert_StefanoBasis(double centreRA, double centreDec, double* xval, double* yval, double* zval);

int          prep_inverseCumulative_nbar();
double       inverse_cumulative_nbar(double arg);

int          CoordinateCalc();

int          JenkinsCoordinates();
int          JenkinsFold(double original[], int lenArray, int axis);
int          ApplyJenkins();


