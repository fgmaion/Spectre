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

// -- Randoms --//
double*      rand_ra     = NULL;
double*      rand_dec    = NULL;
double*      rand_chi    = NULL;
double*      rand_x      = NULL;  // Really need x,y,z?
double*      rand_y      = NULL;
double*      rand_z      = NULL;
double*      rand_weight = NULL;

// -- FKP weighting/normalisation --//
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

// -- n(z) calc. -- //
double nz_smoothRadius; // smoothing length. 

double        cumulative_nbar[400];
double     cumulative_nbar_2d[400];
double    chi_cumulative_nbar[400];



// -- //
double logk_min;
double logk_max;
double logk_interval;


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

int          JenkinsCoordinates();
int          JenkinsFold(double original[], int lenArray, int axis);
int          ApplyJenkins();

