double fracArea; // W%d/Total area fraction. 

// Misc //
int    data_mock_flag;

int order, nrotations, mono_order;

double* xdata;
double****** xtheory;

double* ydata;
double****** ytheory;

double*    kdata = NULL;  // (kVals, kVals) vector for eff. k calc. after covariance diagonalisation. 
double*   dkdata = NULL;

double**     Multipoles;
double**    dMultipoles;

//-- binned pk. --//
int     kbin_no;
double* kVals;

int*   fftlog_indices;

// -- Jenkins -- //
int     jenkins_foldIndex_foldedfile;
int     jenkins_foldIndex_unfoldedfile;
double  jenkins_fold_kjoin;

//-- Likelihood grid --//
int        Res;
double    dRes;

int     Res_ap;
double dRes_ap;

double    paramNumber; // why double?

double               A11Sq;
double             fsigma8; // fsig8.
double             bsigma8;
double           alpha_pad;
double          linearBias;
double         epsilon_pad;
double       velDispersion;

double      fsigma8Interval; // fsig8_interval.
double      bsigma8Interval;
double        sigmaInterval;
double    alpha_padInterval;
double  epsilon_padInterval;
double        A11SqInterval;

double  min_fsigma8;
double  max_fsigma8;

double  min_bsigma8;
double  max_bsigma8;

double    min_velDisperse;
double    max_velDisperse;

double    min_alpha_pad;  // Inclusion of Alcock-Paczynski in Likelihood analysis.
double    max_alpha_pad;

double  min_epsilon_pad;
double  max_epsilon_pad;

double      min_A11Sq;
double      max_A11Sq;

double       minChiSq;
double ChiSq_expected; 

double***** ChiSqGrid;

double    minChiSq_fsigma8; // min_x2_fsig8.
double    minChiSq_A11Sq;
double    minChiSq_sigma;
double    minChiSq_bsigma8;
double    minChiSq_alpha_pad;
double    minChiSq_epsilon_pad;

//-- Model eval. --//
// Spline and squential splint of D^2(k*sigma*mu) Kaiser-Lorentz dispersion factors. 
int     sp_kL_N;
double* sp_kL_ks;

double* sp_kL_mu0;
double* sp_kL_mu2;
double* sp_kL_mu4;
double* sp_kL_mu6;
double* sp_kL_mu8;
double* sp_kL_mu10;

double* sp_kL_mu0_2d;
double* sp_kL_mu2_2d;
double* sp_kL_mu4_2d;
double* sp_kL_mu6_2d;
double* sp_kL_mu8_2d;
double* sp_kL_mu10_2d;

// Functions pre-computed for FFTlog k or r vals. 
double* FFTlog_Pk;

double* FFTlog_W0; // W_0(r) evaluated at FFTlog rvals.
double* FFTlog_W2;
double* FFTlog_W4;
double* FFTlog_W6;

double* FFTlog_W0_joint; // change to FFTlog_jW0
double* FFTlog_W2_joint;
double* FFTlog_W4_joint;
double* FFTlog_W6_joint;

double* FFTlog_Wk0; // tilde W_0(k) etc. 
double* FFTlog_Wk2;


// -- AP correction -- //
double mono_epsilonCorrection_1, mono_epsilonCorrection_2;  // Alcock-Paczynski, epsilon dependent correction terms.
double quad_epsilonCorrection_1, quad_epsilonCorrection_2;

double* spline_lnk; // spline of dln P_R by d lnk. 
double* dlnPR_dlnk;
double* dlnPR_dlnk_2D;


// -- Clipping correction -- //
int           varcalc_index;
int cnvldpk_zeropoint_index;

double         clipmono_amp;
double       clip_distcoeff;

// -- Likelihood -- //
double  app_sigma8;

double  ChiSq_kmax;
double  ChiSq_kmin;

gsl_matrix* Covariance;
gsl_matrix* sigma_norm;

gsl_vector* eval;
gsl_matrix* evec;

gsl_eigen_symmv_workspace* w;

gsl_vector* col;

double*   MeanMultipoles;

double*  eigenVals;
double** eigenVecs;

double  smallestEigenvalue;

int     chiSq_kmaxIndex;
int     chiSq_kminIndex;


// -- Posteriors -- //
double  PosteriorNorm = 0.0; // variable to normalise posteriros. 

double* bsigma8Posterior;
double* fsigma8Posterior;
double*   sigmaPosterior; // velocity dispersion posterior. 



// -- Functions -- //
double splintLinearPk(double k);
double splintMatterPk(double EvalPoint);

double (*pt2Pk)(double k) = NULL;

double kaiserGauss_Monofactor(double ks, double beta);
double kaiserGauss_Quadfactor(double ks, double beta);
double kaiserGauss_Hexfactor(double ks, double beta);

double AP_P0(double kprime, double beta, double sigma, double epsilon, double local_alpha);
double AP_P2(double kprime, double beta, double sigma, double epsilon, double local_alpha);
double AP_P4(double kprime, double beta, double sigma, double epsilon, double local_alpha);

double ChiSqEval();
double ChiSqEval_ap();


// -- Posteriors -- //
double calc_fsigma8Posterior();
double calc_bsigma8Posterior();
double calc_velDispPosterior();
