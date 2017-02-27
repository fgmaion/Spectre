char   logfilepath[200]; // filepath to log file for pair counting.
char    surveyType[200];

int VIPERS_kSpace_multipoles_lineNo;

double* VIPERS_k; // spline and splint \tilde W2(k) 
double* VIPERS_kMono;
double* VIPERS_kMono2D;
double* VIPERS_kQuad;
double* VIPERS_kQuad2D;

int     VIPERS_mask_lineNo_lo; // spline and spline \tilde W2(k).
int     VIPERS_jmask_lineNo_lo;

double* VIPERS_maskr_lo;
double* VIPERS_jmaskr_lo;

double* VIPERS_maskMono_lo;
double* VIPERS_jmaskMono_lo;

double* VIPERS_maskQuad_lo;
double* VIPERS_jmaskQuad_lo;

double* VIPERS_maskHex_lo;
double* VIPERS_jmaskHex_lo;

double* VIPERS_maskOct_lo;
double* VIPERS_jmaskOct_lo;

double* VIPERS_maskDec_lo;
double* VIPERS_jmaskDec_lo;

double* VIPERS_maskMono2D_lo;
double* VIPERS_jmaskMono2D_lo;

double* VIPERS_maskQuad2D_lo;
double* VIPERS_jmaskQuad2D_lo;

double* VIPERS_maskHex2D_lo;
double* VIPERS_jmaskHex2D_lo;

double* VIPERS_maskOct2D_lo;
double* VIPERS_jmaskOct2D_lo;

double* VIPERS_maskDec2D_lo;
double* VIPERS_jmaskDec2D_lo;

double mask_monopolenorm_lo;
double jmask_monopolenorm_lo;

int     VIPERS_mask_lineNo_hi;
int     VIPERS_jmask_lineNo_hi;

double* VIPERS_maskr_hi;
double* VIPERS_maskMono_hi;
double* VIPERS_maskQuad_hi;
double* VIPERS_maskHex_hi;
double* VIPERS_maskOct_hi;
double* VIPERS_maskDec_hi;

double* VIPERS_jmaskr_hi;
double* VIPERS_jmaskMono_hi;
double* VIPERS_jmaskQuad_hi;
double* VIPERS_jmaskHex_hi;
double* VIPERS_jmaskOct_hi;
double* VIPERS_jmaskDec_hi;

double* VIPERS_maskMono2D_hi;
double* VIPERS_maskQuad2D_hi;
double* VIPERS_maskHex2D_hi;
double* VIPERS_maskOct2D_hi;
double* VIPERS_maskDec2D_hi;

double* VIPERS_jmaskMono2D_hi;
double* VIPERS_jmaskQuad2D_hi;
double* VIPERS_jmaskHex2D_hi;
double* VIPERS_jmaskOct2D_hi;
double* VIPERS_jmaskDec2D_hi;

double mask_monopolenorm_hi;
double jmask_monopolenorm_hi;

double  loRes_highRes_join;
double jloRes_highRes_join;

int     VIPERS_mask_lineNo_hihi;
double* VIPERS_maskr_hihi;
double* VIPERS_maskMono_hihi;
double* VIPERS_maskQuad_hihi;
double* VIPERS_maskHex_hihi;
double* VIPERS_maskOct_hihi;
double* VIPERS_maskDec_hihi;

int     VIPERS_jmask_lineNo_hihi;
double* VIPERS_jmaskr_hihi;
double* VIPERS_jmaskMono_hihi;
double* VIPERS_jmaskQuad_hihi;
double* VIPERS_jmaskHex_hihi;
double* VIPERS_jmaskOct_hihi;
double* VIPERS_jmaskDec_hihi;

double* VIPERS_maskMono2D_hihi;
double* VIPERS_maskQuad2D_hihi;
double* VIPERS_maskHex2D_hihi;
double* VIPERS_maskOct2D_hihi;
double* VIPERS_maskDec2D_hihi;

double* VIPERS_jmaskMono2D_hihi;
double* VIPERS_jmaskQuad2D_hihi;
double* VIPERS_jmaskHex2D_hihi;
double* VIPERS_jmaskOct2D_hihi;
double* VIPERS_jmaskDec2D_hihi;

double mask_monopolenorm_hihi;

double loRes_highRes_join;
double hiRes_hihiRes_join;

double jmask_monopolenorm_hihi;

double jloRes_highRes_join;
double jhiRes_hihiRes_join;


//-- Functions --//
double splint_unit_maskMultipoles(double r, int transformOrder);

double splint_VIPERS_maskMono(double r);
double splint_VIPERS_maskQuad(double r);
double splint_VIPERS_maskHex(double r);
double splint_VIPERS_maskOct(double r);
double splint_VIPERS_maskDec(double r);

double splint_VIPERS_maskMultipoles(double r, int transformOrder);

double splint_VIPERS_jmaskMono(double r);
double splint_VIPERS_jmaskQuad(double r);
double splint_VIPERS_jmaskHex(double r);
double splint_VIPERS_jmaskOct(double r);
double splint_VIPERS_jmaskDec(double r);

double splint_VIPERS_jmaskMultipoles(double r, int transformOrder);

double splint_VIPERS_kSpaceMono(double k);
double splint_VIPERS_kSpaceQuad(double k);

double (*pt2maskMultipoles)(double r, int transformOrder) = NULL;
