#define ODD 0
#define EVEN 1
#define PARITY(a) (a)%2 ? ODD : EVEN
#define FFTLog_SWAP(a,b) {FFTLog_TMP = (a); (a) = (b); (b) = FFTLog_TMP;}

typedef struct FFTLog_complex{
  double re;
  double im;
  double amp;
  double arg;
}  FFTLog_complex;


typedef struct FFTLog_config{
  int N;
  
  double min;
  double max;
  double q;
  double mu;
  double kr;
  
  fftw_plan forwardplan;
  fftw_plan backwardplan;
  
  fftw_complex *pk;
  fftw_complex *xi;
  fftw_complex *cm;
  fftw_complex *um;
  fftw_complex *cmum;

  fftw_complex *input;
  fftw_complex *output;
  
  // two column array, first column will correspond to k values, 2nd to r. 
  double** krvals;

} FFTLog_config;

int     FFTlogRes;

double FFTLog_TMP;

void FFTLog(FFTLog_config *fc, fftw_plan p_forward, fftw_plan p_backward);
int  FFTLog_free(FFTLog_config *fc);

FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z);
FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu);

int FastLegendre_setInput(FFTLog_config *fc_mono, FFTLog_config *fc_quad, FFTLog_config *fc_hex, FFTLog_config *fc_oct, double beta, double velDispersion);

// FFTLog_config *fc;

FFTLog_config*  octupole_config;
FFTLog_config*      mono_config;
FFTLog_config*      quad_config;
FFTLog_config*      zero_config;
FFTLog_config*       hex_config;
FFTLog_config*       oct_config;
FFTLog_config*       dec_config;

FFTLog_config*   mono_config_ap;
FFTLog_config*   quad_config_ap;

FFTLog_config*   limmono_config;

FFTLog_config*  clipmono_config;
FFTLog_config*  clipquad_config;

FFTLog_config*    lnnorm_mono_config;
FFTLog_config*    lnnorm_quad_config;

FFTLog_config*    convlmonoCorr;
FFTLog_config*    convlquadCorr;
FFTLog_config*    convlhexCorr;

// Arrays of pre and post factors for xi to pk, pk to xi transforms. 
double* xi_mu_prefactor;
double* xi_mu_postfactor;
double* pk_mu_prefactor;
double* pk_mu_postfactor;

