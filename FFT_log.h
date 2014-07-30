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
  
  fftw_plan pk_forwardplan;
  fftw_plan pk_backwardplan;
  
  fftw_plan xi_forwardplan;
  fftw_plan xi_backwardplan;
  
  fftw_complex *pk;
  fftw_complex *xi;
  fftw_complex *cm;
  fftw_complex *um;
  fftw_complex *cmum;
} FFTLog_config;

double FFTLog_TMP;

void FFTLog(FFTLog_config *fc, fftw_plan p_forward, fftw_plan p_backward);
void FFTLog_free(FFTLog_config *fc);

FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z);
FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu);

// FFTLog_config *fc;

FFTLog_config*  octupole_config;
FFTLog_config*      mono_config;
FFTLog_config*      quad_config;
FFTLog_config*       hex_config;

FFTLog_config*  clipmono_config;
FFTLog_config*  clipquad_config;
