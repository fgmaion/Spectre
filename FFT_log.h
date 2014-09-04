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

double FFTLog_TMP;

void FFTLog(FFTLog_config *fc, fftw_plan p_forward, fftw_plan p_backward);
int FFTLog_free(FFTLog_config *fc);

FFTLog_complex FFTLog_U_mu(double mu, FFTLog_complex z);
FFTLog_config *FFTLog_init(int N, double min, double max, double q, double mu);

// FFTLog_config *fc;

FFTLog_config*  octupole_config;
FFTLog_config*      mono_config;
FFTLog_config*      quad_config;
FFTLog_config*       hex_config;

FFTLog_config*   limmono_config;

FFTLog_config*  clipmono_config;
FFTLog_config*  clipquad_config;

FFTLog_config*    convlmonoCorr;
FFTLog_config*    convlquadCorr;
