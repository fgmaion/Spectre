#define ODD 0
#define EVEN 1
#define PARITY(a) (a)%2 ? ODD : EVEN
#define FFTLog_SWAP(a,b) {FFTLog_TMP = (a); (a) = (b); (b) = FFTLog_TMP;}

FFTLog_config*          eta_perp;
FFTLog_config*          eta_para;
FFTLog_config*        meanInfall;

double sigma2eta;



