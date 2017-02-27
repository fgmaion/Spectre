// #include <omp.h>
#include <stdio.h>
#include <stdlib.h>
//#include <fftw3.h>
// #include <time.h>
#include "/home/mjw/Aux_functions/header.h"


int c2c(int N, double array[]){
  /*                                                                                                                                                                                         
    One-Dimensional DFTs of Real Data                                                                                                                                                        
    -- “Hermitian” redundancy: out[i] is the conjugate of out[n-i].                                                                                                                         
    --  input and output arrays are of different sizes and types                                                                                                                            
    --  input is n real numbers, while the output is n/2 + 1 complex numbers; n/2 modes + 0                                                                                                 
    --  requires slight “padding” of the input array for in-place transforms                                                                                                                 
    --  inverse transform (complex to real) has the side-effect of overwriting its input array, by default                                                                                   
  */
  
  fftw_complex*    in;
  fftw_complex*   out;

  fftw_plan      plan;

  in  = fftw_malloc(N*N*N*sizeof(fftw_complex));
  out = fftw_malloc(N*N*N*sizeof(fftw_complex));
  
  // fftw_init_threads();
  // fftw_plan_with_nthreads(omp_get_max_threads());

  // fftw_import_wisdom_from_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  plan  = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);  // no sign argument. Instead, r2c DFTs are always FFTW_FORWARD and c2r DFTs are always FFTW_BACKWARD

  // fftw_export_wisdom_to_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  for(j=0; j<N*N*N; j++){
    in[j][0] = array[j];
    in[j][1] = 0.0;
  }
  
  fftw_execute(plan);

  for(j=0; j<N*N*N; j++)  printf("\n%+.4lf \t %+.4lf", out[j][0], out[j][1]);

  return 0;
}


int r2c(int N, double array[]){
  /*                                                                                                                                                                                         
    One-Dimensional DFTs of Real Data                                                                                                                                                            -- “Hermitian” redundancy: out[i] is the conjugate of out[n-i].                                                                                                                              --  input and output arrays are of different sizes and types                                                                                                                                 --  input is n real numbers, while the output is n/2 + 1 complex numbers; n/2 modes + 0                                                                                                      --  requires slight “padding” of the input array for in-place transforms                                                                                                                     --  inverse transform (complex to real) has the side-effect of overwriting its input array, by default                                                                                   
  */

  /*
    Multi-Dimensional DFTs of Real Data.

    http://www.fftw.org/fftw3_doc/Multi_002dDimensional-DFTs-of-Real-Data.html#Multi_002dDimensional-DFTs-of-Real-Data
    
    -- complex output array is cut roughly in half and the real array requires padding for in-place transforms (as in 1d, above).


  */
  
  double*          in;
  fftw_complex*   out;

  fftw_plan      plan;

  in  = fftw_malloc(N*N*N*sizeof(double));
  out = fftw_malloc(N*N*(N/2 + 1)*sizeof(fftw_complex));

  // fftw_init_threads();
  // fftw_plan_with_nthreads(omp_get_max_threads());      

  // fftw_import_wisdom_from_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  plan  = fftw_plan_dft_r2c_3d(N, N, N, in, out,  FFTW_ESTIMATE);  // no sign argument. Instead, r2c DFTs are always FFTW_FORWARD and c2r DFTs are always FFTW_BACKWARD

  // fftw_export_wisdom_to_filename("/home/mjw/HOD_MockRun/wisdom/wisdom.dat");

  for(j=0; j<N*N*N; j++)  in[j] = array[j];
  
  fftw_execute(plan);

  for(j=0; j<N*N*(N/2+1); j++)  printf("\n%+.4lf \t %+.4lf", out[j][0], out[j][1]);
  
  return 0;
}


int main(int argc, char *argv[]){
  double begin  = getRealTime();

  int N = 4;
  
  double data[N*N*N];

  for(j=0; j<N*N*N; j++)  data[j] = (double) rand() / (double)RAND_MAX;
  
  c2c(N, data);

  printf("\n\n");

  r2c(N, data);
  
  double   end = getRealTime();

  printf("\n\nWall time: %.6lf", end - begin);
  
  printf("\n\n");
  
  return 0;
}
