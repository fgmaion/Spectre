int fft_qell(regress* inst){    
  int     pp, jj, kk, i;

  int      N=25;
  
  double  rs[N];
  double  Q0[N];
  double  Q2[N];

  for(j=0; j<N; j++){
    rs[j] = pow(10., -2. + j*5./N);
    
    Q0[j] = 0.0;
    Q2[j] = 0.0;
  }

  walltime("Start of Q2_ell calc.");

  for(j=0; j<n0*n1*n2; j++)     overdensity[j] = 0.0;

  double      x2, z2;
  double c_ra, c_dec;

  c_ra    =  CentreRA*(pi/180.);
  c_dec   = CentreDec*(pi/180.);
  
  for(int kk=0; kk<5; kk++){
    #pragma omp parallel for private(j, F, x2, z2) if(thread == 1)
    for(j=0; j<rand_number; j++){
      drand48_r(&randBuffers[omp_get_thread_num()], &F);

      rand_chi[j]    = inverse_cumulative_nbar(F);

      x2             =  rand_chi[j]*cos(rand_ra[j])*cos(rand_dec[j]);
      z2             = -rand_chi[j]*sin(rand_dec[j]);
    
      rand_x[j]      = -sin(c_dec)*x2  - cos(c_dec)*z2;
      rand_z[j]      =  cos(c_dec)*x2  - sin(c_dec)*z2;

      rand_y[j]      =  rand_chi[j]*sin(rand_ra[j])*cos(rand_dec[j]);

      rand_x[j]     +=  stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
      rand_y[j]     +=  stefano_trans_y;
      rand_z[j]     +=  stefano_trans_z;

      rand_weight[j] = 1./(1. + (*pt2nz)(rand_chi[j])*fkpPk);
    
      ngp_assign(rand_x[j], rand_y[j], rand_z[j], rand_weight[j]);
    }
  }
    
  fftw_execute(plan);

  // mu undefined at the origin. 
  inst->kLi[0] = 0.0;
  
  #pragma omp parallel for reduction(+: Q0[:N], Q2[:N]) private(Index, k, j, i, kk, jj, pp, kSq, kmodulus) if(thread == 1)
  for(k=0; k<n0; k++){
    for(j=0; j<n1; j++){
      kk = k;
      jj = j;
        
      if(k > n0/2)  kk -= n0;
      if(j > n1/2)  jj -= n1;
        
      for(i=0; i<nx; i++){
        Index                        = k*n1*nx + j*nx + i;
        
        kSq                          = pow(kk, 2.) + pow(jj, 2.) + pow(i, 2.);
        
        kmodulus                     = fund_kz*pow(kSq, 0.5); // units of fundamental interval.

        H_k[Index][0]               /= inst->kM2[Index];  // Correct mass assignment of randoms; cic = 2, ngp = 1.
        H_k[Index][1]               /= inst->kM2[Index];
          
        H_k[Index][0]                = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.); // where is n0*n1*n2 factor?       
        
        for(pp=0; pp<N; pp++){
          Q0[pp]                    += H_k[Index][0]*gsl_sf_bessel_jl(0, kmodulus*rs[pp]); // L0 is one.
          Q2[pp]                    += H_k[Index][0]*gsl_sf_bessel_jl(2, kmodulus*rs[pp])*inst->kLi[Index];
        }
      }
    }
  }

  for(j=0; j<N; j++)  Q2[j] *= -5.0; // i^q (2q +1).  

  sprintf(filepath,"%s/W1_Spectro_V7_4/Qmultipoles/FFT/W%d_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf.dat", root_dir, fieldFlag, lo_zlim, hi_zlim); 

  output = fopen(filepath, "w");
  
  for(j=0; j<N; j++)  fprintf(output, "%.6le \t %.6le \t %.6le \n", rs[j], Q0[j]/Q0[0], Q2[j]/Q0[0]);

  fclose(output);

  return 0;
}

  

