int fft_qellk(regress* inst){    
  double      x2, z2;
  double       c_dec;
  
  walltime("Start of Q2_ell(k) calc.");

  for(j=0; j<n0*n1*n2; j++)  overdensity[j] = 0.0;

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
    }

    for(j=0; j<rand_number; j++)  ngp_assign(rand_x[j], rand_y[j], rand_z[j], rand_weight[j]);
  }
    
  fftw_execute(plan);

  nosort_MultipoleCalc(inst);

  print_multipoles(inst);
  
  return 0;
}

  

