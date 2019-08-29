int calc_delta(){
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;

  // cic assign the galaxies to overdensity, weighted by sampling, fkp weights and clipping weight.                          
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      cic_assign(1, xCoor[j], yCoor[j], zCoor[j], 1.);
    }
  }

  for(j=0; j<n0*n1*n2; j++){
    overdensity[j][0] /= (*pt2nz)(500.);
      
    overdensity[j][0] -=             1.;
  }

  // <1+\delta> = 1.
  double expected = 0.0;
  
  for(j=0; j<n0*n1*n2; j++)  expected += (1. + overdensity[j][0]); 
  
  expected /= 1.*n0*n1*n2;

  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = (1. + overdensity[j][0])/expected - 1.;
  
  fftw_execute(p);
  
  return 0;
}


double switchup(double kx, double ky, double kz, int pass_i, int pass_j){
  Interim = 1.;

  if(pass_i == 0)  Interim *= kz;
  if(pass_i == 1)  Interim *= ky;
  if(pass_i == 2)  Interim *= kx;

  if(pass_j == 0)  Interim *= kz;
  if(pass_j == 1)  Interim *= ky;
  if(pass_j == 2)  Interim *= kx;

  return Interim;
}


int tidal_tensor(){  
    int ndim = 3;
    
    double lambda_thres = 0.1;

    double T[3][3][256][256][256];

    fftw_complex* dummy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);

    gsl_vector*   eval  = gsl_vector_alloc(ndim);
    gsl_matrix*   evec  = gsl_matrix_alloc(ndim, ndim);

    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(ndim);

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index          = k*n1*n2 + j*n2 + i;
                
                H_k[Index][0] /= n0*n1*n2;
                H_k[Index][1] /= n0*n1*n2; 
            }
        }
    }

    for(kk=0; kk<3; kk++){
        for(jj=0; jj<3; jj++){
            for(k=0; k<n0; k++){
                for(j=0; j<n1; j++){
                    for(i=0; i<n2; i++){
                        k_x = kIntervalx*i;
                        k_y = kIntervaly*j;
                        k_z = kIntervalz*k;

                        if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
                        if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
                        if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

                        Index                              = k*n1*n2 + j*n2 + i;

                        kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);

                        kmodulus                           = pow(kSq, 0.5);

                        mu                                 = k_z/kmodulus;
                        if(kmodulus < 0.000001)       mu   = 0.0;

                        // Tij(k) = k_i x k_j \delta_k / k^2.
                        if(kSq > 0.)    H2_k[Index][0]     = H_k[Index][0]/kSq;
                        if(kSq > 0.)    H2_k[Index][1]     = H_k[Index][1]/kSq;

                        // Smooth density field with gaussian of width 4 h^-1Mpc prior to tidal tensor calc.
                        H2_k[Index][0]                    *= exp(-1.*kSq*0.5*pow(10., 2.));
                        H2_k[Index][1]                    *= exp(-1.*kSq*0.5*pow(10., 2.));

                        H2_k[Index][0]                    *= switchup(k_x, k_y, k_z, kk, jj);
                        H2_k[Index][1]                    *= switchup(k_x, k_y, k_z, kk, jj);
                   }
            }
        }

        int negkIndex;

        // Hermitian condition. One hemi-sphere is independent, e.g. k_z >= 0.                                                     
        for(k=n0-1; k>=n0/2; k--){
            for(j=0; j<n1; j++){
                for(i=0; i<n2; i++){
                    negkIndex          = k*n1*n2 + j*n2 + i;

                    Index              = 0;

                    // zero maps to zero on reflection through the origin.                                                         
                    if(i!=0)  Index   += (n2 - i);
                    if(j!=0)  Index   += (n1 - j)*n2;

                    Index   += (n0 - k)*n1*n2;

                    H2_k[negkIndex][0]  =     H2_k[Index][0];
                    H2_k[negkIndex][1]  = -1.*H2_k[Index][1];

                    if(negkIndex == Index)   H2_k[Index][1] = 0.0; // purely real                                                   
                }
          }
     }

    // in the k_z=0 plane one semi-circle is independent, k_y>0.                                                               
    for(j=n1-1; j>=n1/2; j--){
        for(i=0; i<n2; i++){
            negkIndex          = j*n2 + i;

            Index              = 0;

            // zero maps to zero on reflection through the origin.                                                             
            if(i!=0)  Index   += (n2 - i);
    
            Index   += (n1 - j)*n2;

            H2_k[negkIndex][0]  =     H2_k[Index][0];
            H2_k[negkIndex][1]  = -1.*H2_k[Index][1];

            if(negkIndex == Index)   H2_k[Index][1] = 0.0;
      }
   }

   // on the line k_z=k_y=0, one half is independent, k_x>=0.                                                                 
   for(i=n2-1; i>=n2/2; i--){
     negkIndex          = i;

     Index              = 0;

     // zero maps to zero on reflection through the origin.                                                                 
     Index             += (n2 - i);

     H2_k[negkIndex][0]  =      H2_k[Index][0];
     H2_k[negkIndex][1]  =  -1.*H2_k[Index][1];

     if(negkIndex == Index)     H2_k[Index][1] = 0.0;
   }

   iplan = fftw_plan_dft_3d(n0, n1, n2, H2_k, dummy, FFTW_BACKWARD, FFTW_ESTIMATE);

   fftw_execute(iplan);

   for(k=0; k<n0; k++){
     for(j=0; j<n1; j++){
       for(i=0; i<n2; i++){
         Index = k*n1*n2 + j*n2 + i;

         T[kk][jj][i][j][k] = dummy[Index][0];
       } 
    }
  }

  printf("\n%d \t %d", kk, jj);
 }
}

 double dimension[256][256][256];

 for(k=0; k<n0; k++){
   for(j=0; j<n1; j++){
     for(i=0; i<n2; i++){
       Index = k*n1*n2 + j*n2 + i;

       // all eigenvalues below threshold
       dimension[i][j][k] = 0.;

       double data[] = {  T[0][0][i][j][k],  T[0][1][i][j][k],  T[0][2][i][j][k],
			  T[1][0][i][j][k],  T[1][1][i][j][k],  T[1][2][i][j][k],
			  T[2][0][i][j][k],  T[2][1][i][j][k],  T[2][2][i][j][k],
                       };

       gsl_matrix_view m = gsl_matrix_view_array(data,ndim, ndim);

       gsl_eigen_symmv(&m.matrix, eval, evec, w);
       
       gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
       
       // printf("\n");
       
       for(ii=0; ii<3; ii++){
  	   double eval_i = gsl_vector_get(eval, ii);

	   // printf("%le \t", eval_i);

	   if(eval_i > lambda_thres) dimension[i][j][k] += 1.;
       }
     }
   }
 }


 sprintf(filepath, "%s/Data/500s/hod_cube/environ.dat", root_dir);

 output = fopen(filepath, "w"); 

 double projection[256][256];

 for(j=0; j<256; j++){
   for(i=0; i<256; i++){
     projection[i][j] = 0.;
     
     for(k=100; k<101; k++){
       Index = k*n1*n2 + j*n2 + i;
       
       projection[i][j] = dimension[i][j][k];
     }
 
     fprintf(output, "%le \t", projection[i][j]);  
   }

   fprintf(output, "\n");
 }
  
  fclose(output);



  // ** Density field **
  // Gaussian filter overdensity, result in smooth overdensity. 
  Gaussian_filter(4., 0);

  // <1+\delta> = 1.
  double expected = 0.0;
  
  for(j=0; j<n0*n1*n2; j++)  expected += (1. + smooth_overdensity[j][0]); 
  
  expected /= 1.*n0*n1*n2;

  for(j=0; j<n0*n1*n2; j++)  smooth_overdensity[j][0] = (1. + smooth_overdensity[j][0])/expected - 1.;  
  
  printf("\n<1+\delta> %lf", expected);
  
  
  sprintf(filepath, "%s/Data/500s/hod_cube/density.dat", root_dir);

  output = fopen(filepath, "w"); 

  for(j=0; j<256; j++){
    for(i=0; i<256; i++){
      projection[i][j] = 0.;
     
      for(k=100; k<101; k++){
        Index = k*n1*n2 + j*n2 + i;
       
        projection[i][j] = smooth_overdensity[Index][0];
      }
 
      fprintf(output, "%le \t", projection[i][j]);  
    }

    fprintf(output, "\n");
  }
  
  fclose(output);

  gsl_eigen_symmv_free(w);
 
  gsl_vector_free(eval);
 
  gsl_matrix_free(evec);

  double gal_invoid;

  /*
 // galaxies catalogue by environment.
 sprintf(filepath, "%s/Data/500s/hod_cube/galaxies_invoid.dat", root_dir);

 output = fopen(filepath, "w"); 

 for(j=0; j<Vipers_Num; j++){
    xlabel                  = (int) floor((xCoor[j] - AxisLimsArray[0][2])/xCellSize);    
    ylabel                  = (int) floor((yCoor[j] - AxisLimsArray[0][1])/yCellSize);
    zlabel                  = (int) floor((zCoor[j] - AxisLimsArray[0][0])/zCellSize);
 
    gal_invoid              = dimension[xlabel][ylabel][zlabel];
    
    if(gal_invoid == 0.)  fprintf(output, "%lf \t %lf \t %lf \n", xCoor[j], yCoor[j], zCoor[j]);
 }  
  
 fclose(output);
 */
 
 // randoms catalogue by environment.
 sprintf(filepath, "%s/Data/500s/hod_cube/randoms_insheet.dat", root_dir);

 output = fopen(filepath, "w"); 

 for(j=0; j<rand_number; j++){
    xlabel                  = (int) floor((rand_x[j] - AxisLimsArray[0][2])/xCellSize);    
    ylabel                  = (int) floor((rand_y[j] - AxisLimsArray[0][1])/yCellSize);
    zlabel                  = (int) floor((rand_z[j] - AxisLimsArray[0][0])/zCellSize);
 
    gal_invoid              = dimension[xlabel][ylabel][zlabel];
    
    if(gal_invoid == 1.)  fprintf(output, "%lf \t %lf \t %lf \n", rand_x[j], rand_y[j], rand_z[j]);
 }  
  
 fclose(output);
 
 return 0;
}
