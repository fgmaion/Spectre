// changes: 
// 1) alpha includes fkp_weights. 
// 2) // alpha from including clipping weights to not. 

int ncalc_clippingweights(){
  double                       chi;
  double*             cell_weights;
  double*            rand_occupied;

  double     number_occupied = 0.0;
  
  cell_weights  = malloc(n0*n1*n2*sizeof(*cell_weights));
  rand_occupied = malloc(n0*n1*n2*sizeof(*rand_occupied));

  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;

  for(j=0; j<n0*n1*n2; j++)   rand_occupied[j] = 0.0;

  // initialise to no clipping.                                                                                                                                                                                                            
  for(j=0; j<n0*n1*n2; j++)    cell_weights[j] = 1.0;

  // binary array, is surveyed or not.                                                                                                                                                                                                     
  for(j=0; j<rand_number; j++){
    if(rand_accept[j] == true){
      // bound survey volume.                                                                                                                                                                                                              
      boxlabel                              =                 boxCoordinates(rand_x, rand_y, rand_z, j);

      //cell_chi[boxlabel]                   +=                                               rand_chi[j];                                                                                                                                 

      rand_occupied[boxlabel]               =                                                         1;
    }
  }

  // Number of cells over which we average <1+ delta> rescaling.                                                                                                                                                                           
  for(j=0; j<n0*n1*n2; j++){
    if(rand_occupied[j] > 0.){
      number_occupied += 1.;
    }
  }

  printf("\n\nvolume available: (%.2lf h^-1 Mpc)^3.", pow(CellVolume*number_occupied, 1./3.));

  // assign galaxies using NGP.  Overdensity contains the galaxy counts.                                                                                                                                                                   
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      chi                               =                          interp_comovingDistance(zobs[j]);

      boxlabel                          =                    boxCoordinates(xCoor, yCoor, zCoor, j);

      // Add (1./sampling) in units of nbar.                                                                                                                                                                                           
      overdensity[boxlabel][0]         +=                       pow((*pt2nz)(chi)*sampling[j], -1.);
    }
  }

  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]      /=            CellVolume;

  // Smooth n/<n>.  True/false flag for enforcing zero mean for the filtered field.                                                                                                
  Gaussian_filter(clipping_smoothing_radius, 0);
  
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0]      -=                   1.0;  // -1 or 0 determined by randoms.                                                                                                                           

  // Scale (1 + delta) such that <1+ delta> = 1.; i.e. homogeneous "zero mean constraint".                                                                                                                                                 
  double norm = 0.0;
  
  for(j=0; j<n0*n1*n2; j++){
    // even if randoms don't fully sample the field then still a volume average over a reasonably representative sample.                                                                                                                   
    if(rand_occupied[j] > 0.){
      norm += smooth_overdensity[j][0]; // currently "smooth_overdensity" is really smoothed (1 + delta)                                                                                                                                   
    }
  }

  norm /= number_occupied;

  printf("\n\nMean renormalisation in clipping weights: %.2lf", norm);

  // rescaling of <1+delta>.  by rescaling (1+d), don't slip below d=-1.                                                                                                                                                                   
  for(j=0; j<n0*n1*n2; j++)  smooth_overdensity[j][0]  /=                                         norm;

  // now it's delta.                                                                                                                                                                                                                       
  for(j=0; j<n0*n1*n2; j++)  smooth_overdensity[j][0]  -=                                          1.0;

  for(j=0; j<n0*n1*n2; j++){
    if(smooth_overdensity[j][0] > appliedClippingThreshold){
      cell_weights[j]                   = (1. + appliedClippingThreshold)/(1. + overdensity[j][0]);

      if(rand_occupied[j] > 0.)  fraction_clipped  += 1.0;
    }
  }

  fraction_clipped /= number_occupied;

  printf("\n\napplied clipping threshold: %lf, fraction of cells clipped: %lf", appliedClippingThreshold, fraction_clipped);
  
  double mass_frac = 0.0, frac_norm = 0.0;

  if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_3/mocks_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
  if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_3/data_v1.7/clip_weights/W%d/clip_wghts_d0_%.2lf_z_%.1lf_%.1lf_%d_256_pk.dat", root_dir, fieldFlag, appliedClippingThreshold, lo_zlim, hi_zlim, loopCount);
  
  int count = 0;

  output = fopen(filepath, "w");

  // Clipping weights should not be renormalised, they simply multiply the (fkp_weight/sampling) without further renormalisation.                                                                                                       
  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j]  == true){
	  boxlabel           = boxCoordinates(xCoor, yCoor, zCoor, j);

	  clip_galweight[j]  =                 cell_weights[boxlabel];

	  if(clip_galweight[j] > 1.)  printf("\nspurious weight: %.2lf", clip_galweight[j]);

	  frac_norm         +=                                    1.0;

	  mass_frac         +=                 cell_weights[boxlabel];

	  // if(rand_occupied[boxlabel] == 0.){  printf("%d \t %.3lf \n", count, clip_galweight[j]); count +=1;}                                                                                                                           
    }

    else{clip_galweight[j] = 0.0;}

    fprintf(output, "%le \n", clip_galweight[j]);
  }

  fclose(output);

  printf("\n\nMean weight: %.3lf", mass_frac/frac_norm);
      
  return 0;
}


int nalpha_calc(double daccepted_gals){
  accepted = 0;

  daccepted_gals     = 0.0;

  double fkp_daccepted_gals = 0.0;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      accepted         += 1;

      
      // daccepted_gals   += 1./sampling[j];                                                                                                                                                                                               

      // as of May 13th 2016.                                                                                                                                                                                                              
      daccepted_gals   += clip_galweight[j]/sampling[j];

      fkp_daccepted_gals += (clip_galweight[j]/sampling[j])/(1. + (*pt2nz)(rDist[j])*fkpPk); 
    }
  }

  printf("\n\nTotal number of galaxies on input: %d, accepted: %d, accepted (weighted) %.2lf", Vipers_Num, accepted, daccepted_gals);

  // alpha = 1.*daccepted_gals/accepted_rand;

  alpha = fkp_daccepted_gals/fkp_accepted_rand;

  printf("\n\nalpha %.4f", alpha);

  printf("\n\naccepted. inverted, rotated & translated");

  return 0;
}


int ncalc_overdensity(){
  // Clean overdensity before galaxy assignment.                                                                                                                                                                                           
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
  for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], (1./sampling[j])*fkp_galweight[j]*clip_galweight[j]);
    }
  }

  for(j=0; j<n0*n1*n2; j++)  surveyMask[j] = 0.0;

  for(j=0; j<rand_number; j++){
    // assign randoms to surveyMask weighted by fkp.                                                                                                                                                                                     
    if(rand_accept[j] == true)     cic_assign(0, rand_x[j], rand_y[j], rand_z[j], rand_weight[j]);
  }

  return 0;
}

/*
int ncalc_overdensity2(){
  // Clean overdensity before galaxy assignment. 
  double* mean_z, int* n_gal;

  mean_z = malloc(n0*n1*n2*sizeof(double));
  n_gal  = malloc(n0*n1*n2*sizeof(double));
                                                                                                                                                               
  for(j=0; j<n0*n1*n2; j++){   
    overdensity[j][0] = 0.0;
    overdensity[j][1] = 0.0;

    mean_z[j]         = 0.0;

    n_gal[j]          =   0;
  }

  for(j=0; j<Vipers_Num; j++){
    if(Acceptanceflag[j] == true){
      cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], (1./sampling[j])*clip_galweight[j]/CellVolume);


      boxlabel = boxCoordinates(xCoor, yCoor, zCoor, j);

      mean_z[boxlabel] += zobs[j];

      n_gal[boxlabel]  += 1;
    }
  }

  for(j=0; j<n0*n1*n2; j++){
    if(n_gal[j] > 0)  mean_z[j] /= n_gal[j];
  }

  for(j=0; j<n0*n1*n2; j++)  surveyMask[j] = 0.0;

  for(j=0; j<rand_number; j++){
   // assign randoms to surveyMask weighted by fkp.                                                                                                                                                                                       
    if(rand_accept[j] == true){    
      boxlabel                          =                    boxCoordinates(rand_x, rand_y, rand_z, j);

      surveyMask[boxlabel]              = 1;
    }
  }

  for(j=0; j<n0*n1*n2; j++){
    if(surveyMask[j] > 0.){
      chi                               =                          interp_comovingDistance(zobs[j]);

      overdensity[j][0] /=        (*pt2nz)(interp_comovingDistance(zobs[j]))

    }
  }

  return 0;
}
*/

int nPkCalc(){
   fftw_execute(p);

   // H2_k assigned to hold fft of gals.                                                                                                                                                                                                    
   for(j=0; j<n0*n1*n2; j++){
     H2_k[j][0] = H_k[j][0];
     H2_k[j][1] = H_k[j][1];
   }

   for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = surveyMask[j];

   fftw_execute(p);

   // free overdensity.                                                                                                                                                                                                                     
   // free_grid();                                                                                                                                                                                                                          

   assign2DPkMemory();

   nPkCorrections();

   observedQuadrupole(polar_pkcount);

   return 0;
 }


int nPkCorrections(){
   double pk, WindowFunc;

   double rand_shot = 0.0, gal_shot = 0.0;

   // shot noise from randoms cat. calculation.                                                                                                                                                                                             
   for(j=0; j<rand_number; j++){
     if(rand_accept[j]    == true)  rand_shot += pow(rand_weight[j], 2.);
   }

   // oversight:  could have determined shot noise from randoms aswell.  will scale with alpha.
   rand_shot *= alpha*alpha;

   // shot noise from galaxies cat. calculation, including angular sampling.                                                                                                                                                                
   for(j=0; j<Vipers_Num; j++){
     if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j]/sampling[j], 2.);
   }

   printf("\n\nShot noise contributions: Randoms %.4lf, Galaxies %.4lf", rand_shot, gal_shot);

   //  NOTE:   Clipping is assumed to not alter the galaxy shot noise                                                                                                                                                                       
   // as (a)   The volume affected is very small (~ 1 %)                                                                                                                                                                                    
   //    (b)   Galaxies removed account for linear fluctuations to be                                                                                                                                                                       
   //          smaller, as opposed to a change in number density.                                                                                                                                                                           

   double newZero = H_k[0][0];

   polar_pkcount = 0;

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

	 WindowFunc                         = 1.;

	 // 0/0 - > Taylor expand.                                                                                                                                                                                                    
	 if(k_x != 0.)  WindowFunc         *= sin(pi*k_x*0.5/xNyquistWaveNumber)/(pi*k_x*0.5/xNyquistWaveNumber);

	 if(k_y != 0.)  WindowFunc         *= sin(pi*k_y*0.5/yNyquistWaveNumber)/(pi*k_y*0.5/yNyquistWaveNumber);

	 if(k_z != 0.)  WindowFunc         *= sin(pi*k_z*0.5/zNyquistWaveNumber)/(pi*k_z*0.5/zNyquistWaveNumber);

	 // Correct for mass assignment of randoms.                                                                                                                                                                                   
	 // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2.                                                                                                                                                    
	 H_k[Index][0]                     /= pow(WindowFunc, 2.);
	 H_k[Index][1]                     /= pow(WindowFunc, 2.);

	 // If smoothing the field, do not correct for mass assignment.  Obvious in hindsight.                                                                                                                                        
	 // Cloud in cell. NGP corresponds to WindowFunc rather than WindowFunc^2.                                                                                                                                                    
	 H2_k[Index][0]                    /= pow(WindowFunc, 2.);
	 H2_k[Index][1]                    /= pow(WindowFunc, 2.);

	 // gals          // normalised rands.                                                                                                                                                     
	 H_k[Index][0]                      = H2_k[Index][0] - alpha*H_k[Index][0];
	 H_k[Index][1]                      = H2_k[Index][1] - alpha*H_k[Index][1];

	 // enforce k=0 mode has P(0) = 0.                                                                                                                                                                                            
	 // H_k[Index][0]                      = H2_k[Index][0] - H_k[Index][0]*(H2_k[0][0]/newZero);                                                                                                                                 
	 // H_k[Index][1]                      = H2_k[Index][1] - H_k[Index][1]*(H2_k[0][0]/newZero);                                                                                                                                 

	 pk                                 = pow(H_k[Index][0], 2.) + pow(H_k[Index][1], 2.);

	 pk                                -= rand_shot;

	 if(kmodulus > 0.000001){
	   //// Only half the modes are independent. ////                                                                                                                                                                       
	   if(k_z>0.){
	     // One hemi-sphere is independent, e.g. k_z >= 0.                                                                                                                                                                
	     polar_pk[polar_pkcount][0]   = kmodulus;
	     polar_pk[polar_pkcount][1]   = fabs(mu);
	     polar_pk[polar_pkcount][2]   = pk;

	     twodim_pk[polar_pkcount][0]  = fabs(k_z);
	     twodim_pk[polar_pkcount][1]  = pow(k_y*k_y + k_x*k_x, 0.5);
	     twodim_pk[polar_pkcount][2]  = pk;

	     polar_pkcount               += 1;
	   }

	   else if((k_z == 0.0) && (k_y > 0.0)){
	     // in the k_z=0 plane one semi-circle is independent, k_y>0.                                                                                                                                                         
	     polar_pk[polar_pkcount][0]   = kmodulus;
	     polar_pk[polar_pkcount][1]   = fabs(mu);
	     polar_pk[polar_pkcount][2]   = pk;

	     twodim_pk[polar_pkcount][0]  = fabs(k_z);
	     twodim_pk[polar_pkcount][1]  = pow(k_y*k_y + k_x*k_x, 0.5);
	     twodim_pk[polar_pkcount][2]  = pk;

	     polar_pkcount                += 1;
	   }

	   else if((k_z == 0.0) && (k_y == 0.0) && (k_x > 0.0)){
	     // on the line k_z=k_y=0, one half is independent, k_x>=0.                                                                                                                                                   
	     // in the k_z=0 plane one semi-circle is independent, k_y>0.                                                                                                                         
	     polar_pk[polar_pkcount][0]    = kmodulus;
	     polar_pk[polar_pkcount][1]    = fabs(mu);
	     polar_pk[polar_pkcount][2]    = pk;

	     twodim_pk[polar_pkcount][0]   = fabs(k_z);
	     twodim_pk[polar_pkcount][1]   = pow(k_y*k_y + k_x*k_x, 0.5);
	     twodim_pk[polar_pkcount][2]   = pk;

	     polar_pkcount                += 1;
	   }

	   // if((0.01923825<kmodulus) && (kmodulus<0.02407176))  printf("\n%e \t %e \t %e", kmodulus, mu, pk);                                                                                                             
	 }
       }
     }
   }

   // free_HOD();                                                                                                                                                                                                                           

   return 0;
}
