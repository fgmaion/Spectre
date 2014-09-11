double fiberCollision(double rx, double ry, double dx, double dy){
  if((fabs(rx)>dx) && (fabs(ry)>dy)){
    return 1.;
  }

  else{
    return 0.;
  }
}


double fiberCorr(double r, double theta, double phi){
  return (1.+Pk_powerlaw_truncated_xi(r))*fiberCollision(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), 0.0254, 0.6087);
}


int MonteCarlo_SSPOC(double midr, double dr, double* mono, double* monoerror, double* fiber, double* fibererror, double* quad, double* quaderror){
  
  *monoerror       = pow(10., 6.);
  *mono = 0.;
  
  double monocumulative  = 0.;
  double mono2cumulative = 0.;

  *fibererror = pow(10., 6.);
  *fiber = 0.;

  double fibercumulative = 0.;
  double fiber2cumulative = 0.;

  *quad = 0.;
  *quaderror = pow(10., 6.);

  double quadcumulative = 0.;
  double quad2cumulative = 0.;

  int    daughters       = 0;
 
  double r, mu, phi, theta;
  
  double rx, ry;

  double shellVol = pow(midr, 2.)*dr*4.*pi;

  while((*monoerror/fabs(*mono) > 0.05) || (*fibererror/fabs(*fiber) > 0.05) || (*quaderror/(*quad) > 0.05)){
    for(j=0; j<100000; j++){
      r                = midr + dr*(gsl_rng_uniform(gsl_ran_r) - 0.5);

      phi              = 2.*pi*gsl_rng_uniform(gsl_ran_r);

      theta            = pi*gsl_rng_uniform(gsl_ran_r);

      mu               = cos(theta);

      monocumulative  += 1. +Pk_powerlaw_truncated_xi(r);
    
      mono2cumulative += pow(1. + Pk_powerlaw_truncated_xi(r), 2.);

      fibercumulative += fiberCorr(r, theta, phi);

      fiber2cumulative += pow(fiberCorr(r, theta, phi), 2.);

      quadcumulative  +=  fiberCorr(r, theta, phi)*0.5*(3.*pow(mu, 2.) - 1.);

      quad2cumulative += pow( fiberCorr(r, theta, phi)*0.5*(3.*pow(mu, 2.) - 1.), 2.);

      daughters       += 1;
    }

    *mono    = monocumulative/daughters;

    *fiber   = fibercumulative/daughters;

    *quad    = quadcumulative/daughters;
  
    *monoerror  = sqrt((mono2cumulative/daughters - pow(*mono, 2.))/daughters);

    *fibererror = sqrt((fiber2cumulative/daughters - pow(*fiber, 2.))/daughters);

    *quaderror  = sqrt((quad2cumulative/daughters - pow(*quad, 2.))/daughters);
  }

  return 0;
}


int fprintf_fiberCollision(){
  sprintf(filepath, "%s/Data/SpectralDistortion/toyCorrelationfn_fiberCollisionCorrection.dat", root_dir);                     
  double mono, monoerror, fiber, fibererror, quad, quaderror;                                                                 
  output = fopen(filepath, "w");                                                                                              
  for(jj=1; jj<40; jj++){                                                                                                      printf("\n %d", jj);                                                                                                        
   MonteCarlo_SSPOC(jj*0.25, 0.05, &mono, &monoerror, &fiber, &fibererror, &quad, &quaderror);                                  
   fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", jj*0.25, mono, monoerror, 1. + Pk_powerlaw_truncated_xi(jj*0.25), fiber, fibererror, quad, quaderror);                                                                               }                                                                                                                           
  for(jj=10; jj<150; jj++){                                                                                                    
    printf("\n %d", jj);                                                                                                       
    MonteCarlo_SSPOC(jj*1., 0.1, &mono, &monoerror, &fiber, &fibererror, &quad, &quaderror);                                    
    fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", jj*1., mono, monoerror, 1. + Pk_powerlaw_truncated_xi(jj*1.), fiber, fibererror, quad, quaderror);                                                                               
  }                                                                                                                            
  fclose(output);

  return 0;
}



