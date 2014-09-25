double fiberCollision(double rx, double ry, double dx, double dy){
  if((fabs(rx)<dx) && (fabs(ry)<dy)){
    return 0.;
  }

  else{
    return 1.;
  }
}


double fiberCorr(double r, double theta, double phi){
  return (1.+ (*pt2Xi)(r))*fiberCollision(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), 0.0254, 0.6087);
}


double splint_realSpaceCorrfn(double r){
  double Interim;

  splint(rvals, xi, xi2D, splinexi_N, r, &Interim);

  return Interim;
}


int spline_realSpaceCorrfn(int N){
  splinexi_N = N;

  sprintf(filepath, "%s/MockAvg_realSpaceCorrfn.dat", root_dir);
  
  inputfile = fopen(filepath, "r");

  xi     = malloc(splinexi_N*sizeof(double));
  rvals  = malloc(splinexi_N*sizeof(double));
  xi2D   = malloc(splinexi_N*sizeof(double));

  for(j=0; j<splinexi_N; j++){  
    fscanf(inputfile, "%le \t %\le \n", &rvals[j], &xi[j]);

    // printf("\n%e \t %e", rvals[j], xi[j]);
  }
  
  fclose(inputfile);

  spline(rvals, xi, splinexi_N, 1.0e31, 1.0e31, xi2D);

  return 0;
}


int MonteCarlo_SSPOC(double midr, double dr, double* mono, double* monoerror, double* fiber, double* fibererror, double* quad, double* quaderror){
  
  *monoerror       = pow(10., 6.);
  *mono = 99.;
  
  double monocumulative  = 0.;
  double mono2cumulative = 0.;

  *fibererror = pow(10., 6.);
  *fiber = 99.;

  double fibercumulative = 0.;
  double fiber2cumulative = 0.;

  *quad = 99.;
  *quaderror = pow(10., 6.);

  double quadcumulative = 0.;
  double quad2cumulative = 0.;

  int    daughters       = 0;
 
  double r, mu, phi, theta;
  
  double rx, ry;

  double shellVol = pow(midr, 2.)*dr*4.*pi;

  while(*fibererror > 0.3){
  // while((*monoerror> 0.3) || (*fibererror > 0.3)){
    printf("\n another round, %e \t %e", *monoerror, *fibererror);
  // while((*monoerror/fabs(*mono) > 0.05) || (*fibererror/fabs(*fiber) > 0.05) || (*quaderror/(*quad) > 0.05)){
    for(j=0; j<1000000; j++){
      r                 = midr + dr*(gsl_rng_uniform(gsl_ran_r) - 0.5);

      phi               = 0.5*pi*gsl_rng_uniform(gsl_ran_r);

      theta             = 0.5*pi*gsl_rng_uniform(gsl_ran_r);

      mu                = cos(theta);

      // printf("\n%e \t %e \t %e \t %e \t %e \t %e", r, theta, phi, r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), fiberCorr(r, theta, phi));

      monocumulative   += 1. + (*pt2Xi)(r);
    
      mono2cumulative  += pow(1. + (*pt2Xi)(r), 2.);

      fibercumulative  += fiberCorr(r, theta, phi);

      fiber2cumulative += pow(fiberCorr(r, theta, phi), 2.);

      // quadcumulative  +=  fiberCorr(r, theta, phi)*0.5*(3.*pow(mu, 2.) - 1.);

      // quad2cumulative += pow( fiberCorr(r, theta, phi)*0.5*(3.*pow(mu, 2.) - 1.), 2.);

      daughters       += 1;
    }

    *mono    = monocumulative/daughters;

    *fiber   = fibercumulative/daughters;

    // *quad    = quadcumulative/daughters;
  
    *monoerror  = sqrt((mono2cumulative/daughters - pow(*mono, 2.))/daughters);

    *fibererror = sqrt((fiber2cumulative/daughters - pow(*fiber, 2.))/daughters);

    // printf("\n%e \t %e \t %e \t %e", *mono, *fiber, *monoerror, *fibererror);

    // *quaderror  = sqrt((quad2cumulative/daughters - pow(*quad, 2.))/daughters);
  }

  return 0;
}


int fprintf_fiberCollision(){
  pt2Xi   = &splint_realSpaceCorrfn;

  sprintf(filepath, "%s/splintRealSpaceCorrelationfn_fiberCollision.dat", root_dir);                     

  double mono, monoerror, fiber, fibererror, quad, quaderror;                                                                 

  output = fopen(filepath, "w");                                                                                      

  for(jj=0; jj<55; jj++){
    printf("\n %d", jj);                                                                                                       
    MonteCarlo_SSPOC(rvals[jj], 0.2, &mono, &monoerror, &fiber, &fibererror, &quad, &quaderror);                                    
    fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n", rvals[jj], mono, monoerror, 1. + (*pt2Xi)(rvals[jj]), fiber, fibererror, quad, quaderror);                                                                               
  }                                                                                                                            
  fclose(output);

  return 0;
}



