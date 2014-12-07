double toyCorrfn(double r, double r0, double gamma){
  return pow(r/r0, -gamma);
}


double projectedCorrfn(double sig, double r0, double gamma){
  return sig*pow(r0/sig, gamma)*tgamma(0.5)*tgamma(0.5*(gamma-1.))/tgamma(gamma/2.);
}

int SaundersDeproject(int nbins, double interval){
  // References: google Nick Ross deprojection. Saunders 1992. 
   
  double*    w;
  double* mids;
  double*   xi;

  w    = malloc(nbins*sizeof(double));
  xi   = malloc(nbins*sizeof(double));
  mids = malloc(nbins*sizeof(double));

  for(j=0; j<nbins; j++){ 
    mids[j] = pow(10., -1. + interval*j);

    printf("\n%e", mids[j]);

       w[j] = projectedCorrfn(mids[j], 5.1, 1.57);
  
      xi[j] = 0.0;
  }

  sprintf(filepath, "%s/deprojection.dat", root_dir);

  output = fopen(filepath, "w");

  for(i=0; i<nbins; i++){
    for(j=i; j<nbins-1; j++){
      Interim = (mids[j+1]+sqrt(mids[j+1]*mids[j+1]-mids[i]*mids[i]))/(mids[j] + sqrt(mids[j]*mids[j] - mids[i]*mids[i]));     
       xi[i]  += (w[j+1]-w[j])*pow(mids[j+1]-mids[j],-1.)*log(Interim); 
    }

    xi[i]   *= -1./pi;
 
    fprintf(output, "\n%e \t %e \t %e \t %e", mids[i], w[i], xi[i], toyCorrfn(mids[i], 5.1, 1.57));
  }

  fclose(output);

  return 0;
}

