// Multipoles for the Kaiser redshift space distortion model.
double kaiser_Monofactor(double k, double beta){
  return 1. + (2./3.)*beta + 0.2*beta*beta;
}


double kaiser_Quadfactor(double k, double beta){
  return (4./3.)*beta + (4./7.)*beta*beta;
}


double kaiser_Hexfactor(double k, double beta){
  return (8./35.)*beta*beta;
}


double kaiser_multipole(double k, double beta, int monoQuad){
  switch(monoQuad){
    case 0:
      return kaiser_Monofactor(k, beta);  // Mono, L_0 corresponds to 0. Quad, L_2 corresponds to 2, Hex, L_4 corresponds to 4.
    case 2:
          return kaiser_Quadfactor(k, beta);
    case 4:
      return kaiser_Hexfactor(k, beta);
    }
}


double kaiser_multipole_xifactors(double r, double gamma, int monoQuad){
  switch(monoQuad){
    case 0:
      return 1.0;
    case 2:
      return -gamma*pow(3.-gamma, -1.);
    case 4:
      return gamma*(2.+gamma)*pow(-5.+gamma, -1.)*pow(-3.+gamma, -1.);
  }
}


int kaiser_Multipoles(){
  double beta;
  double waveNumber = 0.0;

  beta = fsigma8/bsigma8;
    
  sprintf(filepath, "%s/Data/Multipoles/KaiserMultipoles_beta_%.2f.dat", root_dir, beta);
    
  output = fopen(filepath, "w");
    
  for(j=0; j<100000; j++){
    waveNumber = pow(10., -5.)*j;

    fprintf(output, "%le \t %le \t %le \t %le \t %le \n", waveNumber, (*pt2Pk)(waveNumber), (*pt2Pk)(waveNumber)*kaiser_Monofactor(waveNumber, beta), (*pt2Pk)(waveNumber)*kaiser_Quadfactor(waveNumber, beta), (*pt2Pk)(waveNumber)*kaiser_Hexfactor(waveNumber, beta));  
  }
    
  fclose(output);
    
  return 0;
}
