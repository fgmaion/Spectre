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

// Print model to file. 

int kaiser_Multipoles(){
    double waveNumber    = 0.0;
    
    sprintf(filepath, "%s/Data/Multipoles/KaiserMultipoles_beta_%.2f.dat", root_dir, beta);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<100000; j++){
        waveNumber = pow(10., -5.)*j;

        fprintf(output, "%le \t %le \t %le \t %le \t %le \n", waveNumber, (*pt2Pk)(waveNumber), (*pt2Pk)(waveNumber)*kaiser_Monofactor(waveNumber, beta), (*pt2Pk)(waveNumber)*kaiser_Quadfactor(waveNumber, beta), (*pt2Pk)(waveNumber)*kaiser_Hexfactor(waveNumber, beta));
        
    }
    
    fclose(output);
    
    return 0;
}
