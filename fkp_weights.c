double calc_fkpweights(){
    double       nbar;
    
    double norm = 0.0;

    // FKP weights for randoms; sets normalisation 
    for(j=0; j<rand_number; j++){
        if(rand_accept[j]    == true){
            nbar              =   interp_nz(rand_chi[j]);
    
            rand_weight[j]    = 1./(1. + fkpPk*nbar);

            norm             += nbar*pow(rand_weight[j], 2.); 
        }
    }
    
    norm                     *= alpha;
    
    norm                      = sqrt(norm);
    
    for(j=0; j<rand_number; j++)  rand_weight[j] /= norm;
    
    printf("\n\nFKP weight normalisation: %.4lf", norm);
    
    return norm;
}


int set_gal_fkpweights(double norm){
    double chi;

    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){
	    chi               = interp_comovingDistance(zobs[j]); 
	    
            fkp_galweight[j]  = 1./(1. + fkpPk*interp_nz(chi));
	    
	    fkp_galweight[j] /= norm; 
        }
    }
       
    return 0;
}
