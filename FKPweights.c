int ApplyFKPweights(){
    TotalFKPweight    = 0.0;

    int NonEmpty;

    float nP;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index               = k*n1*n2 + j*n2 + i;   
                
                if(Cell_VIPERSbools[Index] > 0.0001){
                    Chi                 = Cell_chiVIPERSsystem[Index]; 
                
                    nP                  = interp_nz(Chi)*fkpPk;
                
                    FKPweights[Index]   = nP*pow(1. + nP, -1.);
                
                    TotalFKPweight     += nP*pow(1. + nP, -1.);   

                    NonEmpty           += 1;
                } 
            }
        }
    }
    
    return 0;
}
