int ApplyFKPweights(){
    TotalFKPweight    = 0.0;
    
    int NonEmptyCells = 0;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index               = k*n1*n2 + j*n2 + i;   

                if(booldensity[Index]  != 0.0){
                    NonEmptyCells      += 1;
                    
                    Chi                 = CellSize*pow(k*k + j*j + i*i, 0.5); 
                
                    Interim             = interp_nz(Chi)*fkpPk;
                
                    FKPweights[Index]   = Interim*pow(1. + Interim, -1.);
                
                    TotalFKPweight     += Interim*pow(1. + Interim, -1.);    
                }
            }
        }
    }
    
    return 0;
}
