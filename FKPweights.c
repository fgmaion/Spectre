int ApplyFKPweights(double sampling){
    // Mean sampling rate of n bar, e.g. 100% or 40 %. 

    double nP;
    
    double minFKP  = 99.0;
    double maxFKP  =  0.0;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index                                    = k*n1*n2 + j*n2 + i;   

                Chi                                      = Cell_chiVIPERSsystem[Index]; 
                
                nP                                       = (*pt2nz)(Chi)*sampling*fkpPk;
                
                
                if(Cell_AppliedWindowFn[Index] > 0.){
                // Cells of the embedding volume close to the origin, chi << 1, cause problems when extapolating n(z) 
                // to these distances. Limit to only those in the survey.  
                 
                    Cell_AppliedWindowFn[Index]             *= nP*pow(1. + nP, -1.);
                    
                    if(nP*pow(1. + nP, -1.) > maxFKP) maxFKP = nP*pow(1. + nP, -1.);
                    if(nP*pow(1. + nP, -1.) < minFKP) minFKP = nP*pow(1. + nP, -1.);
                }
            }
        }
    }
    
    printf("\n\nMean sampl.: %.2e, FKP Pk: %.2e", meanSampling, fkpPk);
    printf("\nmin  FKP:    %.2e",    minFKP);
    printf("\nmax  FKP:    %.2e",    maxFKP);
    printf("\nvol. ratio:  %.2e", SumDoubleArray(Cell_AppliedWindowFn, n0*n1*n2)*CellVolume/TotalSurveyedVolume);
    
    return 0;
}
