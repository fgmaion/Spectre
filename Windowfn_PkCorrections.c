int CalcWfCorrections(){  
    fkpWeightedVolume  = 0.0;
    fkpSqWeightsVolume = 0.0;
    
    // eqn. 16.120, pg 524., P(k) overestimated due to larger density of states. 
    for(j=0; j<n0*n1*n2; j++) fkpWeightedVolume  += Cell_AppliedWindowFn[j]*CellVolume;

    for(j=0; j<n0*n1*n2; j++) fkpSqWeightsVolume += pow(Cell_AppliedWindowFn[j], 2.)*CellVolume;
    
    printf("\n\nFKP     weighted volume:     %e    [Total Volume]", fkpWeightedVolume/TotalVolume);
    printf("\nFKP Sq. weighted volume:     %e    [Total Volume]", fkpSqWeightsVolume/TotalVolume);
    
    printf("\nShot noise correction:  %f", ShotNoise());
    
    return 0;
}
