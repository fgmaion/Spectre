int calc_clippingweights(){
    // ** Calculate clippling weights. Must have no folding. **//  
    if(Jenkins_foldfactor > 1.){  printf("\n\nError in clipping weights calc., set Jenkins fold factor to 1."); return 1;}
    
    double                       chi;
    double    cell_weights[n0*n1*n2];
    
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;
    
    // initiliase to no clipping. 
    for(j=0; j<n0*n1*n2; j++)  cell_weights[j] = 1.0;
    
    // assign galaxies using NGP.  overdensity contains the galaxy counts. 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j] == true){ 
            chi                               =                          interp_comovingDistance(zobs[j]);
        
            boxlabel                          =                    boxCoordinates(xCoor, yCoor, zCoor, j);
            
            // Add (1./sampling) in units of nbar. 
            overdensity[boxlabel][0]         +=                      pow(interp_nz(chi)*sampling[j], -1.);
        }
    }
    
    for(j=0; j<n0*n1*n2; j++){
        overdensity[j][0]                    /=                                                CellVolume;
        
        overdensity[j][0]                    -=                                                       1.0;
    }
    
    
    // True/false flag for enforcing zero mean for the filtered field.
    Gaussian_filter(clipping_smoothing_radius, 0);
   
    
    for(j=0; j<n0*n1*n2; j++){           
        if(smooth_overdensity[j][0] > appliedClippingThreshold){
            cell_weights[j]                   = (1. + appliedClippingThreshold)/(1. + overdensity[j][0]);
        
            // clipped volume fraction.
            fraction_clipped                 += 1.0;
        }
    }
    
    fraction_clipped /= n0*n1*n2;
    
    printf("\n\napplied clipping threshold: %lf, fraction of cells unclipped: %lf", appliedClippingThreshold, 1. - fraction_clipped);
    
    //  Clipping weights should not be renormalised, they simply multiply the (fkp_weight/sampling) without further renormalisation. 
    for(j=0; j<Vipers_Num; j++){
        if(Acceptanceflag[j]  == true){ 
            boxlabel           = boxCoordinates(xCoor, yCoor, zCoor, j);

            clip_galweight[j]  =                 cell_weights[boxlabel];
        }
        
        else{clip_galweight[j] = 0.0;}
    }

    return 0;
}
