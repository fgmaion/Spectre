int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][0])/CellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/CellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][2])/CellSize);
    
    boxlabel                = (int)                        xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){
    printf("\n\nCalculating overdensity...");

    double  chi;    
    double  app_mean;

    // assign memory for overdensity.
    prep_grid();

    for(j=0; j<n0*n1*n2; j++)           overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)           overdensity[j][1] = 0.0;

    for(j=0; j<Vipers_Num; j++){
        // magnitude and redshift selection. 
        if(Acceptanceflag[j] == true){ 
            boxlabel                    = boxCoordinates(xCoor, yCoor, zCoor, j);
            
            overdensity[boxlabel][0]   += 1;
            // overdensity[boxlabel][0]   += sampling[j];
            
            // cic_assign(xCoor[j], yCoor[j], zCoor[j], gal_z[j], 1.0);
        }
    }
    
    double x, y, z;
    
    for(k=0; k<n0; k++){
        // printf("\n%.2lf percentage complete.", 100.*k/n0);
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
               Index = i + n2*j + n2*n1*k;
            
               if(overdensity[Index][0] > 0.0){ 
                 x = (i + 0.5)*CellSize;
                 y = (j + 0.5)*CellSize;
                 z = (k + 0.5)*CellSize;
               
                 chi   = invert_StefanoBasis(CentreRA, CentreDec, &x, &y, &z);
            
                 overdensity[Index][0]   /= CellVolume*interp_nz(chi);
             
                 overdensity[Index][0]   -= 1.;  
               }
               
               // no galaxies. delta = -1. zeros handled by mask. 
               else overdensity[Index][0]   = -1.; 
           }
        }
    }
    
    free_HOD();
    
    printf("\nOverdensity calculation complete.");
    
    return 0;
}

// outdated.
double CubeMeanNumberDensity(double chi){
    // Galaxies in cube/comoving volume. 
    
    return Vipers_Num/TotalVolume; 
}
