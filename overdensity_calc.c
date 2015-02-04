int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][0])/CellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/CellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][2])/CellSize);
    
    boxlabel                = (int)                        xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){
    double  chi;    
    double* meanCellRedshift;
    double  app_mean;
    
    meanCellRedshift = (double *) malloc(n0*n1*n2*sizeof(*meanCellRedshift));

    // assign memory for overdensity.
    prep_grid();

    for(j=0; j<n0*n1*n2; j++)           overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)           overdensity[j][1] = 0.0;
    
    for(j=0; j<n0*n1*n2; j++)      meanCellRedshift[j]    = 0.0;

    for(j=0; j<Vipers_Num; j++){
        // magnitude and redshift selection. 
        if(Acceptanceflag[j] == true){ 
            boxlabel                    = boxCoordinates(xCoor, yCoor, zCoor, j);
            
            overdensity[boxlabel][0]   += 1;
            
            meanCellRedshift[boxlabel] += gal_z[j];
            
            // cic_assign(xCoor[j], yCoor[j], zCoor[j], gal_z[j], 1.0);
        }
    }
    
    for(j=0; j<n0*n1*n2; j++){
        if(overdensity[j][0] > 0.0){
	        // Currently densityArray contains solely galaxy counts per cell.
            meanCellRedshift[j] /= overdensity[j][0];
        
            chi                  = interp_comovingDistance(meanCellRedshift[j]);
        
            overdensity[j][0]   /= CellVolume*interp_nz(chi);
            
            overdensity[j][0]   -= 1.;  
        }
        
        // no galaxies. delta = -1. zeros handled by mask. 
        else overdensity[j][0]   = -1.;
    }
    
    free(meanCellRedshift);
    
    free_HOD();
    
    return 0;
}

// outdated.
double CubeMeanNumberDensity(double chi){
    // Galaxies in cube/comoving volume. 
    
    return Vipers_Num/TotalVolume; 
}
