int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][0])/CellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/CellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][2])/CellSize);
    
    boxlabel                = (int)                        xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){
    double chi;

    for(j=0; j<n0*n1*n2; j++)  densityArray[j]     = 0.0;
    for(j=0; j<n0*n1*n2; j++)  meanCellRedshift[j] = 0.0;

    for(j=0; j<Vipers_Num; j++){
        // magnitude and redshift selection. 
        if(Acceptanceflag[j] == true){ 
            boxlabel                    = boxCoordinates(xCoor, yCoor, zCoor, j);
           
            densityArray[boxlabel]     += 1;
		    
		    meanCellRedshift[boxlabel] += zUtilized[j];
        }
    }
    
    for(j=0; j<n0*n1*n2; j++){
        if(densityArray[j] > 0.0){
	        // Currently densityArray contains solely galaxy counts per cell.
            meanCellRedshift[j] /= densityArray[j];
        
            chi                  = interp_comovingDistance(meanCellRedshift[j]);
        
            densityArray[j]     /= CellVolume*interp_nz(chi);
            
            densityArray[j]     -= 1.;  
        
            // printf("\n%e \t %e", meanCellRedshift[j], densityArray[j]);
        }
        
        // no galaxies. delta = -1. zeros handled by mask. 
        else densityArray[j]     = -1.;
    
        densityArray[j]         *= Cell_AppliedWindowFn[j];
    }
    
    printf("\n\noverdensity, min: %e, max: %e", arrayMin(densityArray, n0*n1*n2), arrayMax(densityArray, n0*n1*n2));
    
    return 0;
}


// outdated.
double CubeMeanNumberDensity(double chi){
    // Galaxies in cube/comoving volume. 
    
    return Vipers_Num/TotalVolume; 
}
