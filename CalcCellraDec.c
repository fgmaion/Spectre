int CalcCellraDec(){
    double xCell, yCell, zCell;
    double rCell;

    double raCell, polarCell, decCell;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
	           xCell      = AxisLimsArray[0][0] + CellSize*(i+0.5);
	           yCell      = AxisLimsArray[0][1] + CellSize*(j+0.5);
	           zCell      = AxisLimsArray[0][2] + CellSize*(k+0.5);

               Index      = k*n1*n2 + j*n2 + i;

               Cell_rotatedXvals[Index] = xCell;
               Cell_rotatedYvals[Index] = yCell;
               Cell_rotatedZvals[Index] = zCell;
            }
        }
    }
    
    Celestialbasis(34.5, -5.10, Cell_rotatedXvals, Cell_rotatedYvals, Cell_rotatedZvals, n0*n1*n2);

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index       = k*n1*n2 + j*n2 + i;

                xCell       = Cell_rotatedXvals[Index];
                yCell       = Cell_rotatedYvals[Index];
                zCell       = Cell_rotatedZvals[Index]; 

                rCell       = pow(xCell*xCell + yCell*yCell + zCell*zCell, 0.5);     
		
                // acos returns radians, argument must be in the range -1 to 1. Returned value is between 0 and pi inclusive. 
                polarCell   = acos(zCell/rCell);

                // radians. 
                decCell     = pi/2. - polarCell;

                // returns the arc tangent in radians of y/x based on the signs of both values to determine the correct quadrant. both x and y cannot be zero, returned value is in the range !-pi to pi!
                raCell      = atan2(yCell, xCell);

                if(raCell < 0.0){
                    // right ascension in the range [0, 2*pi]
                    raCell += 2.*pi;
                }

                // conversion to degrees. 
                decCell    *= 180./pi;
                raCell     *= 180./pi;

                // right ascension in degrees.
	            Cell_raVIPERSsystem[Index]  = raCell;
		
                // declination in degrees.
                Cell_decVIPERSsystem[Index] = decCell;
                
                // Radial comoving distance in h^-1 Mpc. 
                Cell_chiVIPERSsystem[Index] = rCell;
            }
        }
    }
    
    sprintf(filepath, "%s/Data/ra_decCells/ra_dec_degs.dat", root_dir);
    output = fopen(filepath, "wb");
 
    fwrite(Cell_raVIPERSsystem, sizeof(double), n0*n1*n2, output);
    fwrite(Cell_decVIPERSsystem, sizeof(double), n0*n1*n2, output);
    fclose(output);

    /* // Read in the VIPERS weights. 
    
    sprintf(filepath, "%s/Data/ra_decCells/ra_dec_wght.dat", root_dir);
    inputfile = fopen(filepath, "rb");

    fread(Cell_raVIPERSsystem,  sizeof(double), n0*n1*n2, inputfile);
    fread(Cell_decVIPERSsystem, sizeof(double), n0*n1*n2, inputfile);
    fread(Cell_VIPERSweights,   sizeof(double), n0*n1*n2, inputfile);

    for(j=0; j<n0*n1*n2; j++){
        if(Cell_VIPERSweights[j] != 0.0){
             Cell_VIPERSbools[j]  = 1.00;
        }
    }
    
    fclose(inputfile);
    */
    
    for(j=0; j<n0*n1*n2; j++){
        if((LowerRAlimit<Cell_raVIPERSsystem[j]) && (Cell_raVIPERSsystem[j]<UpperRAlimit) && (LowerDecLimit<Cell_decVIPERSsystem[j]) && (Cell_decVIPERSsystem[j]<UpperDecLimit) 
            && (LowerChiLimit < Cell_chiVIPERSsystem[j]) && (Cell_chiVIPERSsystem[j] < UpperChiLimit)){
                Cell_SurveyLimitsMask[j] = 1.0;
        } 
    }
    
    // projectVIPERSsystem();

    return 0;
}


int projectVIPERSsystem(){
    sprintf(filepath, "%s/Data/ra_decCells/xyzwNonemptyCells.dat", root_dir);
    output = fopen(filepath, "w");

    double xCell, yCell, zCell;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index       = k*n1*n2 + j*n2 + i;

                xCell       = AxisLimsArray[0][0] + CellSize*(i+0.5);
                yCell       = AxisLimsArray[0][1] + CellSize*(j+0.5);
                zCell       = AxisLimsArray[0][2] + CellSize*(k+0.5);

                if(Cell_VIPERSweights[Index] > 0.0001){
                    fprintf(output, "%f \t %f \t %f \t %f\n", xCell, yCell, zCell, Cell_VIPERSweights[Index]);
                }
            }
        }
    }

    fclose(output);

    return 0;
}
