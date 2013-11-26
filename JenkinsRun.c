int JenkinsCoordinates(){
    // Jenkin's run to beat aliasing of P(k) near the Nyquist wavenumber. 
    AxisLimsArray[0][0]         *=     1.0;                                    // h^-1 Mpc
    AxisLimsArray[1][0]         /=     JenkinsScalefactor;                     // h^-1 Mpc

    AxisLimsArray[0][1]         *=     1.0;                                    // h^-1 Mpc
    AxisLimsArray[1][1]         /=     JenkinsScalefactor;                     // h^-1 Mpc

    AxisLimsArray[0][2]         *=     1.0;                                    // h^-1 Mpc
    AxisLimsArray[1][2]         /=     JenkinsScalefactor;                     // h^-1 Mpc
         
    CellSize                    /=     JenkinsScalefactor;                     // Cell size, comoving distance, h^-1 Mpc
    
    kbinInterval                *=     JenkinsScalefactor;
    return 0;
}


int JenkinsFold(float original[], int lenArray, int axis){      // Upper limit           // Lower limit
    for(j=0; j<lenArray; j++)  original[j] = fmod(original[j], (AxisLimsArray[1][axis] - AxisLimsArray[0][axis]));               

    return 0;
}


int ApplyJenkins(){
    // Jenkins run to beat aliasing. 
    JenkinsFold(xCoor, Vipers_Num, 0);
    JenkinsFold(yCoor, Vipers_Num, 1);
    JenkinsFold(zCoor, Vipers_Num, 2);  

    printf("\n\nAfter Jenkin's folding...");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));

    return 0;
}
