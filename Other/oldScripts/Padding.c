int padVolume(double zfactor, double yfactor, double xfactor){
    double xextent = 0.;
    double yextent = 0.;
    double zextent = 0.;
    
    xextent = AxisLimsArray[1][0] - AxisLimsArray[0][0];
    yextent = AxisLimsArray[1][1] - AxisLimsArray[0][1];
    zextent = AxisLimsArray[1][2] - AxisLimsArray[0][2];
    
    // x dir
    AxisLimsArray[1][0] += 0.5*xfactor*xextent;
    AxisLimsArray[0][0] -= 0.5*xfactor*xextent;
    
    // y dir
    AxisLimsArray[1][1] += 0.5*yfactor*yextent;
    AxisLimsArray[0][1] -= 0.5*yfactor*yextent;
    
    // z dir 
    AxisLimsArray[1][2] += 0.5*zfactor*zextent;
    AxisLimsArray[0][2] -= 0.5*zfactor*zextent;
     
    return 0;
}
