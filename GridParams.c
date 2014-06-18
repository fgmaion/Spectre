int EvaluateGridParameters(){
    // Galaxy distribution plus padding.
    n0                    = (int) floor((AxisLimsArray[1][2] - AxisLimsArray[0][2])/CellSize);
    n1                    = (int) floor((AxisLimsArray[1][1] - AxisLimsArray[0][1])/CellSize);
    n2                    = (int) floor((AxisLimsArray[1][0] - AxisLimsArray[0][0])/CellSize);

    // Dimensions should be divisible by 2. 
    if(n0%2 == 1)  n0    += 1;
    if(n1%2 == 1)  n1    += 1;
    if(n2%2 == 1)  n2    += 1;

    printf("\n\nDimensions:  %d \t %d \t %d", n0, n1, n2);

    CellVolume            = pow(CellSize, 3);                                         // h^-3 Mpc^3
    TotalVolume           = n0*n1*n2*CellVolume;                                      // h^-3 Mpc^3
    
    clippedVolume         = 0.0;

    LowerChiLimit         = interp_comovingDistance(redshiftLowLimit);
    UpperChiLimit         = interp_comovingDistance(redshiftHiLimit);
    
    printf("\nRedshift limits, lower bound: %e \t %e, \n\t\t upper bound: %e \t %e", redshiftLowLimit, LowerChiLimit, redshiftHiLimit, UpperChiLimit);
    
    // FFTw calc assignment.
    xNyquistIndex         = n2/2 + 1;
    yNyquistIndex         = n1/2 + 1;
    zNyquistIndex         = n0/2 + 1;

    NyquistWaveNumber     = pi/CellSize;                                               // k = 2*pi x Nyquist frequency       
        
    kIntervalx            = 2.*pi*pow(n2, -1)*pow(CellSize, -1);
    kIntervaly            = 2.*pi*pow(n1, -1)*pow(CellSize, -1);
    kIntervalz            = 2.*pi*pow(n0, -1)*pow(CellSize, -1);
    
    return 0;
}


int assignbinninginterval(double interval){
    kbinInterval      = interval;
    
    kBinNumb          =  (int) ceil(modkMax/kbinInterval);   
    
    loskBinNumb       =  (int) ceil(modkMax/kbinInterval); 
    
    perpkBinNumb      =  loskBinNumb;

    printf("\n\nkbinterval:  %f", kbinInterval);
    printf("\nNumber of k bins:  %d", kBinNumb);
    
    printf("\nLos bin number:  %d", loskBinNumb);
    printf("\nPerp bin number: %d", perpkBinNumb);

    return 0;
}


double SolidAngleCalc(double decLowerBound, double decUpperBound, double raInterval){
    double SolidAngle    = 0.0;
    
    double thetaLowerRad = pi/2. - decUpperBound*pi/180.;
    double thetaUpperRad = pi/2. - decLowerBound*pi/180.;
    
    SolidAngle = (raInterval*pi/180.)*(-cos(thetaUpperRad) + cos(thetaLowerRad));

    printf("\nSolid angle observed in a given mock: %e steradians.", SolidAngle);
    printf("\n                                    : %e sq degrees.", steradians2sqdegs(SolidAngle));
 
    return SolidAngle;
}


double steradians2sqdegs(double inSteradians){
    return pow(180./pi, 2.)*inSteradians;
}


double sqdegs2steradians(double inSqdegs){
  return pow(pi/180., 2.)*inSqdegs;
}
