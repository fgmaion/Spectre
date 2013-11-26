int EvaluateGridParameters(){
    n0                    = (int) floor((AxisLimsArray[1][2] - AxisLimsArray[0][2])/CellSize);
    n1                    = (int) floor((AxisLimsArray[1][1] - AxisLimsArray[0][1])/CellSize);
    n2                    = (int) floor((AxisLimsArray[1][0] - AxisLimsArray[0][0])/CellSize);

    // Dimensions should be divisible by 2. 
    if(n0%2 == 1)  n0    += 1;
    if(n1%2 == 1)  n1    += 1;
    if(n2%2 == 1)  n2    += 1;

    printf("\nDimensions:  %d \t %d \t %d", n0, n1, n2);

    CellVolume            = pow(CellSize, 3);                                         // h^-3 Mpc^3
    TotalVolume           = n0*n1*n2*CellVolume;                                      // h^-3 Mpc^3
    
    // FFTw calc assignment.
    xNyquistIndex         = n2/2 + 1;
    yNyquistIndex         = n1/2 + 1;
    zNyquistIndex         = n0/2 + 1;

    NyquistWaveNumber     = pi/CellSize;                                               // k = 2*pi x Nyquist frequency        
    
    kIntervalx            = 2*pi*pow(n2, -1)*pow(CellSize, -1);
    kIntervaly            = 2*pi*pow(n1, -1)*pow(CellSize, -1);
    kIntervalz            = 2*pi*pow(n0, -1)*pow(CellSize, -1);
    
    return 0;
}


int assignbinninginterval(float interval){
    kbinInterval = interval;
    
    kBinNumb     =  (int) ceil(modkMax/kbinInterval);   

    printf("\n\nkbinterval:  %f", kbinInterval);
    
    printf("\nNumber of k bins:  %d", kBinNumb);

    return 0;
}

float SolidAngleCalc(float decLowerBound, float decUpperBound, float raInterval){
    float SolidAngle;
    
    float thetaLowerRa = pi/2. - decUpperBound*pi/180.;
    float thetaUpperRa = pi/2. - decLowerBound*pi/180.;
    
    SolidAngle = (raInterval*pi/180.)*(-cos(thetaUpperRa) + cos(thetaLowerRa));

    printf("\nSolid angle observed in a given mock: %f", SolidAngle);

    return SolidAngle;
}
