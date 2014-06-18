float Age_Integrand(float x){
    float  Interim_yVal;
    
    Interim_yVal = pow(1.0 - Om_m + Om_m*pow(e, -3.0*x) + Om_r*pow(e, -4.0*x) - (Om_tot -1.0)*pow(e, -2.0*x), -0.5);
    
    return Interim_yVal;
}


float f_Om_545(float x){
    float InterimInterp_yVal;
    
    splint(lnAarray, f_Om_mOfa545, f_Om_mOfa545_2derivs, linearGrowth_nPoints, x, &InterimInterp_yVal);
   
    return InterimInterp_yVal;
}


int UniverseAge(){
    H_0             = 100.0*h;                                                            // km s^-1 Mpc^-1 
    
    H_0inPerSec     = H_0*pow(3.08567758, -1.)*pow(10.0, -19.0);                          // s^{-1}                    
               
    HubbleTime      = pow(H_0inPerSec, -1.)*pow(OneYearInSeconds, -1.)*pow(10., -9);      // Hubble time in Gyears

    pt2AgeIntegrand = &Age_Integrand;
    
    
    linearGrowth_nPoints        = 14000;
    
    AgeOftheUniverse            = realloc(AgeOftheUniverse,             (linearGrowth_nPoints+1)*sizeof(*AgeOftheUniverse));
    Age2derivatives             = realloc(Age2derivatives,              (linearGrowth_nPoints+1)*sizeof(*Age2derivatives));

    lnAarray                    = realloc(lnAarray,                     (linearGrowth_nPoints+1)*sizeof(*lnAarray));

    HubbleCnstWithTime          = realloc(HubbleCnstWithTime,           (linearGrowth_nPoints+1)*sizeof(*HubbleCnstWithTime));
    HubbleConstant2derivatives  = realloc(HubbleConstant2derivatives,   (linearGrowth_nPoints+1)*sizeof(*HubbleConstant2derivatives));

    Om_mOfa                     = realloc(Om_mOfa,                      (linearGrowth_nPoints+1)*sizeof(*Om_mOfa));
    f_Om_mOfa545                = realloc(f_Om_mOfa545,                 (linearGrowth_nPoints+1)*sizeof(*f_Om_mOfa545));
    f_Om_mOfa545_2derivs        = realloc(f_Om_mOfa545_2derivs,         (linearGrowth_nPoints+1)*sizeof(*f_Om_mOfa545_2derivs));
    
    linear_growthfactor         = realloc(linear_growthfactor,          (linearGrowth_nPoints+1)*sizeof(*linear_growthfactor));
    approx2linear_growthfactor  = realloc(approx2linear_growthfactor,   (linearGrowth_nPoints+1)*sizeof(*approx2linear_growthfactor));
    SplineParams_ofgrowthfactor = realloc(SplineParams_ofgrowthfactor,  (linearGrowth_nPoints+1)*sizeof(*SplineParams_ofgrowthfactor));
    
    for(i=0; i<linearGrowth_nPoints+1; i++){
        lnAarray[i]             =  -11.0 + pow(1000., -1.0)*i;
    
        AgeOftheUniverse[i]     =          qromb(pt2AgeIntegrand, -13.5, lnAarray[i]);                  // In units of H_0
        
        HubbleCnstWithTime[i]   =  (float) H_0*pow(1. - Om_m + Om_m*pow(e, -3.0*lnAarray[i]), 0.5);     // In units of km s^-1 Mpc^-1
        
        Om_mOfa[i]              =          Om_m*pow(H_0/HubbleCnstWithTime[i], 2.0)*pow(e, -3.0*lnAarray[i]);
        
        f_Om_mOfa545[i]         =          pow(Om_mOfa[i], 0.545);
    
        // printf("\n%f \t %f \t %f \t %f", lnAarray[i], AgeOftheUniverse[i], HubbleCnstWithTime[i], Om_mOfa[i]);
    }

    spline(lnAarray, f_Om_mOfa545, linearGrowth_nPoints, 1.0e31, 1.0e31, f_Om_mOfa545_2derivs);

    printf("\nAge of the Universe today: %f Gyr\n\n", HubbleTime*qromb(pt2AgeIntegrand, -13.5, 0.0));  
    
    return 0;
}
