double age_Integrand(double x, void* p){
  (void) p;
  
  return pow(1.0 - Om_m + Om_m*pow(e, -3.0*x) + Om_r*pow(e, -4.0*x) - (Om_tot -1.0)*pow(e, -2.0*x), -0.5);
}

double HubbleInterp(double x){
  double HubbleVal;

  splint(lnAarray, HubbleCnstWithTime, HubbleConstant2derivatives, linearGrowth_nPoints, x, &HubbleVal);

  return HubbleVal;
}

double AgeInterp(double x){
  double Interim;

  splint(lnAarray, AgeOftheUniverse, Age2derivatives, linearGrowth_nPoints, x, &Interim);

  return Interim;
}

double f_Om_545(double x, void* p){
  double Interim;

  (void) p;
  
  splint(lnAarray, f_Om_mOfa545, f_Om_mOfa545_2derivs, linearGrowth_nPoints, x, &Interim);
   
  return Interim;
}


int UniverseAge(){
  H_0                         = 100.0*h;                                                            // km s^-1 Mpc^-1 
  H_0inPerSec                 = H_0*pow(3.08567758, -1.)*pow(10.0, -19.0);                          // s^{-1}                    
  HubbleTime                  = pow(H_0inPerSec, -1.)*pow(OneYearInSeconds, -1.)*pow(10., -9.);     // Hubble time in Gyears
    
  linearGrowth_nPoints        = 1120;

  lnAarray                    = realloc(lnAarray,                    linearGrowth_nPoints*sizeof(*lnAarray));
 
  AgeOftheUniverse            = realloc(AgeOftheUniverse,            linearGrowth_nPoints*sizeof(*AgeOftheUniverse));
  Age2derivatives             = realloc(Age2derivatives,             linearGrowth_nPoints*sizeof(*Age2derivatives));

  HubbleCnstWithTime          = realloc(HubbleCnstWithTime,          linearGrowth_nPoints*sizeof(*HubbleCnstWithTime));
  HubbleConstant2derivatives  = realloc(HubbleConstant2derivatives,  linearGrowth_nPoints*sizeof(*HubbleConstant2derivatives));

  Om_mOfa                     = realloc(Om_mOfa,                     linearGrowth_nPoints*sizeof(*Om_mOfa));
  f_Om_mOfa545                = realloc(f_Om_mOfa545,                linearGrowth_nPoints*sizeof(*f_Om_mOfa545));
  f_Om_mOfa545_2derivs        = realloc(f_Om_mOfa545_2derivs,        linearGrowth_nPoints*sizeof(*f_Om_mOfa545_2derivs));

  linear_growthfactor         = realloc(linear_growthfactor,         linearGrowth_nPoints*sizeof(*linear_growthfactor));
  approx2linear_growthfactor  = realloc(approx2linear_growthfactor,  linearGrowth_nPoints*sizeof(*approx2linear_growthfactor));
  SplineParams_ofgrowthfactor = realloc(SplineParams_ofgrowthfactor, linearGrowth_nPoints*sizeof(*SplineParams_ofgrowthfactor));
  
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double  result, error;
  gsl_function F;

  F.function = &age_Integrand;
    
  for(i=0; i<linearGrowth_nPoints; i++){
    lnAarray[i]             =  -11.0 + pow(100., -1.0)*i;

    gsl_integration_qags(&F, -13.5, lnAarray[i], 0, 1e-7, 1000, w, &result, &error);
    
    AgeOftheUniverse[i]     =   result;                 // In units of H_0
    
    // matter and vacuum energy. 
    HubbleCnstWithTime[i]   =   H_0*pow(1. - Om_m + Om_m*pow(e, -3.0*lnAarray[i]), 0.5);     // In units of km s^-1 Mpc^-1
    
    Om_mOfa[i]              =   Om_m*pow(H_0/HubbleCnstWithTime[i], 2.0)*pow(e, -3.0*lnAarray[i]);
        
    // GR
    f_Om_mOfa545[i]         =   pow(Om_mOfa[i], 0.545);
    
    // if(i%10 == 0)  printf("\n%d \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf", i, lnAarray[i], HubbleTime*AgeOftheUniverse[i], HubbleCnstWithTime[i], Om_mOfa[i], f_Om_mOfa545[i]);
  }

  spline(lnAarray, f_Om_mOfa545, linearGrowth_nPoints, 1.0e31, 1.0e31, f_Om_mOfa545_2derivs);

  gsl_integration_qags(&F, -13.5, 0.0, 0, 1e-7, 1000, w, &result, &error);
  
  printf("Age of the Universe today: %.3lf Gyr\n", HubbleTime*result);  

  spline(lnAarray,   AgeOftheUniverse,  linearGrowth_nPoints, 1.0e31, 1.0e31,            Age2derivatives);
  spline(lnAarray, HubbleCnstWithTime,  linearGrowth_nPoints, 1.0e31, 1.0e31, HubbleConstant2derivatives); // Natural cubic spline of H(a), zero second   derivative at both endpoints.
  
  return 0;
}
