double HubbleInterp(double x){
    double InterimInterp_HubbleVal;
 
    splint(lnAarray, HubbleCnstWithTime, HubbleConstant2derivatives, linearGrowth_nPoints, x, &InterimInterp_HubbleVal);
 
    return InterimInterp_HubbleVal;
}

/*
float AgeInterp(float x){
    float InterimInterp_TimeVal;
 
    splint(lnAarray, AgeOftheUniverse, Age2derivatives, linearGrowth_nPoints, x, &InterimInterp_TimeVal);
 
    return InterimInterp_TimeVal;
}*/


double linearGrowth_factor(double x){ 
    double Interim;
  
    splint(lnAarray, approx2linear_growthfactor, SplineParams_ofgrowthfactor, linearGrowth_nPoints, x, &Interim);
  
    return Interim;
}

/*
void  derivs(float x, float y[], float dydx[]){
// Input are the dependent variable vector y[1..n] and its derivative dydx[1..n] at the starting value of the independent variable x.
// Test case:                 Delta = (t/t_0)^(2/3)
//
// differential eqn. :        d2Delta/dt2 + t^(-1)*dDelta/dt + t^(-2)*Delta = (13/9)*t^(-4/3)*t_0^(-2/3)
//                  
//                            delta(t_0) = 1.
//                            dDelta/dt = alpha(t_0) = (2/3)*t_0^(-1)

//  dydx[1] = y[2];
//  dydx[2] = -1.*pow(x, -1.)*y[2] - 1.*pow(x,-2.)*y[1]+(13./9.)*pow(x,-4./3.)*pow(InitialStartTime, -2./3.);}  

  float Interim_HubbleParameter;
 
  Interim_HubbleParameter = (float) HubbleInterp(x);

  dydx[1] = y[2]*pow(Interim_HubbleParameter, -1.0);
  dydx[2] = -2.0*y[2] + (3.0/2.0)*pow(H_0, 2.0)*Om_m*y[1]*pow(Interim_HubbleParameter, -1.0)*pow(e, -3.0*x);
}*/


int linearGrowthRate(){
    /*
    #include  "/disk1/mjw/Aux_functions/RungeKutta_stepper.c"
    #include  "/disk1/mjw/Aux_functions/RungeKutta_driver.c"
    
    pt2derivs          = &derivs;

    pt2rkqs            = &rkqs;
    
    pt2f_Om_545        = &f_Om_545;
    
    yStartArray        = realloc(yStartArray,        (nVar+1)*sizeof(*yStartArray));
    yFinalArray        = realloc(yFinalArray,        (nVar+1)*sizeof(*yFinalArray));
    y2derivsStartArray = realloc(y2derivsStartArray, (nVar+1)*sizeof(*y2derivsStartArray));
    
    
    spline(lnAarray, AgeOftheUniverse,    linearGrowth_nPoints, 1.0e31, 1.0e31, Age2derivatives);

    spline(lnAarray, HubbleCnstWithTime,  linearGrowth_nPoints, 1.0e31, 1.0e31, HubbleConstant2derivatives); // Natural cubic spline of H(a), zero second   derivative at both endpoints.

    InitialStartTime      =  AgeOftheUniverse[10];                  
    Initial_lnScalefactor =          lnAarray[10];
    
    // delta                                                                    // y start vec is a 2 component column vector with components y[1] = g(InitialStartTime) = 1
    yStartArray[1]   = pow(10.0, -0.1);                                         // and y[2] = alpha(InitialStartTime) = dgdt(InitialStartTime) = (2./3.)*(1./InitialStartTime)
    
    yStartArray[2]   = yStartArray[1]*(2.0/3.0)*pow(InitialStartTime, -1.0);    // corresponding to assuming a matter dominated growing mode initial condition, delta(t) propto
                                                                                // pow(t, 2./3.).
                                                          
    derivs(Initial_lnScalefactor, yStartArray, y2derivsStartArray);
    
    // yFinalArray initialised to yStartArray as the variables are updated in place. 
    for(j=1; j<=nVar; j++) yFinalArray[j] = yStartArray[j];
    
    odeint(yFinalArray, nVar, Initial_lnScalefactor, 0.0, eps, defaultStepSize, MinAllowedStepSize, &nok, &nbad, pt2derivs, pt2rkqs);
    
    // arbitrary normalisation for growth factor. 
    growthfactor_today        = yFinalArray[1];
    */

    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double error, result;

    double alpha = 1.0;

    gsl_function F;

    F.function = &f_Om_545;

    gsl_integration_qags(&F, -9.5, 0.0, 0, 1e-7, 1000, w, &result, &error);

    // approx_growthfactor_today = exp(qromb(&f_Om_545, -9.5, 0.0));
    approx_growthfactor_today = exp(result);    

    // printf("\n\nLinear growth factor. \n");
    
    for(k=151; k<linearGrowth_nPoints; k++){ 
        // for(j=1; j<=nVar; j++) yFinalArray[j] = yStartArray[j];
            
        // odeint(yFinalArray, nVar, Initial_lnScalefactor, lnAarray[k], eps, defaultStepSize, MinAllowedStepSize, &nok, &nbad, pt2derivs, pt2rkqs);
    
        // linear_growthfactor[k]             = yFinalArray[1]/growthfactor_today;
        
        // approx2linear_growthfactor[k]         = qromb(&f_Om_545, -9.5, lnAarray[k]);
        // approx2linear_growthfactor[k]         = exp(approx2linear_growthfactor[k])/approx_growthfactor_today;
        
        gsl_integration_qags(&F, -9.5, lnAarray[k], 0, 1e-7, 1000, w, &result, &error);

        approx2linear_growthfactor[k]         = exp(result)/approx_growthfactor_today;

        // printf("%d \t %.4e \t %.4e \t %.4e\n", k, 1./exp(lnAarray[k]) - 1., f_Om_545(lnAarray[k]), approx2linear_growthfactor[k]);        
        // printf("%.4e \t %.4e \t %.4e \t %.4e \n", lnAarray[k], AgeInterp(lnAarray[k]), linear_growthfactor[k], approx2linear_growthfactor[k]);
    }
    
    spline(lnAarray, approx2linear_growthfactor, linearGrowth_nPoints, 1.0e31, 1.0e31, SplineParams_ofgrowthfactor);
    
    // gsl_integration_qags(&F, -9.5, -log(1.325), 0, 1e-7, 1000, w, &result, &error);

    // printf("\n\n%.4e \t %.4e \n", 0.325, exp(result)/approx_growthfactor_today);

    // Runge-Kutta driver with adaptive stepsize control.  Integrate starting values ystart[1..nvar] from x1 to x2 with accuracy eps, storing intermediate
    // results in global variables.  h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize(can be zero). On output nok and nbad
    // are the number of good and bad(but retried and fixed) steps taken, and ystart is replaced by values at the end of the integration interval. derivs is 
    // the user-supplied routine for calculating the right-hand side derivative, while rkqs is the name of the stepper routine to be used.

    // On output nok and nbad are the number of good and bad(but retried and fixed) steps taken
    
    return 0;
}
