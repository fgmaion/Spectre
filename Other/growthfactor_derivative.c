#define NR_END 1
#define FREE_ARG char*

float ln_linearGrowth_factor(float x){
    Interim = linearGrowth_factor(x);
 
    return log(Interim);
}


float linearGrowth_factor_firstDeriv(float x){
    float Interim_InterpVal_deriv;
  
    splint(xVals, gdot, SplineParams_ofgdot, nDerivs, x, &Interim_InterpVal_deriv);
  
    return Interim_InterpVal_deriv;
}


float linearGrowth_factor_scndDeriv(float x){
    float Interim_InterpVal_2deriv;
  
    splint(xVals, Second_derivs, SplineParams_ofgdotdot, nDerivs, x, &Interim_InterpVal_2deriv);
  
    return Interim_InterpVal_2deriv;
}


int growthfactor_derivative(){
    linearGrowthRate();
    
    nDerivs                       = 8000;
    
    pt2lnlinearGrowthRate         = &ln_linearGrowth_factor;
    pt2linearGrowthRate_deriv     = &linearGrowth_factor_firstDeriv;
    pt2linearGrowthRate_2deriv    = &linearGrowth_factor_scndDeriv;
    
    gdot                          = realloc(gdot, (nDerivs+1)*sizeof(*gdot));        
    
    xVals                         = realloc(xVals, (nDerivs+1)*sizeof(*xVals));
    length_scales                 = realloc(length_scales, (nDerivs+1)*sizeof(*length_scales));
    derivs_error                  = realloc(derivs_error, (nDerivs+1)*sizeof(*derivs_error));
    
    Second_derivs                 = realloc(Second_derivs, (nDerivs+1)*sizeof(*Second_derivs));
    Second_derivs_error           = realloc(Second_derivs_error, (nDerivs+1)*sizeof(*Second_derivs_error));
    Second_deriv_lengthscales     = realloc(Second_deriv_lengthscales, (nDerivs+1)*sizeof(*Second_deriv_lengthscales));
    
    SplineParams_ofgdot           = realloc(SplineParams_ofgdot, (nDerivs+1)*sizeof(*SplineParams_ofgdot));
    SplineParams_ofgdotdot        = realloc(SplineParams_ofgdotdot, (nDerivs+1)*sizeof(* SplineParams_ofgdotdot));

    for(j=1; j<=nDerivs; j++){
        // Calculate the derivative at the last 4000 points of lnAarray.  In the vacuum dominated era, growth rate freezes, calculating derivative is difficult.
        // Shift back in time by 1000. 
        xVals[j]                     = lnAarray[linearGrowth_nPoints - nDerivs - 1000 + j];
  
        length_scales[j]             = 0.5;
        
        Second_deriv_lengthscales[j] = 0.5;
    }
    
    for(j=1; j<=nDerivs; j++){ 
        gdot[j]                   = dfridr(pt2lnlinearGrowthRate, xVals[j], length_scales[j], &derivs_error[j]);    
        
        for(;;){ 
            if(derivs_error[j] > 0.0001){
  	            length_scales[j] /= 2.0;
  	            
  	            gdot[j] = dfridr(pt2lnlinearGrowthRate, xVals[j], length_scales[j], &derivs_error[j]);}
            else{
                break;
            }
        }
    }
    
    spline(xVals, gdot, nDerivs, 1.0e31, 1.0e31, SplineParams_ofgdot);
    
    printf("\nxVals max:  %e", lnAarray[linearGrowth_nPoints] - 1000);
    
    printf("\nvalue of the growth factor today:  %f",                    linearGrowth_factor(0.0));
    printf("\n\nvalue of the growth factor at redshift 9:  %f",          linearGrowth_factor(-2.3));
    printf("\n\nfirst derivative of ln D_{+} wrt ln(a) today:  %f",      linearGrowth_factor_firstDeriv(0.0));
    printf("\n\nf today,    for a MB=-20. sample:  %f",                  linearGrowth_factor_firstDeriv(0.0));
    printf("\n\nf at z=0.8, for a MB=-20. sample:  %f",                  linearGrowth_factor_firstDeriv(-0.587786665));
    
    for(j=1;j<=nDerivs;j++){ 
        Second_derivs[j] = dfridr(pt2linearGrowthRate_deriv, xVals[j], Second_deriv_lengthscales[j], &Second_derivs_error[j]);
      
        for(;;){ 
            if(Second_derivs_error[j] > 0.0001){
	            Second_deriv_lengthscales[j] /= 2.0;
	            Second_derivs[j] = dfridr(pt2linearGrowthRate_deriv, xVals[j], Second_deriv_lengthscales[j], &Second_derivs_error[j]);
	        }

    	    else{
    	        break;
       	    }
        }
    }
    
    spline(xVals, Second_derivs, nDerivs, 1.0e31, 1.0e31, SplineParams_ofgdotdot);

    printf("\n\nsecond derivative of g wrt ln(a) today:  %f\n",   linearGrowth_factor_scndDeriv(0.0));
    
    sprintf(filepath, "%s/Data/zComovingDistance/f_sigma8.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    double fsigma8;
    
    for(k=linearGrowth_nPoints - nDerivs - 999; k<linearGrowth_nPoints - 1000; k++){
        fsigma8 = sigma_8*linear_growthfactor[k]*linearGrowth_factor_firstDeriv(lnAarray[k]);
    
        fprintf(output, "%f \t %f \t %f \t %f \t %f \t %f \t %f \n", lnAarray[k], linear_growthfactor[k], approx2linear_growthfactor[k], linearGrowth_factor_firstDeriv(lnAarray[k]), linearGrowth_factor_firstDeriv(lnAarray[k])/1.495903, f_Om_mOfa545[k], fsigma8);
    }
    
    fclose(output);

    return 0;
}
