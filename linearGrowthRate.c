#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

double linearGrowth_factor(double x){ 
  double Interim;
  
  splint(lnAarray, approx2linear_growthfactor, SplineParams_ofgrowthfactor, linearGrowth_nPoints, x, &Interim);
    
  return Interim;
}

int  derivs(double x, double y[], double dydx[], void* p){
  // Input are the dependent variable vector y[1..n] and its derivative dydx[1..n] at the starting value of the independent variable x.
  // Test case:                 Delta = (t/t_0)^(2/3)
  //
  // differential eqn. :        d2Delta/dt2 + t^(-1)*dDelta/dt + t^(-2)*Delta = (13/9)*t^(-4/3)*t_0^(-2/3)
  //                  
  //                            delta(t_0) = 1.
  //                            dDelta/dt = alpha(t_0) = (2/3)*t_0^(-1)

  //  dydx[1] = y[2]; 
  //  dydx[2] = -1.*pow(x, -1.)*y[2] - 1.*pow(x,-2.)*y[1]+(13./9.)*pow(x,-4./3.)*pow(InitialStartTime, -2./3.);}  

  (void) x; /* avoid unused parameter warning */
  (void) p;
  
  double HubParam = HubbleInterp(x);

  dydx[0] =      y[1]/HubParam;
  dydx[1] = -2.0*y[1] + (3.0/2.0)*pow(H_0, 2.0)*Om_m*y[0]*pow(HubParam, -1.0)*pow(e, -3.0*x);

  return GSL_SUCCESS;
}

int jac(double x, const double y[], double *dfdy, double dfdt[], void *params){
  (void)      x; /* avoid unused parameter warning */
  (void)      y;
  (void) params;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
  gsl_matrix* m            = &dfdy_mat.matrix;

  double HubParam = HubbleInterp(x);
  
  gsl_matrix_set(m, 0, 0, 0.0);
  gsl_matrix_set(m, 0, 1, 1.0/HubParam);

  gsl_matrix_set(m, 1, 0, (3.0/2.0)*pow(H_0, 2.0)*Om_m*pow(HubParam, -1.0)*pow(e, -3.0*x));
  gsl_matrix_set(m, 1, 1, -2.);

  dfdt[0] = 0.0;
  dfdt[1] = 0.0;

  return GSL_SUCCESS;
}

int linearGrowthRate(){
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  gsl_function null, F;
  double  result, error;

  null.function =      NULL;
  F.function    = &f_Om_545;
  
  // gsl_odeiv2_system sys = {derivs, jac, 2, jac};
  // gsl_odeiv2_driver* d  = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-3, 1e-8, 1e-8);

  int      s;
  void*    v;
  double y[2], ydash[2];

  double init_t, init_lna;                  

  init_t    = HubbleTime*AgeOftheUniverse[10];
  init_lna  =                    lnAarray[10];
    
  y[0]      = pow(10.0, -0.1);                          // y[1] = alpha(InitialStartTime) = dgdt(InitialStartTime) = (2./3.)*(1./InitialStartTime)
  y[1]      = y[0]*(2.0/3.0)*pow(init_t, -1.0);         // corresponding to assuming a matter dominated growing mode initial condition,
                                                        // delta(t) \propto pow(t, 2./3.).                                                          
  derivs(init_lna, y, ydash, v);
  
  // s = gsl_odeiv2_driver_apply(d, &init_lna, 0., y);
  
  // if(s != GSL_SUCCESS)  printf ("error: driver returned %d\n", s);
    
  growthfactor_today = y[0];                            // arbitrary normalisation for growth factor.

  // gsl_odeiv2_driver_reset(d);
  
  gsl_integration_qags(&F, -9.5, 0.0, 0, 1e-7, 1000, w, &result, &error);
  
  approx_growthfactor_today = exp(result);    
  
  // printf("\n\nln(a) \t\t Age \t\t H(z) \t\t Om(a) \t\t f(a) \t\t D+(a) \t\t approx. D+(a) \n");
  
  for(k=151; k<linearGrowth_nPoints; k++){
    init_lna  = lnAarray[10];

    y[0]      = pow(10.0, -0.1);                        // y[1] = alpha(InitialStartTime) = dgdt(InitialStartTime) = (2./3.)*(1./InitialStartTime)
    y[1]      = y[0]*(2.0/3.0)*pow(init_t, -1.0);
    
    // s         = gsl_odeiv2_driver_apply(d, &init_lna, lnAarray[k], y);    

    linear_growthfactor[k] = y[0]/growthfactor_today;
                
    gsl_integration_qags(&F, -9.5, lnAarray[k], 0, 1e-7, 1000, w, &result, &error);

    approx2linear_growthfactor[k] = exp(result)/approx_growthfactor_today;

    /*
    if((k > 1050) && (i%10 == 0)){
      printf("%+.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \t %.4lf \n", lnAarray[k], HubbleTime*AgeInterp(lnAarray[k]), HubbleCnstWithTime[k], Om_mOfa[k], f_Om_mOfa545[k],
                                                                                linear_growthfactor[k], approx2linear_growthfactor[k]);
    }
    */    

    // gsl_odeiv2_driver_reset(d);
  }
    
  spline(lnAarray, approx2linear_growthfactor, linearGrowth_nPoints, 1.0e31, 1.0e31, SplineParams_ofgrowthfactor);
    
  gsl_integration_qags(&F, -9.5, -log(1.325), 0, 1e-7, 1000, w, &result, &error);
  
  // printf("\n\nNormalised growth factor @ z=%.3lf is %.3lf \n", 0.325, exp(result)/approx_growthfactor_today);
  
  // gsl_odeiv2_driver_free(d);
  
  return 0;
}
