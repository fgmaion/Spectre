double nbar_dV(double chi){
  return pow(chi, 2.)*interp_nz(chi);
}


int prep_inverseCumulative_nbar(){
  // Calculate cumulative nbar and splint its inverse.
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double norm, error, result;

  double alpha = 1.0;

  gsl_function F;

  F.function = &nbar_dV;
  
  gsl_integration_qags(&F, loChi - 1., hiChi + 1., 0, 1e-7, 1000, w, &norm, &error);

  for(ii=0; ii<400; ii++){
    chi_cumulative_nbar[ii] = loChi - 1. + ii*(hiChi - loChi + 2.)/399.;

    // function, low limit, high limit, abs err., rel. err., max step number, workspace (memory), result, error.
    gsl_integration_qags(&F, loChi - 1., chi_cumulative_nbar[ii], 0, 1e-7, 1000, w, &result, &error);

    // cumulative_nbar[i]  = qromb(&nbar_dV, loChi-100., chi_cumulative_nbar[i]);
    cumulative_nbar[ii]     = result;

    cumulative_nbar[ii]    /= norm;

    printf("\n%.4lf \t %.4lf", chi_cumulative_nbar[ii], cumulative_nbar[ii]);
  }

  spline(cumulative_nbar, chi_cumulative_nbar, 400, 1.0e31, 1.0e31, cumulative_nbar_2d);

  return 0;
}


double inverse_cumulative_nbar(double arg){
  // if F is the cumulative distribution of nbar.
  // return the inverse of F, for random generation.

  double result;

  splint(cumulative_nbar, chi_cumulative_nbar, cumulative_nbar_2d,  400, arg, &result);

  return result;
}
