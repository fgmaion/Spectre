double splintLinearPk(double k){
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);

  else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}

double sigma8_integrand(double k, void* params){
  return pow(2.*pi, -3.)*splintLinearPk(k)*4.*pi*k*k*spherical_tophat(k, 8.)*spherical_tophat(k, 8.);
}

double sigma8_calc(){
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  double result, error;

  gsl_function F;

  F.function = &sigma8_integrand;

  gsl_integration_qags(&F, 0.0, 200., 0, 1e-7, 1000, w, &result, &error);

  // double precision handles %.15lf for initialisation, e.g. p = PI.
  printf("\n\nsig_8: %.10lf (GSL)", sqrt(result));

  return sqrt(result);
}

int get_linearsig8(){
  // e.g.  /home/mjw/CAMB/models/linear_sigma8_Om_cdm_0.262_Om_v_0.690_H0_67.30_z_0.607.dat
  sprintf(filepath, "/home/mjw/CAMB/models/linear_sigma8_Om_cdm_%.3lf_Om_v_%.3lf_H0_%.2lf_z_%.3lf.dat", Om_m - Om_b, Om_v, 100.*h, z_eff);

  while((inputfile = fopen(filepath, "r")) == NULL){
    camb_call(0, z_eff);      // nonlinear/linear flag, redshift;
  }

  fscanf(inputfile, "%lf", &camb_sig8);

  fclose(inputfile);

  // printf("\n\nLinear sigma_8 retrieved: %.4lf; z_eff: %.4lf", camb_sig8, z_eff);

  return 0;
}

double linearPk_Gaussian(double k){
  return splintLinearPk(k)*exp(-0.5*pow(3.*k, 2.));
}

int linear_pk(){
  char buffer[200];
  
  get_linearsig8();
  
  pt2Pk = &splintLinearPk;

  sprintf(model_flag, "linear");
  
  sprintf(filepath, "/home/mjw/CAMB/models/camb_matter_pk_linearity_0_Om_cdm_%.3lf_Om_v_%.3lf_H0_%.2lf_z_%.3lf.dat", Om_m - Om_b, Om_v, 100.*h, z_eff);

  printf("\nFor z_eff of %.3lf, linear sigma_8 is %.6lf, f sigma_8 is %.3lf, and loading: \n%s", z_eff, camb_sig8, print_fsig8(), filepath);

  while((inputfile = fopen(filepath, "r")) == NULL){
    camb_call(0, z_eff);      // nonlinear/linear flag, redshift;
  }

  linecount_header(inputfile, 1, &pk_lineNo);

  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
  sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.

  fgets(buffer, 200, inputfile);  

  for(j=0; j<pk_lineNo; j++)  fscanf(inputfile, "    %lE    %lE \n", &sdltk[j], &sdltPk[j]);

  fclose(inputfile);

  // camb_sig8 = sigma8_calc();

  for(j=0; j<pk_lineNo; j++)  sdltPk[j] /= pow(camb_sig8, 2.); // normalised to sigma_8 of unity.

  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  print_fsig8();  // needs camb sig_8(z).
  
  return 0;
}

double print_fsig8(){
  void*  params;
 
  double lna   = log(1./(1. + z_eff));
  
  double Dplus = linearGrowth_factor(lna);

  double fsig8 = camb_sig8*f_Om_545(lna, params); // camb returns sig8(z).

  return fsig8;
}
