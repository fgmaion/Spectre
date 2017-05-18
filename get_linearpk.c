double splintLinearPk(double k){
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);

  else if(k<0.0001)  return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}

double linearPk_Gaussian(double k){
  return splintLinearPk(k)*exp(-0.5*pow(3.*k, 2.));
}

int inputLinearPk(){
  sprintf(model_flag, "linear");

  pt2Pk = &splintLinearPk;

  if((lo_zlim == 0.6) && (hi_zlim == 0.8))                          sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.593_z_0.75.dat", root_dir);
  else if ((lo_zlim == 0.6) && (hi_zlim == 0.9))                    sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.593_z_0.75.dat", root_dir);

  else if((lo_zlim == 0.8) && (hi_zlim == 1.0))                     sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.518_z_1.05.dat", root_dir);
  else if((lo_zlim == 0.9) && (hi_zlim == 1.2))                     sprintf(filepath, "%s/W1_Spectro_V7_2/pkmodels/linear_matter_pk_sig8_0.518_z_1.05.dat", root_dir);

  else{
    printf("\n\nError during input of real-space P(k) model.");

    exit(EXIT_FAILURE);
  }

  inputfile  = fopen(filepath, "r");

  line_count(inputfile, &pk_lineNo);

  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));   // Interpolated theoretical P(k) on a regular grid in k.
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));
  sdlt2d = realloc(sdlt2d, pk_lineNo*sizeof(*sdlt2d));  // Second derivates of HOD P(k) for cubic spline.

  for(j=0; j<pk_lineNo; j++)  fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);

  fclose(inputfile);

  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for sigma_8 calc.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  camb_sig8 = sigma8_calc();

  for(j=0; j<pk_lineNo; j++)  sdltPk[j] /= pow(camb_sig8, 2.); // normalised to sigma_8 of unity.
  /*                                                                                                                                                                                        
  for(j=0; j<pk_lineNo; j++){                                                                                                                                                                 
    if((0.02 < sdltk[j]) && (sdltk[j] < 0.4))  printf("\n%.6le \t %.6le", sdltk[j], sdltPk[j]);                                                                                               
  }                                                                                                                                                                                          
  */
  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon); // Add power laws for k<<1 and k>>1 for FFTlog calcs.
  powerlaw_regression(pk_lineNo,    8.0,  10.0, 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  // if((lo_zlim==0.6) && (hi_zlim<1.0))       app_sigma8 = 0.593139;
  // else if((lo_zlim==0.8) && (hi_zlim<1.3))  app_sigma8 = 0.518062;  // CHANGED FROM LOZ=0.9 AND app_sigma8 = 0.518 on 22 Mar 2017.

  // printf("\n\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 0.75, 0.593, f_Om_545(log(1./(1. + 0.75))), 0.593*f_Om_545(log(1./(1. + 0.75))));
  // printf("\nz_eff: %.2f.  sigma_8(z_eff): %.4f. f(z_eff): %.4f. f(z_eff)*sigma_8(z_eff): %.4f.", 1.05, 0.518, f_Om_545(log(1./(1. + 1.05))), 0.518*f_Om_545(log(1./(1. + 1.05))));

  return 0;
}
