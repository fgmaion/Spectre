double splint_maskedRSD_pk(double k){
  // Interpolated matter power spectrum evaluated at mod(k_vec - q_vec).
  if(k>10.)          return pk_hiA*pow(k, 3. + pk_hin);

  else if(k<0.001)   return pk_loA*pow(k, 3. + pk_lon);

  else{
    double Interim;

    splint(sdltk, sdltPk, sdlt2d, pk_lineNo, k, &Interim);

    return Interim;
  }
}


int set_maskedRSDpaper_pk(){
  sprintf(filepath, "%s/Data/HODTheoryPk/outdated_cosmology/cambExtendedPk_hod_20.00.dat", root_dir);

  inputfile  = fopen(filepath, "r");

  ch         = 0;
  pk_lineNo  = 0;

  do{
    ch = fgetc(inputfile);
    if(ch == '\n')
      pk_lineNo += 1;
  } while (ch != EOF);

  rewind(inputfile);

  // Interpolated theoretical P(k) on a regular grid in k.
  sdltk  = realloc(sdltk,  pk_lineNo*sizeof(*sdltk));
  sdltPk = realloc(sdltPk, pk_lineNo*sizeof(*sdltPk));

  // Second derivates of HOD P(k) for cubic spline.
  sdlt2d = realloc(sdlt2d,         pk_lineNo*sizeof(*sdlt2d));

  for(j=0; j<pk_lineNo; j++)       fscanf(inputfile, "%le \t %le \n", &sdltk[j], &sdltPk[j]);

  fclose(inputfile);

  spline(sdltk, sdltPk, pk_lineNo, 1.0e31, 1.0e31, sdlt2d);

  powerlaw_regression(pk_lineNo, 0.0001, 0.001, 1., sdltk, sdltPk, &pk_loA, &pk_lon);

  powerlaw_regression(pk_lineNo,    8.,  10., 1., sdltk, sdltPk, &pk_hiA, &pk_hin);

  pt2Pk = &splint_maskedRSD_pk;

  return 0;
}
