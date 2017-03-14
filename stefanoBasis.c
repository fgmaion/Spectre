int StefanoBasis(int Num, double ra[], double dec[], double rDist[], double xCoor[], double yCoor[], double zCoor[]){
  printf("\n\nAngular limits of galaxies: %.4lf < ra < %.4lf, %.4lf < dec < %.4lf", AcceptedMin(ra,  Acceptanceflag, Vipers_Num),  AcceptedMax(ra, Acceptanceflag, Vipers_Num),
                                                                                    AcceptedMin(dec, Acceptanceflag, Vipers_Num), AcceptedMax(dec, Acceptanceflag, Vipers_Num));  
  // Reflection applied in coordinate calc.  
  // StefanoReflection(Vipers_Num, CentreRA, CentreDec, xCoor, yCoor, zCoor);
  
  // Rotate the input co-ordinates such that the ra direction is aligned more or less with the y axis, dec direction with x, and redshift along z.
  StefanoRotated(Vipers_Num, CentreRA, CentreDec, xCoor, yCoor, zCoor);
  
  printf("\n\nAccepted, inverted, rotated & translated");

  printf("\nx min:  %.3f \t x max:  %.3f", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
  printf("\ny min:  %.3f \t y max:  %.3f", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
  printf("\nz min:  %.3f \t z max:  %.3f", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));
  
  return 0;
}


double invert_StefanoBasis(double centreRA, double centreDec, double* xval, double* yval, double* zval){
  double x,  y,   z;
  double x1, y1, z1;
  double x2, y2, z2;

  double c_ra, c_dec;

  double rmod;

  x       = *xval;
  y       = *yval;
  z       = *zval;

  c_ra    =  centreRA*(pi/180.);
  c_dec   = centreDec*(pi/180.);

  // reverse translation.
  x  -= stefano_trans_x;
  y  -= stefano_trans_y;
  z  -= stefano_trans_z;

  // invert R2, x'' to x'
  x2         =  -sin(c_dec)*x + cos(c_dec)*z;
  y2         =              y;
  z2         =  -cos(c_dec)*x - sin(c_dec)*z;

  // invert R1, x' to x.
  x1         =  cos(c_ra)*x2 - sin(c_ra)*y2;
  y1         =  sin(c_ra)*x2 + cos(c_ra)*y2;
  z1         =  z2;

  // invert z inversion through xy plane.
  *xval      =  x1;
  *yval      =  y1;
  *zval      = -z1;

  rmod       = sqrt(x1*x1 + y1*y1 + z1*z1);

  return rmod;
}


int StefanoReflection(int Number, double centreRA, double centreDec, double xCoors[], double yCoors[], double zCoors[]){
  for(j=0; j<Number; j++)  zCoors[j] *= -1.0;  // inversion through the xy plane.
  
  return 0;
}


int StefanoRotated(int Number, double centreRA, double centreDec, double xCoors[], double yCoors[], double zCoors[]){
  double x1,  y1, z1;
  double x2,  y2, z2;
  double c_ra, c_dec;

  c_ra    =  centreRA*(pi/180.);
  c_dec   = centreDec*(pi/180.);
  
  // basis formed by: normal spherical co-ordinates subject to inversion through xy plane, then R1 and finally R2.
  for(j=0; j<Number; j++){
    // R1: rotation about z such that in the new basis, (x',y',z'), x' hat lies in x-y plane at an angle centreRA to x.
    x1  =     cos(c_ra)*xCoors[j] + sin(c_ra)*yCoors[j];
    y1  =    -sin(c_ra)*xCoors[j] + cos(c_ra)*yCoors[j];
    z1  =               zCoors[j];
    
    // R2: rotation about y such that in the new basis, (x'', y'', z''), z'' hat lies in (x', z') plane at an angle -CentreDec to x' hat.
    x2  = -sin(c_dec)*x1  - cos(c_dec)*z1;
    y2  =  y1;
    z2  =  cos(c_dec)*x1  - sin(c_dec)*z1;

    xCoors[j] = x2 + stefano_trans_x;  // Translate to fit in the box. P(k) unaffected.
    yCoors[j] = y2 + stefano_trans_y;
    zCoors[j] = z2 + stefano_trans_z;
  }

  return 0;
}
