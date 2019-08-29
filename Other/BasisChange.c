int VIPERSbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len){
    // Starting from Celestial co-ordinates, compute the co-ordinates in a new basis in which ^e1 lies along the VIPERS los.                                                                 
    double theta  = pi/2.0 - (pi/180.0)*centerDec;
    double phi    = centerRA*pi/180.0;

    double gamma1, gamma2, gamma3;

    // V = (x, y, z)^T = (y1, y2, y3)^T for {^x, ^y, ^z} and {^e1, ^e2, ^e3} respectively. 
    
    // where ^z corresponds to the Celestial pole and ^e1 corresponds to the VIPERS los.  
    
    // for t=theta, p=phi, S=sin(), C=cos(), e.g. sin(phi) = Sp

    // ^e1 = ( StCp,  StSp,  Ct)^T
    // ^e2 = (  -Sp,    Cp,   0)^T
    // ^e3 = (-CtCp, -CtSp,  St)^T

    for(j=0; j<len; j++){
        gamma1 =  xCoors[j]*sin(theta)*cos(phi)   + yCoors[j]*sin(theta)*sin(phi)   + zCoors[j]*cos(theta);
        gamma2 = -xCoors[j]*sin(phi)              + yCoors[j]*cos(phi)              + 0.;
        gamma3 = -xCoors[j]*cos(theta)*cos(phi)   - yCoors[j]*cos(theta)*sin(phi)   + zCoors[j]*sin(theta);

        xCoors[j] = gamma1;
        yCoors[j] = gamma2;
        zCoors[j] = gamma3;
    }

    /*
    printf("\nIn the VIPERS basis..");
    printf("\nx max:  %f \t x min:  %f", arrayMax(xCoor, Vipers_Num), arrayMin(xCoor, Vipers_Num));
    printf("\ny max:  %f \t y min:  %f", arrayMax(yCoor, Vipers_Num), arrayMin(yCoor, Vipers_Num));
    printf("\nz max:  %f \t z min:  %f", arrayMax(zCoor, Vipers_Num), arrayMin(zCoor, Vipers_Num));
    */
  
    return 0;
}


int Celestialbasis(double centerRA, double centerDec, double xCoors[], double yCoors[], double zCoors[], int len){
  // Starting from VIPERS los co-ordinates, compute the co-ordinates in a new basis in which ^e1 lies along the Celestial pole.

  double theta  = pi/2.0 - (pi/180.0)*centerDec;
  double phi    = centerRA*pi/180.0;

  double x, y, z;

  for(j=0; j<len; j++){
      x = xCoors[j]*sin(theta)*cos(phi) - yCoors[j]*sin(phi) - zCoors[j]*cos(theta)*cos(phi);
      y = xCoors[j]*sin(theta)*sin(phi) + yCoors[j]*cos(phi) - zCoors[j]*cos(theta)*sin(phi);
      z = xCoors[j]*cos(theta)          + 0.                 + zCoors[j]*sin(theta);

      xCoors[j] = x;
      yCoors[j] = y;
      zCoors[j] = z;
  }

  return 0;
}
