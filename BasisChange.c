int VIPERSbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len){
  // Rotate the z Cartesian axis to line along the LOS.                                                                 

  float theta  = pi/2.0 - (pi/180.0)*centerDec;
  float phi    = centerRA*pi/180.0;

  float gamma1, gamma2, gamma3;

  for(j=0; j<len; j++){
    gamma1 =  xCoors[j]*sin(theta)*cos(phi)   + yCoors[j]*sin(theta)*sin(phi)   + zCoors[j]*cos(theta);
    gamma2 = -xCoors[j]*sin(phi)              + yCoors[j]*cos(phi)              + 0.;
    gamma3 = -xCoors[j]*cos(theta)*cos(phi)   - yCoors[j]*cos(theta)*sin(phi)   + zCoors[j]*sin(theta);

    xCoors[j] = gamma1;
    yCoors[j] = gamma2;
    zCoors[j] = gamma3;
  }

  return 0;
}


int Celestialbasis(float centerRA, float centerDec, float xCoors[], float yCoors[], float zCoors[], int len){
  // Rotate the z Cartesian axis to line along the LOS.                                                                                              
  float theta  = pi/2.0 - (pi/180.0)*centerDec;
  float phi    = centerRA*pi/180.0;

  float x, y, z;

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