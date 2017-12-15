int ngp_assign(double x, double y, double z, double p){
  // Assigns an object at (x, y, z) into a 3D density grid with weight p using Nearest-Grid-Point assignment. 
  // e.g xCellSize is the cell size in the x direction.  

  int indx, indy, indz;
  
  double nuCellSize, newlim;
  
  switch(fold){
    case 0:
      x   /= xCellSize;
      y   /= yCellSize;
      z   /= zCellSize;

      // same as "box label" co-ordinates.
      indx = (int) trunc(x); // as opposed to floor, due to behaviour for -ve x. https://www.gnu.org/software/libc/manual/html_node/Rounding-Functions.html
      indy = (int) trunc(y);
      indz = (int) trunc(z);

      // printf("\n%.6lf \t %.6lf \t %d", x, trunc(x), (int) trunc(x));
      
      // overdensity[indx + n2*indy + n1*n2*indz] += p;
      overdensity[indx + n2*indy + n1*n2*indz][0] += p;
      
      break;

      
    case 1:
      newlim     = AxisLimsArray[1][2]/FOLDFACTOR;
      nuCellSize =           xCellSize/FOLDFACTOR;
      
      x = fmod(x, newlim)/nuCellSize;
      y = fmod(y, newlim)/nuCellSize;
      z = fmod(z, newlim)/nuCellSize;

      // same as "box label" co-ordinates.
      indx = (int) trunc(x);
      indy = (int) trunc(y);
      indz = (int) trunc(z);

      // overdensity[indx + n2*indy + n1*n2*indz] += p;
      overdensity[indx + n2*indy + n1*n2*indz][0] += p;
      
      break;

  case 2:
    newlim     = AxisLimsArray[1][2]/(FOLDFACTOR*FOLDFACTOR);
    nuCellSize =           xCellSize/(FOLDFACTOR*FOLDFACTOR);

    x = fmod(x, newlim)/nuCellSize;
    y = fmod(y, newlim)/nuCellSize;
    z = fmod(z, newlim)/nuCellSize;

    // same as "box label" co-ordinates.
    indx = (int) trunc(x);
    indy = (int) trunc(y);
    indz = (int) trunc(z);

    // overdensity[indx + n2*indy + n1*n2*indz] += p;
    overdensity[indx + n2*indy + n1*n2*indz][0] += p;
    
    break;
  }
    
  return 0;
}
