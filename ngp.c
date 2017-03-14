#define FOLDFACTOR 2.0

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
      indx = (int) floor(x);
      indy = (int) floor(y);
      indz = (int) floor(z);
     
      overdensity[indx + n2*indy + n1*n2*indz] += p;
      
    case 1:
      newlim     = AxisLimsArray[1][2]/FOLDFACTOR;
      nuCellSize =           xCellSize*FOLDFACTOR;
      
      x = fmod(x, newlim)/nuCellSize;
      y = fmod(y, newlim)/nuCellSize;
      z = fmod(z, newlim)/nuCellSize;

      // same as "box label" co-ordinates.
      indx = (int) floor(x);
      indy = (int) floor(y);
      indz = (int) floor(z);

      overdensity[indx + n2*indy + n1*n2*indz] += p;
  }
    
  return 0;
}
