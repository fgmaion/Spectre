int cic_assign(int rand_ordelta, double x, double y, double z, double p){
  // Assigns an object at (x, y, z) into a 3D density grid with weight p using Cloud-in-Cell assignment. 
  // e.g xCellSize is the cell size in the x direction.  

  int              indx,  indy,  indz;
  int             iindx, iindy, iindz;
  
  double                   dx, dy, dz;
  double hx0, hy0, hz0, hx1, hy1, hz1;

  x = (x - AxisLimsArray[0][2])/xCellSize;
  y = (y - AxisLimsArray[0][1])/yCellSize;
  z = (z - AxisLimsArray[0][0])/zCellSize;

  // same as "box label" co-ordinates. 
  indx = (int) floor(x);
  indy = (int) floor(y);
  indz = (int) floor(z);
  
  // distance from centre of cell.
  dx = x - (double) indx - 0.5; // 0.5 for centre of cell labelled by (indx, indy, indz).
  dy = y - (double) indy - 0.5;
  dz = z - (double) indz - 0.5;

  hx0 = 1.0 - fabs(dx);
  hx1 =       fabs(dx);

  hy0 = 1.0 - fabs(dy);
  hy1 =       fabs(dy);

  hz0 = 1.0 - fabs(dz);
  hz1 =       fabs(dz);

  // periodic boundaries.
  if(indx>=n2) indx -= n2;
  if(indy>=n1) indy -= n1;
  if(indz>=n0) indz -= n0;
  
  // right of cell centre -> up a cell, otherwise down.
  if(dx>0) iindx=indx+1; else iindx=indx-1;
  if(dy>0) iindy=indy+1; else iindy=indy-1;
  if(dz>0) iindz=indz+1; else iindz=indz-1;
   
  // periodic boundaries. 
  if(iindx>=n2) iindx -= n2;
  if(iindy>=n1) iindy -= n1;
  if(iindz>=n0) iindz -= n0;
  
  if(iindx<0)   iindx +=n2;
  if(iindy<0)   iindy +=n1;
  if(iindz<0)   iindz +=n0;

  // Assignment of 8 corners of cube in which the particle sits,
  // Weights are ------
  //                    
  //             ------
  //
  // W = product_i of W_i(x_i - x_cell) = 1. - |x_i - x_cell|
  // Note: e.g. weight at bottom left   = 1. - |dx|, at bottom right = 1. - (1. - dx) 

  /*
  overdensity[iindx+n2*iindy+n1*n2*iindz][0] += hx1*hy1*hz1*p;
  overdensity[indx +n2*iindy+n1*n2*iindz][0] += hx0*hy1*hz1*p;
  overdensity[iindx+n2*indy +n1*n2*iindz][0] += hx1*hy0*hz1*p;
  overdensity[indx +n2*indy +n1*n2*iindz][0] += hx0*hy0*hz1*p;
  overdensity[iindx+n2*iindy+n1*n2*indz ][0] += hx1*hy1*hz0*p;
  overdensity[indx +n2*iindy+n1*n2*indz ][0] += hx0*hy1*hz0*p;
  overdensity[iindx+n2*indy +n1*n2*indz ][0] += hx1*hy0*hz0*p;
  overdensity[indx +n2*indy +n1*n2*indz ][0] += hx0*hy0*hz0*p;
  */
  
  overdensity[iindx + n2*iindy + n1*n2*iindz] += hx1*hy1*hz1*p;
  overdensity[indx  + n2*iindy + n1*n2*iindz] += hx0*hy1*hz1*p;
  overdensity[iindx + n2*indy  + n1*n2*iindz] += hx1*hy0*hz1*p;
  overdensity[indx  + n2*indy  + n1*n2*iindz] += hx0*hy0*hz1*p;
  overdensity[iindx + n2*iindy + n1*n2*indz ] += hx1*hy1*hz0*p;
  overdensity[indx  + n2*iindy + n1*n2*indz ] += hx0*hy1*hz0*p;
  overdensity[iindx + n2*indy  + n1*n2*indz ] += hx1*hy0*hz0*p;
  overdensity[indx  + n2*indy  + n1*n2*indz ] += hx0*hy0*hz0*p;
  
  return 0;
}
