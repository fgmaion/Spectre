#define NNCHAR 40
#define NNFIELD 200
#define NNVERTICES 100

double PEPS = 1e-12;

NodeP *init_poly(char *file, double *x0){
  FILE *fic = fopen(file, "r");
  
  double xmin[2], xmax[2];
  
  NodeP *polyTree = readPolygonFileTree(fic, xmin, xmax);
  
  fclose(fic);

  x0[0] = xmin[0] - 1.0;
  x0[1] = xmin[1] - 1.0;

  return polyTree;
}


int inPoly(NodeP *polyTree, double *x0, double ra, double dec){
  int id;
  
  double xin[2];

  xin[0] =  ra*180/M_PI;
  xin[1] = dec*180/M_PI;
  
  return insidePolygonTree(polyTree, x0, xin, &id);
}


int idPoly(NodeP *polyTree, double *x0, double ra, double dec){
  int id;
  
  double xin[2];

  xin[0] =  ra*180/M_PI;
  xin[1] = dec*180/M_PI;
  
  insidePolygonTree(polyTree, x0, xin, &id);
  
  return id;
}


double distAngSpher(const double RA1, double DEC1, double RA2, double DEC2){
  // Returns the angular distance between points with RA,DEC (in degree) 
  double sin2_ra  = 0.5*(1.0 - cos(RA1*PI/180.0)*cos(RA2*PI/180.0)  - sin(RA1*PI/180.0)*sin(RA2*PI/180.0));
  double sin2_dec = 0.5*(1.0 - cos(DEC1*PI/180.0)*cos(DEC2*PI/180.0)- sin(DEC1*PI/180.0)*sin(DEC2*PI/180.0));
  
  return 2.0*asin(sqrt(MAX(PEPS/100.0, sin2_dec + cos(DEC1*PI/180.0)*cos(DEC2*PI/180.0)*sin2_ra)))*180.0/PI;
}


void rotate(double x0, double y0, double x, double y, double *xrot, double *yrot, double angle, int spherical){
  // Takes positions (in degree if spherical = 1) and returns the rotated coordinates.
  if(spherical){
    double sign  = x - x0 < 0.0 ? -1.0 : 1.0;
    
    *yrot = -sign*(distAngSpher(x0, y0, x, y0))*sin(angle) + (y-y0)*cos(angle)+y0;
    
    *xrot =  sign*(distAngSpher(x0, y0, x, y0))*cos(angle) + (y-y0)*sin(angle)+x0;
    
    *xrot =  2.0*asin(sin((*xrot-x0)*PI/180.0/2.0)/cos(*yrot*PI/180.0))*180.0/PI+x0;
  }
  
  else{
    *xrot = +(x-x0)*cos(angle) - (y-y0)*sin(angle)+x0;
    
    *yrot =  (x-x0)*sin(angle) + (y-y0)*cos(angle)+y0;
  }
}


int insidePolygonTree(NodeP *polyTree, double x0[2], double x[2], int *poly_id){
  // Returns 1 if the point (x,y) is inside one of the polygons in polys. Returns 0 if the object is oustide of any polygon or outside the 
  // mask limits. See insidePolygon() for the algorithm explanations.
  
  int i, j, k, Ncross, result;
  
  double s, t, D; 
  
  if(polyTree->Npolys == 0){
    *poly_id = -1;
    
    return 0;
  }

  if(polyTree->type != LEAF){
    // pass to left child. 
    if(x[polyTree->SplitDim] < polyTree->SplitValue){ 
        result = insidePolygonTree(polyTree->Left,  x0, x, poly_id);
    }
    
    // pass to right child. 
    else{ 
        result = insidePolygonTree(polyTree->Right, x0, x, poly_id);
    }
  }
  
  else{
    // only evaluated for leaves. no need to store poly ids for all nodes, only leaves.  
    Polygon *polys = (Polygon *) polyTree->polysAll;
    
    for(k=0; k<polyTree->Npolys; k++){
      i = polyTree->poly_id[k];
      
      // x is the angular position of the galaxy to be counted. 
      if((polys[i].xmin[0] < x[0]) && (x[0] < polys[i].xmax[0]) && (polys[i].xmin[1] < x[1]) && (x[1] < polys[i].xmax[1])){
	    // The object is inside the rectangle formed from extreme points of the polygon.
	    Ncross = 0;
	
	    for(j=0; j<polys[i].N; j++){
	        if(j<polys[i].N-1){
	            D = (polys[i].x[j+1]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[j+1]-polys[i].y[j])*(x[0]-x0[0]);
	            s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
	            t = ((polys[i].x[j]-x[0])*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[j+1]-polys[i].x[j]))/D;
	        }
	        
	        else{
	            D = (polys[i].x[0]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[0]-polys[i].y[j])*(x[0]-x0[0]);
	            s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
	            t = ((polys[i].x[j]-x[0])*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[0]-polys[i].x[j]))/D;
	        } 
	  
	        if((0.0 < s) && (s < 1.0 + PEPS) && (0.0 < t) && (t < 1.0 + PEPS)) Ncross++;	
	    }
	
	    if(GSL_IS_ODD(Ncross)){
	        *poly_id = i;
	
	        // returns on first polygon.   
	        return 1;
	    }
      }
    }
    
    // point is not within any polygon. 
    *poly_id = -1;
    
    return 0;
  }
  
  return result;
}


int addgal_overlapPoly_counts(NodeP *polyTree, double *x0, double ra, double dec, double* weights){
  int id;
  
  double xin[2];

  xin[0] =  ra*180/M_PI;
  xin[1] = dec*180/M_PI;
  
  overlapping_polygonCounts(polyTree, x0, xin, &id, weights);
  
  return 0;
}


int overlapping_polygonCounts(NodeP *polyTree, double x0[2], double x[2], int *poly_id, double* weights){
  // given a galaxy, and a list of weights[Npolygon], add 1 to weight[i] if gal is in polygon i. 
  // gal can be in more than one polygon. 
  
  int i, j, k, Ncross, result;
  
  double s, t, D; 
  
  // no polygons in tree. 
  if(polyTree->Npolys == 0){
    *poly_id = -1;
    
    return 0;
  }


  if(polyTree->type != LEAF){
    // leaf if contains no polygons or area < min. area
    if(x[polyTree->SplitDim] < polyTree->SplitValue){ 
        result = overlapping_polygonCounts(polyTree->Left,  x0, x, poly_id, weights);
    }
    
    else{ 
        result = overlapping_polygonCounts(polyTree->Right, x0, x, poly_id, weights);
    }
  }
  
  
  else{
    Polygon *polys = (Polygon *) polyTree->polysAll;
    
    for(k=0; k<polyTree->Npolys; k++){
      // for each polygon associated to the node, (angular patch on the sky). 
      i = polyTree->poly_id[k];
      
      if((polys[i].xmin[0] < x[0]) && (x[0] < polys[i].xmax[0]) && (polys[i].xmin[1] < x[1]) && (x[1] < polys[i].xmax[1])){
	    // The object is inside the rectangle formed from extreme points of the polygon, i.e. first test that it's inside the polygon. 
	    Ncross = 0;
	
	    for(j=0; j<polys[i].N; j++){
	        if(j<polys[i].N-1){
	            D = (polys[i].x[j+1]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[j+1]-polys[i].y[j])*(x[0]-x0[0]);
	            s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
	            t = ((polys[i].x[j]-x[0])*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[j+1]-polys[i].x[j]))/D;
	        }
	        
	        else{
	            D = (polys[i].x[0]-polys[i].x[j])*(x[1]-x0[1])-(polys[i].y[0]-polys[i].y[j])*(x[0]-x0[0]);
	            s = ((x[0]-x0[0])*(polys[i].y[j]-x[1])-(x[1]-x0[1])*(polys[i].x[j]-x[0]))/D;
	            t = ((polys[i].x[j]-x[0])*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-x[1])*(polys[i].x[0]-polys[i].x[j]))/D;
	        } 
	  
	        if((0.0 < s) && (s < 1.0 + PEPS) && (0.0 < t) && (t < 1.0 + PEPS)) Ncross++;	
	    }
	
	    if(GSL_IS_ODD(Ncross)){
	        *poly_id = i;
	
	        // gal lies inside polygon. i is an absolute id for the polygon, ie relative to original catalogue, not list
	        // corresponding to the leaf (which would be k).
	        weights[i] += 1;
	    }
      }
    
        // repeat for all polygons associated to the leaf in which the galaxy lies.  
    }
    
    // point is not within any polygon. 
    *poly_id = -1;
    
    return 0;
  }
  
  return result;
}


int insideOnePolygon(NodeP *polyTree, int i, double *X, double x, double y){
  int j, Ncross;
  
  double s, t, D, x0, y0;
  
  Polygon *polys = (Polygon *) polyTree->polysAll;

  x0 = X[0];
  y0 = X[1];
  
  x *= 180/M_PI;
  y *= 180/M_PI;

  // i = polyTree->poly_id[i];

  if(polys[i].xmin[0] < x && x < polys[i].xmax[0] && polys[i].xmin[1] < y && y < polys[i].xmax[1]) {
    Ncross=0;
  
    for(j=0;j<polys[i].N;j++){
      if(j<polys[i].N-1){
	    D = (polys[i].x[j+1]-polys[i].x[j])*(y-y0)-(polys[i].y[j+1]-polys[i].y[j])*(x-x0);
	    s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
	    t = ((polys[i].x[j]-x)*(polys[i].y[j+1]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[j+1]-polys[i].x[j]))/D;
      }
      
      else{
	    D = (polys[i].x[0]-polys[i].x[j])*(y-y0)-(polys[i].y[0]-polys[i].y[j])*(x-x0);
	    s = ((x-x0)*(polys[i].y[j]-y)-(y-y0)*(polys[i].x[j]-x))/D;
	    t = ((polys[i].x[j]-x)*(polys[i].y[0]-polys[i].y[j])-(polys[i].y[j]-y)*(polys[i].x[0]-polys[i].x[j]))/D;
      } 
      
      if(0.0 < s && s < 1.0 + PEPS && 0.0 < t && t < 1.0 + PEPS) Ncross++;	
    }
    
    if(GSL_IS_ODD(Ncross)) return 1;
  }

  return 0;
}


Polygon *readPolygonFile(FILE *fileIn, int *Npolys, NodeP *polyTree){
  Polygon *result = malloc(sizeof(Polygon));
  return result;
}


NodeP *readPolygonFileTree(FILE *fileIn, double xmin[2], double xmax[2]){
  // Reads the file file_in and returns the polygons tree.
  char   line[NNFIELD*NNCHAR], item[NNFIELD*NNCHAR], *str_begin, *str_end;
  int    i, j;
  size_t N, NpolysAll;
  double minArea;
  int    SplitDim = 0, firstCall = 1;
  
  NpolysAll = 0;
  
  // Read the entire file and count the total number of polygons, NpolysAll.
  while(fgets(line, NNFIELD*NNCHAR, fileIn) != NULL) if(strstr(line, "polygon") != NULL) NpolysAll += 1;
  
  rewind(fileIn);
  
  // memory for NpolysAll polygons. 
  // Polygon* polysAll = (Polygon *) malloc(NpolysAll*sizeof(Polygon));
  polysAll = (Polygon *) realloc(polysAll, NpolysAll*sizeof(Polygon));
  
  for(i=0; i<NpolysAll; i++){
    polysAll[i].x    = (double *) malloc(NVERTICES*sizeof(double));
    polysAll[i].y    = (double *) malloc(NVERTICES*sizeof(double));
    
    polysAll[i].xmin = (double *) malloc(2*sizeof(double));
    polysAll[i].xmax = (double *) malloc(2*sizeof(double));
  }
  
  i = 0;
  
  // Read the file and fill the array with polygons.
  while(fgets(line, NNFIELD*NNCHAR, fileIn) != NULL){
    if(strstr(line,"polygon") != NULL){
      str_begin = strstr(line,"(") + sizeof(char);
      str_end = strstr(line, ")");
      
      // char * strcpy ( char * destination, const char * source );
      // Copies the C string pointed by source into the array pointed by destination, including the terminating null character (and stopping at that point).
      strcpy(str_end, "\n\0");
      
      getStrings(str_begin, item, ",", &N);
      
      // get all coordinates separated by comas.
      polysAll[i].N = N/2;
      
      if(N/2 > NVERTICES){
	    fprintf(stderr,"%s: %zd = too many points for polygon %d (%d maxi). Exiting...\n", MYNAME, N/2, i, NVERTICES);
	    exit(EXIT_FAILURE);
      }
      
      polysAll[i].id      = i;
      
      // initilisation. 
      polysAll[i].xmin[0] = atof(item);
      polysAll[i].xmax[0] = atof(item);
      
      // printf("\n%e \t %e", polysAll[i].xmin[0], polysAll[i].xmax[0]);
      
      polysAll[i].xmin[1] = atof(item+NNCHAR);
      polysAll[i].xmax[1] = atof(item+NNCHAR);
      
      // printf("\n%e \t %e", polysAll[i].xmin[1], polysAll[i].xmax[1]);
      
      for(j=0; j<N/2; j++){
        // read N/2 tuples of (x, y).
	    polysAll[i].x[j] = atof(item+NNCHAR*2*j);
	    polysAll[i].y[j] = atof(item+NNCHAR*(2*j+1));
	    
	    // vertices at minimum and maximum right ascension. 
	    polysAll[i].xmin[0] = MIN(polysAll[i].xmin[0], polysAll[i].x[j]);
	    polysAll[i].xmax[0] = MAX(polysAll[i].xmax[0], polysAll[i].x[j]);
	    
	    // vertices at minimum and maximum declination. 
	    polysAll[i].xmin[1] = MIN(polysAll[i].xmin[1], polysAll[i].y[j]);
	    polysAll[i].xmax[1] = MAX(polysAll[i].xmax[1], polysAll[i].y[j]);
      }
      
      // printf("\n%e \t %e", polysAll[i].xmin[0], polysAll[i].xmax[0]);
      // printf("\n%e \t %e", polysAll[i].xmin[1], polysAll[i].xmax[1]);
      
      // for each polygon. 
      i++;
    }
  }
  
  if(i==0){
    fprintf(stderr,"%s: 0 polygon found, check input file. Exiting...\n",MYNAME);
    
    exit(EXIT_FAILURE);
  } 
  
  else fprintf(stderr,"%d polygon(s) found\n", i);
  
  xmin[0] = polysAll[0].xmin[0];
  xmax[0] = polysAll[0].xmax[0];
  
  xmin[1] = polysAll[0].xmin[1];
  xmax[1] = polysAll[0].xmax[1];
  
  // assuming polygon is a rectangle on the sky. 
  minArea = (polysAll[0].xmax[0] - polysAll[0].xmin[0])*(polysAll[0].xmax[1] - polysAll[0].xmin[1]);
  
  for(i=1; i<NpolysAll; i++){
    // vertices at minimum and maximum right ascension, minimum/maximum of all polygons. 
    xmin[0]  = MIN(xmin[0], polysAll[i].xmin[0]);
    xmax[0]  = MAX(xmax[0], polysAll[i].xmax[0]);
    
    // vertices at minimum and maximum declination, minimum/maximum of all polygons. 
    xmin[1]  = MIN(xmin[1], polysAll[i].xmin[1]);
    xmax[1]  = MAX(xmax[1], polysAll[i].xmax[1]);
    
    minArea += (polysAll[i].xmax[0] - polysAll[i].xmin[0])*(polysAll[i].xmax[1] - polysAll[i].xmin[1]);
  }
  
  // polygon average area, assuming they are rectangular.
  minArea /= 1.0*(double) NpolysAll;
  
  printf("\n\nminimum area (rectangular) of polygons, %e", minArea);
  
  NodeP* result;
  
  result = createNodeP(polysAll, NpolysAll, minArea, SplitDim, xmin, xmax, firstCall); 
  
  return result;
}


NodeP *createNodeP(Polygon* polys, size_t Npolys, double minArea, int SplitDim, double xmin[2], double xmax[2], int firstCall){  
  // allocate memory for this node
  NodeP* result = (NodeP *) malloc(sizeof(NodeP));
  
  static size_t countNodes, NpolysAll;
  
  static void   *root, *polysAll;
  
  // Root & node id
  if(firstCall){
    root         = result;
    countNodes   = 0;
    
    // data for all polygons. 
    polysAll     = polys;
    NpolysAll    = Npolys;   // total number of polygons. 
  }
  
  result->root      = root;
  result->id        = countNodes;
  result->Npolys    = Npolys;
  result->NpolysAll = NpolysAll;
  
  countNodes++;
  
  // Copy address of the complete polygon sample.
  result->polysAll = polysAll;
  
  // save ids of polygons inside the node.  moved to leaf only save. 
  
  // area of this node.  patch on the sky. 
  double area = (xmax[0] - xmin[0])*(xmax[1] - xmin[1]);
  
  // leaf: either no polygon or cell smaller than minArea  
  if(result->Npolys == 0 || area < minArea){
    result->type     = LEAF;
    
    // no children. 
    result->Left     = NULL;
    result->Right    = NULL;
    
    result->poly_id  = (int *) malloc(Npolys*sizeof(int));
  
    // list of polygons associated to this leaf. 
    for(i=0; i<Npolys; i++)  result->poly_id[i] = polys[i].id;
  }
  
  
  else{    
    result->type       = NODE;
    
    // dimension along which the rectangle on the sky is to be sub-divided. 
    result->SplitDim   = SplitDim;
    
    // split midway. 
    result->SplitValue = (xmax[result->SplitDim] + xmin[result->SplitDim])/2.0;
    
    // each 'Node', a rectangle on the sky (resulting from recursively sub-diving an inital patch) has an associated list of polygons. 
    // those whose minimum point in a given ra or dec lies to the left of splitValue. 
    Polygon* polysChild = (Polygon *) malloc(Npolys*sizeof(Polygon));
    
    double xminChild[2], xmaxChild[2];
    
    for(i=0; i<Npolys; i++){
      // ra and dec for all polygon vertices. 
      polysChild[i].x    = (double *) malloc(NNVERTICES*sizeof(double));
      polysChild[i].y    = (double *) malloc(NNVERTICES*sizeof(double));
      
      // min and max ra and dec. 
      polysChild[i].xmin = (double *) malloc(2*sizeof(double));
      polysChild[i].xmax = (double *) malloc(2*sizeof(double));
    }
    
    for(i=0; i<2; i++){
      // set angular limits to same as parent. 
      xminChild[i] = xmin[i];
      xmaxChild[i] = xmax[i];
    }
    
    // new splitDim for children
    SplitDim++;
    
    if(SplitDim>1) SplitDim = 0;
    
    // subdivide. first for the left side, set new right limit.
    xmaxChild[result->SplitDim] = result->SplitValue;
    
    j=0;
    
    for(i=0; i<Npolys; i++){
      // minimum value of polygon along split dimension lies to left of split. add to list of polygons associated to child. 
      if(polys[i].xmin[result->SplitDim] < result->SplitValue){
        // copies polygon from polys, to polysChild.
	    cpyPolygon(&polysChild[j], &polys[i]);
	    
	    j++;
      }
    }
    
    // pass list of polygons, number of polygons in left child = j. 
    result->Left = createNodeP(polysChild, j, minArea, SplitDim, xminChild, xmaxChild, 0);
    
    // for right child... restore right limit
    xmaxChild[result->SplitDim] = xmax[result->SplitDim];
    
    // and set new left limit.
    xminChild[result->SplitDim] = result->SplitValue;
    
    j=0;
    
    for(i=0;i<Npolys;i++){
      // polys[i].xmin[result->splitDim] must be >= splitValues. so therefore must be polys[i].xmax.
      if(polys[i].xmax[result->SplitDim] > result->SplitValue){
	    cpyPolygon(&polysChild[j], &polys[i]);
	    j++;
      }
    } 
    
    result->Right = createNodeP(polysChild, j, minArea, SplitDim, xminChild, xmaxChild, 0);
    
    // done with polysChild. free.
    for(kk=0; kk<Npolys; kk++){
        free(polysChild[kk].x);
        free(polysChild[kk].y);
    
        free(polysChild[kk].xmin);
        free(polysChild[kk].xmax);
    }
  }
  
  result->Nnodes = countNodes;
  
  return result;
}


void free_NodeP(NodeP *node){
  if(node->type == LEAF){
    // free list of polygon ids for each node. 
    free(node->poly_id);
    free(node);
  }
  
  else{
    free_NodeP(node->Left);
    free_NodeP(node->Right);
    
    free(node);
  }
  
  return;
}


void cpyPolygon(Polygon *a, Polygon *b){
  // Copies b into a
  size_t i;
  
  a->N  = b->N;
  a->id = b->id;
  
  for(i=0; i<NVERTICES; i++){
    a->x[i] = b->x[i];
    a->y[i] = b->y[i];
  }
  
  for(i=0;i<2;i++){
    a->xmin[i] = b->xmin[i];
    a->xmax[i] = b->xmax[i];
  }
}
