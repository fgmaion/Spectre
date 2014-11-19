#define NDIM 3

/* -----------------------------------------------------------------------------*
 *  Creates a tree for the set of particles "particleSet" and returns a pointer *
 *  on the tree root.                                                           *
 * -------------------------------------------------------------------------    */

double minimum_twodoubles(double first, double second){
    if(first<=second){ 
        return first;
    }
    
    else{ 
        return second;
    }
}


double maximum_twodoubles(double first, double second){
    if(first>=second){ 
        return first;
    }
    
    else{ 
        return second;
    }
}


Node* buildTree(Particle *particleSet, int N){  
  // Computes the range of x, y and z co-ordinates for given set and assigns to xmin[3], xmax[3].
    
  double xmin[NDIM], xmax[NDIM];
    
  // Initialise xmin, xmax with first particle info.
  for(j=0; j<NDIM; j++){
      xmin[j] = particleSet[0].x[j];
	  xmax[j] = particleSet[0].x[j];
  }
  
  // printf("\n\n");
  
  for(i=1; i<N; i++){
    for(j=0; j<NDIM; j++){
        // printf("%e \t", particleSet[i].x[j]);        
        
        xmin[j] = minimum_twodoubles(particleSet[i].x[j], xmin[j]);    // MIN(particleSet[i].x[j], xmin[j]);
	    xmax[j] = maximum_twodoubles(particleSet[i].x[j], xmax[j]);    // MAX(particleSet[i].x[j], xmax[j]);
    }
    
    // printf("\n");
  }
  
  // bit of leeway.
  for(j=0; j<NDIM; j++) xmin[j] *= 0.99;
  for(j=0; j<NDIM; j++) xmax[j] *= 1.01;
  
  // Decide which dimension to split along first, seems to be trying to choose the Dimension with the largest input range. fixed.
  tree_labelCount = firstSplitDim  = 0;
  
  // Beginning new tree, first node at level zero. 
  treelevel = 1;
  tree_N    = N;
  
  for(j=1; j<NDIM; j++){ 
    if((xmax[j]-xmin[j]) > (xmax[firstSplitDim] - xmin[firstSplitDim])){ 
        firstSplitDim = j;
    }
  }
  
  printf("\n\nPlanting tree.");
  
  printf("\n\nSplitting tree along dimension: %d\n", firstSplitDim);
  
  printf("\n%e \t %e \t %e", xmin[0], xmin[1], xmin[2]);
  printf("\n%e \t %e \t %e", xmax[0], xmax[1], xmax[2]);
  
  printf("\n\n");
  
  // printf("\n\nTree has %d levels, first split occured along %d", treelevel, firstSplitDim);
  
  return createNode(particleSet, N, firstSplitDim, xmin, xmax);
}


int switchSplit(int N, int* NLeft, int* NRight){
    switch(N%2){
      case 0:
      // Even number of particles, split evenly left and right
        *NLeft  = N/2;
        *NRight = N/2;
        
        return 0;
        
      case 1:
      // Odd number of particles, split N-1 particles evenly, add extra particle to left hand side. 
        *NLeft  = (N+1)/2;
        *NRight = (N-1)/2;
    
        return 0;
    }
}


int printSorted(Particle *particleSet, int N, int NLeft, int NRight, double SplitValue){
    printf("\n\nN: %d, N left: %d, N right: %d, split value: %e", N, NLeft, NRight, SplitValue);

    for(j=0; j<NLeft; j++)  printf("\n%d \t %d \t %e", j, particleSet[j].index, particleSet[j].x[sortDim]);

    printf("\n\n");
    
    for(j=NLeft; j<N; j++)  printf("\n%d \t %d \t %e", j, particleSet[j].index, particleSet[j].x[sortDim]);

    printf("\n\n");

    return 0;
}

/* ------------------------------------------------------------------------- *
 *    Computes SplitValue and returns a pointer on the node, recursively.    *
 *    If N<Nmin returns a node with                                        *
 *    no child => node.Left = NULL and node.Right = NULL <=> leaf.           *
 *    To build a tree, just give a point, SplitDim, xmin and xmax.  *
 * ------------------------------------------------------------------------- */

Node* createNode(Particle *particleSet, int N, int SplitDim, double xmin[NDIM], double xmax[NDIM]){  
  // for(j=0; j<NDIM; j++)  printf("\n%d \t %d \t %e \t %e", tree_labelCount, N, xmin[j], xmax[j]);
  
  // Allocates memory for THIS node
  Node *result        = (Node *) malloc(sizeof(Node));
  
  // printSorted();
  
  // Writes N, SplitDim, and SplitValue
  result->N           = N;               // Number of particles in Node.
  
  result->treelevel   = treelevel;
  
  result->label       = tree_labelCount;
  
  tree_labelCount    += 1;
  
  result->particle    = particleSet;
  
  // Split along x,y,z. now go back to x for finer nodes. 
  if(SplitDim>(NDIM-1))  SplitDim   = 0;
  
  result->SplitDim    = SplitDim;        // Next split dimension.
  
  // Limits
  for(j=0; j<NDIM; j++){
    result->xmin[j]   = xmin[j];
    result->xmax[j]   = xmax[j];
  }
    
  // Less than Nmin particles, don't subdivide.
  if(N<=Nmin){
      result->Left        = NULL;
      result->Right       = NULL;
      
      result->NLeft       =    0;
      result->NRight      =    0;
      
      result->Children    =    0;
      
      result->SplitValue  = pow(10., 99.);  // No further split.
  } 
  
  else{
    // Children!
    // Contains the info on all particles grouped to be 'left' after the first split. Similarly for the right hand side. 

    result->Children    =      1;

    sortDim = SplitDim;

    // Arrange particles by their co-ordinate along dimenion sortDim.
    qsort(particleSet, N, sizeof(particleSet[0]), sortby_position_alongDim_splitDim);
    
    switchSplit(N, &NLeft, &NRight);
    
    result->NLeft       =  NLeft;
    result->NRight      = NRight;
  
    // Split value is the co-ordinate along dimension Split Dimension for which there are Nleft particles to the left, Nright particles to the right. Half value bewteen the rightmost left particle
    // and the leftmost right particle. 
    SplitValue          = 0.5*(particleSet[NLeft-1].x[SplitDim] + particleSet[NLeft].x[SplitDim]); 
  
    // printSorted(particleSet, N, NLeft, NRight, SplitValue);
  
    result->SplitValue  = SplitValue;
    
    double rmax[NDIM], rmin[NDIM], lmax[NDIM], lmin[NDIM];
    
    Particle *lset      = &particleSet[0];
    Particle *rset      = lset + NLeft;
    
    // Limits, initially same as the parent but ..
    for(j=0; j<NDIM; j++){
        rmin[j] = xmin[j];
        rmax[j] = xmax[j];
        
        lmin[j] = xmin[j];
        lmax[j] = xmax[j];
    }
    
    // Limits.
    lmax[SplitDim]      = particleSet[NLeft-1].x[SplitDim];
    rmin[SplitDim]      =   particleSet[NLeft].x[SplitDim];
    
    // printf("\n%e \t %e \t %e \t %e \t %e \t %e", lmin[0], lmin[1], lmin[2], lmax[0], lmax[1], lmax[2]);        
    // printf("\n%e \t %e \t %e \t %e \t %e \t %e", rmin[0], rmin[1], rmin[2], rmax[0], rmax[1], rmax[2]);

    // printf("\n\nBoundary particles: %e \t %e", particleSet[NLeft-1].x[SplitDim], particleSet[NLeft].x[SplitDim]);
    
    // printf("\n\n\n\n");
    
    if(SplitDim==firstSplitDim){
        treelevel += 1;
    }
    
    // Nice, recursive call.  Having split along one dimension, repeat along the next. 
    result->Left                           = createNode(lset, result->NLeft,  SplitDim + 1, lmin, lmax);
    
    result->Right                          = createNode(rset, result->NRight, SplitDim + 1, rmin, rmax);
  }
  
  return result;
}


Node* Create_toyChildNode(double xlo, double ylo, double zlo, double xhi, double yhi, double zhi){
    Node *toyNode = (Node *) malloc(sizeof(Node));

    toyNode->N    = 50;               // Number of particles in Node.
  
    toyNode->particle    = &point_gals[50];
    
    // Limits
    toyNode->xmin[0]   = xlo;
    toyNode->xmin[1]   = ylo;
    toyNode->xmin[2]   = zlo;
    
    toyNode->xmax[0]   = xhi;
    toyNode->xmax[1]   = yhi;
    toyNode->xmax[2]   = zhi;
    
    toyNode->Left          = NULL;
    toyNode->Right         = NULL;
      
    toyNode->NLeft         =    0;
    toyNode->NRight        =    0;
      
    toyNode->Children      =    0;
      
    toyNode->SplitValue    = pow(10., 99.);  // No further split. 
  
    return toyNode;
}


int sortby_position_alongDim_splitDim(const void *a, const void *b){
    if((*(Particle *)a).x[sortDim] >= (*(Particle *)b).x[sortDim]){  
        return  1;
    }
    
    else{
        return -1;
    }
}


double log10_minimum_modDisplacementBetweenNodes(Node *node1, Node *node2){  
  double dr2 = 0.0;
  
  // Sum of (minimum) projected distances^2
  // contribution to |displacement| is zero for a given dimension if the nodes overlap.
  for(i=0; i<NDIM; i++){
      // 2  left of 1
      if(node1->xmin[i] > node2->xmax[i]){       
          dr2   += pow(node1->xmin[i] - node2->xmax[i], 2.);
      }  
        
      // 2 right of 1
      else if(node2->xmin[i] > node1->xmax[i]){  
          dr2   += pow(node2->xmin[i] - node1->xmax[i], 2.);
      }
  
      // overlap
      else{ 
          dr2   += 0.;
      }
  }  
  
  // log_10(0.) will return inf. 
  return 0.5*log10(maximum_twodoubles(dr2, pow(10., -99.)));
}


double log10_maximum_modDisplacementBetweenNodes(Node *node1, Node *node2){
  double dr2 = 0.0;
  
  // This is an approximation to the maximum displacement, as that given by taking the maximum of dx, dy, dz individually.. and is therefore an overestimate. 
  for(i=0; i<NDIM; i++)  dr2 += pow(maximum_twodoubles(fabs(node1->xmin[i] - node2->xmax[i]), fabs(node1->xmax[i] - node2->xmin[i])), 2.);
  
  return 0.5*log10(dr2);
}


int free_tree(Node *t){
  free(t->particle);
  
  if(t->Left !=NULL) free_tree(t->Left);
  if(t->Right!=NULL) free_tree(t->Right);
  
  free(t);
  
  return 0;
}
