/* -----------------------------------------------------------------------------*
 *  Creates a tree for the set of particles "particleSet" and returns a pointer *
 *  on the tree root.                                                           *
 * -------------------------------------------------------------------------    */

Node* buildTree(Particle *particleSet, int N, int Ndim){  
  // Computes the range of x, y and z co-ordinates for given set and assigns to xmin[3], xmax[3].
    
  // Initialise xmin, xmax with first particle info.
  for(j=0; j<Ndim; j++){
      xmin[j] = particleSet[0].x[j];
	  xmax[j] = particleSet[0].x[j];
  }
  
  for(i=1; i<N; i++){
    for(j=0; j<NDIM; j++){
        xmin[j] = MIN(particleSet[i].x[j], xmin[j]);
	    xmax[j] = MAX(particleSet[i].x[j], xmax[j]);
    }
  }
  
  // Decide which dimension to split along first, seems to be trying to choose the Dimension with the largest input range. fixed.
  tree_labelCount = firstSplitDim  = 0;
  
  // Beginning new tree, first node at level zero. 
  treelevel = 1;
  tree_N    = N;
  
  for(j=1; j<NDIM; j++) if((xmax[j]-xmin[j]) > (xmax[firstSplitDim] - xmin[firstSplitDim])) firstSplitDim = j;
  
  printf("\n\nPlanting tree.");
  
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


int printSorted(){
    for(j=0; j<10; j++)  printf("\n%d", point_gals[j].index);

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
  
      // printf("\t %d \t %d", result->NLeft, result->NRight);
  
      // printf("\n%d", result->label);
  } 
  
  else{
    // Contains the info on all particles grouped to be 'left' after the first split. Similarly for the right hand side. 

    result->Children    =      1;

    sortDim = SplitDim;

    // Arrange particles by their co-ordinate along dimenion splitDim.
    qsort(particleSet, N, sizeof(particleSet[0]), sortby_position_alongDim_splitDim);
    
    switchSplit(N, &NLeft, &NRight);
    
    result->NLeft       =  NLeft;
    result->NRight      = NRight;
  
    // printf("\t %d \t %d", result->NLeft, result->NRight);
  
    // Split value is the co-ordinate along dimension Split Dimension for which there are Nleft particles to the left, Nright particles to the right. Half value bewteen the rightmost left particle
    // and the leftmost right particle. 
    SplitValue          = 0.5*(particleSet[NLeft-1].x[SplitDim] + particleSet[NLeft].x[SplitDim]); 
  
    result->SplitValue  = SplitValue;
    
    Particle *lset      = &particleSet[0];
    Particle *rset      = lset + NLeft;
    
    // Limits, initially same as the parent but ..
    for(j=0; j<NDIM; j++){
        rmin[j] = result->xmin[j];
        rmax[j] = result->xmax[j];
        
        lmin[j] = result->xmin[j];
        lmax[j] = result->xmax[j];
    }
    
    // Limits.
    lmax[SplitDim]      = particleSet[NLeft-1].x[SplitDim];
    rmin[SplitDim]      =   particleSet[NLeft].x[SplitDim];
    
    if(SplitDim==firstSplitDim){
        treelevel += 1;
    }
    
    // Nice, recursive call.  Having split along one dimension, repeat along the next. 
    result->Left                           = createNode(lset, result->NLeft,  SplitDim + 1, lmin, lmax);
    
    result->Right                          = createNode(rset, result->NRight, SplitDim + 1, rmin, rmax);
  }
  
  return result;
}


int free_tree(Node *t){
  free(t->particle);
  
  if(t->Left !=NULL) free_tree(t->Left);
  if(t->Right!=NULL) free_tree(t->Right);
  
  free(t);
  
  return 0;
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
  Interim = 0.0;
  
  // Zero contribution to |displacement| if the nodes overlap in a given dimension.
  for(i=0; i<NDIM; i++){
      // 2  left of 1
      if(node1->xmin[i] > node2->xmax[i])       Interim += log10(node1->xmin[i] - node2->xmax[i]);

      // 2 right of 1
      else if(node2->xmin[i] > node1->xmax[i])  Interim += log10(node2->xmin[i] - node1->xmax[i]);
  
      else Interim = -99.;
  }  
  
  return Interim;
}


double log10_maximum_modDisplacementBetweenNodes(Node *node1, Node *node2){
  Interim = 0.0;
  
  // Zero contribution to |displacement| if the nodes overlap in a given dimension.
  for(i=0; i<NDIM; i++){
      // printf("\n %e", MAX(node1->xmin[i] - node2->xmax[i], node1->xmax[i] - node2->xmin[i]));
      
        Interim += pow(MAX(fabs(node1->xmin[i] - node2->xmax[i]), fabs(node1->xmax[i] - node2->xmin[i])), 2.);
  }
  
  return 0.5*log10(Interim);
}
