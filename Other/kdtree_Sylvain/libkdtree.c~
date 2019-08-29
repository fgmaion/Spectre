int growTrees(){
  assignLeafValues( point_gals,   xCoor,  yCoor,  zCoor,  Vipers_Num);

  assignLeafValues(point_rands,  rand_x, rand_y, rand_z, rand_number);

  galTree  = buildTree(point_gals,   Vipers_Num);

  randTree = buildTree(point_rands,  rand_number);

  galTree->label  = 0;

  randTree->label = 0;

  return 0;
}


int CountPairs_rMu(double **C, double **r, double **mu, Node* firstTree, Node* secndTree, int sameTree){
  findSuitableNodePairs_bruteforcePairCount(C, r, mu, firstTree, secndTree, sameTree);
  
  for(i=0; i<nlogbins; i++){
    for(j=0; j<nlinbins; j++){
         r[i][j] /= C[i][j];
        
        mu[i][j] /= C[i][j];
    }
  }
  
  sumPairs(C);

  return 0;
}


double computeNorm(int Nfirst, int Nsecnd){
    double norm;
    
    if(Nfirst == Nsecnd) norm = 0.5*Nfirst*(Nsecnd - 1.);
    else                 norm = 1.0*Nfirst*Nsecnd;
    
    return norm;
}


int assignLeafValues(Particle* cat, double xCoors[], double yCoors[], double zCoors[], int N){
    // Initialise tree with the co-ordinates of each particle. 
    
    for(i=0; i<N; i++){      	    
	    cat[i].x[0]        = xCoors[i];
	    cat[i].x[1]        = yCoors[i];
	    cat[i].x[2]        = zCoors[i];
	    
	    cat[i].index       = i;
	    cat[i].weight      = 1.0;
    }
    
    // for(i=0; i<10; i++) printf("\n %e", cat[i].x[0]);
    
    return 0;
}


int sumPairs(double** C){
    double totalpairs = 0.0;
  
    for(j=0; j<nlogbins; j++){
        for(i=0; i<nlinbins; i++){
            totalpairs += C[j][i];
        }
    }

    printf("\n\nTotal pairs: %e", totalpairs);

    return 0;
}


int findSuitableNodePairs_bruteforcePairCount(double **C, double **r, double **mu, Node *node1, Node *node2, int sameTree){
    if((node1->label%2 == 0) && (node2->label%500 == 0))  printf("\n%d \t %d", node1->label, node2->label);
    
    // Given two nodes, is their maximum displacement smaller than the smallest bin? is their minimum displacement larger than the largest bin -> don't bother finding their children for
    // counting pairs.    
    
    if( (log10_maximum_modDisplacementBetweenNodes(node1, node2) < 0.99*zerolog)){ 
         printf("\n %d \t %d \t %e \t %e \t %e \t %e", node1->label, node2->label, node1->xmin[2], node1->xmax[2], node2->xmin[2], node2->xmax[2]);   
    
        return 0; 
    }
    
    if( (log10_minimum_modDisplacementBetweenNodes(node1, node2) > 1.10*maxlog))  return 0;
    
    // Only count distinct pairs for two different children, take node1>node2 for instance. 
    if((sameTree == 1) && (node1->Children == 0) && (node2->Children == 0) && (node1->label<node2->label)) return 0;
    
    // Both node 1 and node 2 are leaves, their sub-divisons contain less than 200 particles. Brute force count pairs between node 1 and node 2.        
    if((node1->Children == 0) && (node2->Children == 0))  return bruteforceCountpairs_betweenNodes(C, r, mu, node1, node2, sameTree);
    
    else{
	    if((node1->Children == 0)){
	      // node 1 is a leaf, node 2 is not. Sub-divide node 2 and reevaluate for its children. 
	      findSuitableNodePairs_bruteforcePairCount(C, r, mu, node1,  node2->Left, sameTree);
	      findSuitableNodePairs_bruteforcePairCount(C, r, mu, node1,  node2->Right, sameTree);
	        
	      return 0;
	    } 
	
	    else{
	      // last scenario.. node 2 might still be a leaf, node 1 is not. Sub-divide node 1 and reevaluate for its children. 
	      findSuitableNodePairs_bruteforcePairCount(C, r, mu, node1->Left,  node2, sameTree);
	      findSuitableNodePairs_bruteforcePairCount(C, r, mu, node1->Right, node2, sameTree);
	    
	      return 0;
	    }
    }
    
    return 0;
}


int bruteforceCountpairs_betweenNodes(double **C, double **r, double **mu, Node *node1, Node *node2, int sameTree){    
    // Only called for children.      
    if((sameTree == 1) && (node1->label == node2->label)){
        // Same tree, same child. count distinct pairs.
        for(ii=0; ii<node1->N; ii++){
            for(jj=ii+1; jj<node2->N; jj++){
  		        // logarithmic binning in r.	    
    	        indi           = (int) floor(  (log10_particleSeparation(node1->particle[ii], node2->particle[jj]) - zerolog)/logbinsz);
	    
	            // linear binning in mu.
	            indj           = (int) floor(  (                pair_zmu(node1->particle[ii], node2->particle[jj]) - zerolin)/linbinsz);

                // Quicker to test for failures, less if evaluations. change this. 
                if((indi>=0) && (indj >=0) && (indi<nlogbins) && (indj<nlinbins)){
        	        C[indi][indj]  += 1.0;
    	            
        	        r[indi][indj]  += pow(10., log10_particleSeparation(node1->particle[ii], node2->particle[jj])); 
    	            
        	        mu[indi][indj] += pair_zmu(node1->particle[ii], node2->particle[jj]);
                }
            }
	    }
    }	    
	
	else{
	   // Different children. 
	   for(ii=0; ii<node1->N; ii++){
            for(jj=0; jj<node2->N; jj++){
          		// logarithmic binning in r.	    
            	indi           = (int) floor(  (log10_particleSeparation(node1->particle[ii], node2->particle[jj]) - zerolog)/logbinsz);
	    
	            // linear binning in mu.
	            indj           = (int) floor(  (                pair_zmu(node1->particle[ii], node2->particle[jj]) - zerolin)/linbinsz);
                
                // Quicker to test for failures, less if evaluations. change this. 
                if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0)){
                    C[indi][indj] += 1.0;
	                
	                r[indi][indj] += pow(10., log10_particleSeparation(node1->particle[ii], node2->particle[jj])); 
	                
	                mu[indi][indj]+= pair_zmu(node1->particle[ii], node2->particle[jj]);
	            }
		    }
	    }
	}
    
    return 0;
}


int bruteforce_nonodes(double **C, double **r, double **mu, Particle* cat, Particle* cat2, int N, int N2, int sameTree){    
    if(sameTree == 1){
        for(ii=0; ii<N; ii++){
            for(jj=ii+1; jj<N; jj++){
                // logarithmic binning in r.	    
                indi           = (int) floor(  (log10_particleSeparation(cat[ii], cat[jj]) - zerolog)/logbinsz);
	    
	            // linear binning in mu.
	            indj           = (int) floor(  (                pair_zmu(cat[ii], cat[jj]) - zerolin)/linbinsz);
   
    
                if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0))  C[indi][indj] += 1.0;
           }
        }
    }
    
    else{
        for(ii=0; ii<N; ii++){
            for(jj=0; jj<N2; jj++){
                // logarithmic binning in r.	    
                indi           = (int) floor(  (log10_particleSeparation(cat[ii], cat2[jj]) - zerolog)/logbinsz);
	    
	            // linear binning in mu.
		indj           = (int) floor(  (                pair_zmu(cat[ii], cat2[jj]) - zerolin)/linbinsz);
    
                if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0))  C[indi][indj] += 1.0;
           }
        }
    }
    
    sumPairs(C);

    return 0;
}


double log10_particleSeparation(Particle a, Particle b){
  double r2 = pow(10., zerolog);

  for(i=0; i<NDIM; i++){ 
      // printf("\n %e \t %e", point_gals[0].x[i], point_rands[0].x[i]);
    
    r2 += pow(a.x[i] - b.x[i], 2.);
  }
  
  // log10(r^2) -> 2.0 log_10(r)
  return 0.5*log10(r2);
}


double pair_zmu(Particle a, Particle b){    
  // assuming z is the polar axis, return |cos(theta)| of pair separation.
  return fabs(b.x[2] - a.x[2])/pow(10., log10_particleSeparation(a, b));
}
