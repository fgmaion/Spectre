int grow_randTree(){
  assignLeafValues(point_rands,  rand_x, rand_y, rand_z, rand_weight, rand_number);
    
  randTree = buildTree(point_rands,  rand_number);

  randTree->label = 0;

  printf("\n\nRand tree is fully grown! (%d nodes).", tree_labelCount);
  
  return 0;
}

int grow_galTree(){
  assignLeafValues(point_gals,   xCoor,  yCoor,  zCoor, zCoor, Vipers_Num);
  
  galTree  = buildTree(point_gals,   Vipers_Num);

  galTree->label  = 0;

  printf("\n\nRand tree is fully grown! (%d nodes).", tree_labelCount);
  
  return 0;
}

double computeNorm(int Nfirst, int Nsecnd){
    double norm;
    
    if(Nfirst == Nsecnd) norm = 0.5*Nfirst*(Nsecnd - 1.);
    else                 norm = 1.0*Nfirst*Nsecnd;
    
    return norm;
}

int assignLeafValues(Particle* cat, double xCoors[], double yCoors[], double zCoors[], double weight[], int N){
  // Initialise tree with the co-ordinates of each particle. 
  // #pragma omp parallel for private(i) if (thread == 1)
  for(i=0; i<N; i++){      	    
    cat[i].x[0]        = xCoors[i];
    cat[i].x[1]        = yCoors[i];
    cat[i].x[2]        = zCoors[i];
	
    cat[i].index       =         i;
    cat[i].weight      = weight[i];
  }
    
  return 0;
}

int sumPairs(double* C){
  double totalpairs = 0.0;
  
  for(i=0; i<nlogbins; i++){
    for(j=0; j<nlinbins; j++){
      Index       = j + nlinbins*i;

      totalpairs += C[Index];
    }
  }

  printf("\n\nTotal pairs: %e", totalpairs);

  return 0;
}

int print_tree(Node *node1){
  if(node1->Children != 0){
    print_tree(node1->Left);
        
    print_tree(node1->Right);
  }
  
  if(node1->Children == 0){
    Sumof_Childrenparticles += node1->N;
    
    printf("\n%d \t %d", node1->label, node1->N);
  }
  
  return 0;
}

int nodes_printParticles_print(){
  sprintf(filepath, "%s/Data/stacpolly/Node_particle_positions.dat", root_dir);
    
  output = fopen(filepath, "w");
    
  for(j=0; j<Vipers_Num; j++)  fprintf(output, "%e \t %e \t %e \n", xCoor[j], yCoor[j], zCoor[j]);

  fclose(output);

  return 0;
}

int prep_nodeLimits_print(){
  sprintf(filepath, "%s/Data/stacpolly/kdtree_nodeLimits_children.dat", root_dir);
    
  output = fopen(filepath, "w");
    
  print_nodeLimits(galTree);

  fclose(output);

  return 0;
}

int print_nodeLimits(Node *node){
  if(node->Children != 0){
    print_nodeLimits(node->Left);
        
    print_nodeLimits(node->Right);
  }
  
  if(node->Children == 0){
    fprintf(output, "%d \t %e \t %e \t %e \t %e \t %e \t %e \n", node->N, node->xmin[0], node->xmin[1], node->xmin[2], node->xmax[0], node->xmax[1], node->xmax[2]);
  }

  return 0;
}

int findSuitableNodePairs_bruteforcePairCount(double* n, double* wn, double* r, double* mu, Node* node1, Node* node2, int sameTree){
  if((node1->label >= nodeone_progresscount) && (node2->label >= nodetwo_progresscount)){      
    printf("\n%d \t %d", node1->label, node2->label);
      
    postprocesspairs(n, wn, r, mu, node1, node2);
      
    nodeone_progresscount = node1->label + (long) floor(tree_labelCount/1000.);
    nodetwo_progresscount = node2->label + (long) floor(tree_labelCount/1000.);
  }
    
  // Given two nodes, is their maximum displacement smaller than the smallest bin? is their minimum displacement larger than the largest bin
  // Don't bother counting pairs. Otherwise, brute force count between the nodes. 
  if((log10_minimum_modDisplacementBetweenNodes(node1, node2) > 1.01*maxlog))   return 0;
    
  // Only count distinct pairs, given two different children, take node1 >= node2 for instance. 
  if((sameTree == 1) && (node1->Children == 0) && (node2->Children == 0) && (node1->label < node2->label)) return 0;

  // if((node1->label < nodeone_savedcount) && (node2->label < nodetwo_savedcount))  return 0; // fast forward to saved counts.
  
  // Both node 1 and node 2 are leaves, their sub-divisons contain less than Nmin particles. Brute force count pairs between node 1 and node 2.        
  if((node1->Children == 0) && (node2->Children == 0))  return bruteforceCountpairs_betweenChildren(n, wn, r, mu, node1, node2, sameTree);
    
  else{
    if((node1->Children == 0)){
      // node 1 is a leaf, node 2 is not. Sub-divide node 2 and reevaluate for its children. 
      findSuitableNodePairs_bruteforcePairCount(n, wn, r, mu, node1,  node2->Left,  sameTree);
      findSuitableNodePairs_bruteforcePairCount(n, wn, r, mu, node1,  node2->Right, sameTree);
	        
      return 0;
    } 
	
    else{
      // last scenario.. node 2 might still be a leaf, node 1 is not. Sub-divide node 1 and reevaluate for its children. 
      findSuitableNodePairs_bruteforcePairCount(n, wn, r, mu, node1->Left,  node2, sameTree);
      findSuitableNodePairs_bruteforcePairCount(n, wn, r, mu, node1->Right, node2, sameTree);
	    
      return 0;
    }
  }
    
  return 0;
}

int bruteforceCountpairs_betweenChildren(double* n, double* wn, double* r, double* mu, Node *node1, Node *node2, int sameTree){   
  double  log10_r, weight;
  double  pair_mu, pair_mu2;
  
  // Only called for children.      
  if((sameTree == 1) && (node1->label == node2->label)){  // Same tree, same child. count distinct pairs.
   for(ii=0; ii<node1->N; ii++){
      for(jj=ii+1; jj<node2->N; jj++){
        log10_r                = log10_particleSeparation(node1->particle[ii], node2->particle[jj]); 
                    
        if((log10_r > zerolog) && (log10_r < maxlog)){
          pair_mu              = pair_zmu(node1->particle[ii], node2->particle[jj], log10_r);
          
          indi                 = (int) floor(  (log10_r - zerolog)/logbinsz);  // logarithmic binning in r.	      
          indj                 = (int) floor(  (pair_mu - zerolin)/linbinsz);  // linear binning in mu.
          
          if((indi < nlogbins) && (indj < nlinbins) && (indi >= 0) && (indj >= 0)){                      
            weight             = node1->particle[ii].weight*node2->particle[jj].weight;

            Index              =     indj + nlinbins*indi;
            
             n[Index]         +=                      1.0;
            wn[Index]         +=                   weight; 
             r[Index]         +=           log10_r*weight;                     // l=2, (2l+1) L_2 -> C2[indi][indj] +=(3.*pair_mu2 - 1.);
            mu[Index]         +=           pair_mu*weight;                     // l=4, (2l+1) L_4 -> C4[indi][indj] += (35.*pair_mu2*pair_mu2 - 30.*pair_mu2 + 3.);
          }
        }
                
        else{
          printf("\nUnderflow, reduce lower zerolog (%e) to: %e", zerolog, log10_r);
        }
      }
    }
  }	    
  
  else{
    #pragma omp parallel for reduction(+: n[:NBINS], wn[:NBINS], r[:NBINS], mu[:NBINS]) private(ii, jj, log10_r, pair_mu, indi, indj, weight, Index) if(thread==1) 
    for(ii=0; ii<node1->N; ii++){      // Either different tree, or different children with node1 > node2. all pairs are distinct at this point.
      for(jj=0; jj<node2->N; jj++){
        log10_r               = log10_particleSeparation(node1->particle[ii], node2->particle[jj]); 
                        
        if((log10_r > zerolog) && (log10_r < maxlog)){
          pair_mu             = pair_zmu(node1->particle[ii],  node2->particle[jj], log10_r);
          
          indi                = (int) floor(  (log10_r - zerolog)/logbinsz);  // logarithmic binning in r.	    
          indj                = (int) floor(  (pair_mu - zerolin)/linbinsz);  // linear binning in mu.
          
          if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0)){                      
            weight             = node1->particle[ii].weight*node2->particle[jj].weight;

            Index              =     indj + nlinbins*indi;

             n[Index]         +=                      1.0;
            wn[Index]         +=                   weight;
             r[Index]         +=           log10_r*weight;                     // l=2, (2l+1) L_2 -> C2[indi][indj] +=(3.*pair_mu2 - 1.);
            mu[Index]         +=           pair_mu*weight;                    // l=4, (2l+1) L_4 -> C4[indi][indj] += (35.*pair_mu2*pair_mu2 - 30.*pair_mu2 + 3.);    
          }
        }
        
        else{
          printf("\nUnderflow, reduce lower zerolog (%e) to: %e", zerolog, log10_r);
        }
      }
    }
  }

  return 0;
}


int bruteforce_nonodes(double* C0, double* C2, double* C4, double* C6, double* C8, double* C10, double* r, double* mu, Particle* cat, Particle* cat2, int N, int N2, int sameTree){
  double  log10_r, weight;
  double  pair_mu, pair_mu2, pair_mu4;
  
  if(sameTree == 1){
    // #pragma omp parallel for reduction(+: C0[:NBINS],C2[:NBINS],C4[:NBINS],C6[:NBINS],C8[:NBINS], Distinct_pairCount) private(ii, jj, log10_r, indi, indj, pair_mu, pair_mu2, pair_mu4, weight) if(thread==1)
    for(ii=0; ii<N; ii++){  // Counting distinct pairs.      
      for(jj=ii+1; jj<N; jj++){
        Distinct_pairCount += 1;
        
        log10_r             = log10_particleSeparation(cat[ii], cat[jj]); 
        
        pair_mu             = pair_zmu(cat[ii], cat[jj], log10_r);
                            
        if(log10_r>zerolog){  
          indi           = (int) floor(  (log10_r - zerolog)/logbinsz);  // logarithmic binning in r.
          indj           = (int) floor(  (pair_mu - zerolin)/linbinsz);  // linear binning in mu.

          if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0)){
            pair_mu2          = pair_mu*pair_mu;
            pair_mu4          = pair_mu2*pair_mu2;

            weight            = cat[ii].weight*cat[jj].weight;

            Index             =     indj + nlinbins*indi;

            C0[Index]        +=                   weight;
            C2[Index]        +=          pair_mu2*weight;              // l=2, (2l+1) L_2 -> C2[indi][indj] +=(3.*pair_mu2 - 1.);
            C4[Index]        +=          pair_mu4*weight;              // l=4, (2l+1) L_4 -> C4[indi][indj] += (35.*pair_mu2*pair_mu2 - 30.*pair_mu2 + 3.);
            C6[Index]        += pair_mu4*pair_mu2*weight;
            C8[Index]        += pair_mu4*pair_mu4*weight;
          }
        }
            
        else{
          printf("\nUnderflow, reduce lower zerolog (%e) to: %e", zerolog, log10_r);
        }
      }
    }
  }
    
  else{
    // #pragma omp parallel for reduction(+: C0[:NBINS],C2[:NBINS],C4[:NBINS],C6[:NBINS],C8[:NBINS], Distinct_pairCount) private(ii, jj, log10_r, indi, indj, pair_mu, pair_mu2, pair_mu4, weight) if(thread==1)
    for(ii=0; ii<N; ii++){      
      for(jj=0; jj<N2; jj++){
        Distinct_pairCount   += 1;
                
        log10_r               = log10_particleSeparation(cat[ii], cat2[jj]); 
            
        pair_mu               = pair_zmu(cat[ii], cat2[jj], log10_r);
                        
        if(log10_r>zerolog){
          indi              = (int) floor(  (log10_r - zerolog)/logbinsz);  // logarithmic binning in r.
          indj              = (int) floor(  (pair_mu - zerolin)/linbinsz);  // linear binning in mu.
                    
          if((indi<nlogbins) && (indj<nlinbins) && (indi>=0) && (indj >=0)){
            pair_mu2          = pair_mu*pair_mu;
            pair_mu4          = pair_mu2*pair_mu2;

            weight            = cat[ii].weight*cat[jj].weight;
            
            Index             =     indj + nlinbins*indi;

            C0[Index]        +=                   weight;
            C2[Index]        +=          pair_mu2*weight;              // l=2, (2l+1) L_2 -> C2[indi][indj] +=(3.*pair_mu2 - 1.);
            C4[Index]        +=          pair_mu4*weight;              // l=4, (2l+1) L_4 -> C4[indi][indj] += (35.*pair_mu2*pair_mu2 - 30.*pair_mu2 + 3.);
            C6[Index]        += pair_mu4*pair_mu2*weight;
            C8[Index]        += pair_mu4*pair_mu4*weight;
          }
        }
        
        else{
          printf("\nUnderflow, reduce lower zerolog (%e) to: %e", zerolog, log10_r);
        }
      }
    }
  }
    
  sumPairs(C0);

  printf("\n\nDistinct pair count: %d (%d)", Distinct_pairCount, rand_number*(rand_number - 1)/2);
    
  return 0;
}

double log10_particleSeparation(Particle a, Particle b){
  double r2 = 0.0;

  for(int i=0; i<NDIM; i++)  r2 += pow(a.x[i] - b.x[i], 2.);
  
  // printf("\n\n%e \t %e \t %e \t %e \t %e \t %e \t %e", a.x[0], a.x[1], a.x[2], b.x[0], b.x[1], b.x[2], sqrt(r2));
  
  // log10(r^2) -> 2.0 log_10(r). deal with log_10(0.)
  return 0.5*log10(maximum_twodoubles(r2, pow(10., -99.)));
}

double pair_zmu(Particle a, Particle b, double log10r){    
  // assuming z is the polar axis, return |cos(theta)| of pair.
  return fabs(b.x[2] - a.x[2])*pow(10., -log10r);
}
