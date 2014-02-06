int apodiseWindowfn(){
    int    NewIndex;
    double Distance;
    double NewWeight;

    int r, s, p;
    
    for(j=0; j<n0*n1*n2; j++)  Cell_ApodiseWeights[j] *= Cell_AppliedWindowFn[j];
 
    for(k=4; k<n0-4; k++){
        for(j=4; j<n1-4; j++){
            for(i=4; i<n2-4; i++){
	            end_nested_loop: Index = k*n1*n2 + j*n2 + i;
                
                if((Cell_AppliedWindowFn[Index] > 0.5) && (Cell_SurveyEdge[Index] < 0.5)){
                    // Find all cells which are adjacent (3D) to an empty cell. Select these cells.
                    
                    for(r=-1; r<2; r++){
                        for(s=-1; s<2; s++){
                            for(p=-1; p<2; p++){
                                NewIndex = (k + r)*n1*n2 + (j + s)*n2 + (i + p);
                        
                                if(Cell_AppliedWindowFn[NewIndex] < 0.5){
                                    Cell_SurveyEdge[Index] = 1.0;
						            goto end_nested_loop;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Cell_SurveyEdge should now contain all edge cells only. 
    printf("\nFound edge cells.");
    
    int bubbleSize;
    
    bubbleSize = 1 + (int) ceil(GibbsSkinDepth);
    
    // Send out a bubble from each edge cell, assign an apodise weight to each cell in the bubble, keep only the largest weight (corresponding to shortest distance to survey edge).
    for(k=4; k<n0-4; k++){
        for(j=4; j<n1-4; j++){
            for(i=4; i<n2-4; i++){
                Index = k*n1*n2 + j*n2 + i;
            
                // Send out a bubble from each edge cell only. 
                if(Cell_SurveyEdge[Index] > 0.5){
                    for(r=-bubbleSize; r<1+bubbleSize; r++){
                        for(s=-bubbleSize; s<1+bubbleSize; s++){
                            for(p=-bubbleSize; p<1+bubbleSize; p++){
                                NewIndex  = (k + r)*n1*n2 + (j + s)*n2 + (i + p);
    
                                Distance  = CellSize*pow(r*r + s*s + p*p, 0.5);
                                
                                if(Cell_ShortDist2edge[NewIndex] > Distance){
                                    Cell_ShortDist2edge[NewIndex] = Distance;
                                }   
                                
                                if(Cell_ShortDist2edge[NewIndex]<=GibbsSkinDepth){
                                     Cell_ApodiseWeights[NewIndex] = 1.0 - cos(0.5*pi*Cell_ShortDist2edge[NewIndex]/GibbsSkinDepth);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //  Now apply the apodisation weights to the applied window function. 
    for(kk=0; kk<n0*n1*n2; kk++)  Cell_ApodiseWeights[kk] *= Cell_AppliedWindowFn[kk];   

    sprintf(filepath, "%s/Data/ra_decCells/xyz_aposideWeights.dat", root_dir);
    output = fopen(filepath, "w");
    
    for(j=0; j<n0*n1*n2; j=j+10){
        if(Cell_AppliedWindowFn[j] > 0.0){
            fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \n", Cell_rotatedXvals[j],  Cell_rotatedYvals[j], Cell_rotatedZvals[j], Cell_ApodiseWeights[j], Cell_ShortDist2edge[j], Cell_SurveyEdge[j]);
        }
    }
    
    fclose(output);

    //  Now apply the apodisation weights to the applied window function. 
    for(kk=0; kk<n0*n1*n2; kk++)  Cell_AppliedWindowFn[kk]  *= Cell_ApodiseWeights[kk];   
    
    return 0;
}
