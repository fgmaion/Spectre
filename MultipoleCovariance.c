int CovarianceMatrix(int mockNumber){
    // Multipoles is a [CatalogNumber][hiMultipoleOrder][chiSq_kmaxIndex] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 

    // Retrieve the number of rows necessary before k value is greater than kmax for chi sq. evaluation. 
    sprintf(filepath, "%s/Data/Multipoles/zCube_xVel/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_001.dat", root_dir);

    inputfile       = fopen(filepath, "r");
    
    // Retrieve the number of lines in the input file. 
    ch         = 0;
    lineNo     = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch    == '\n')
       	  lineNo += 1;
    } while (ch  != EOF);

    rewind(inputfile);

    for(i=0; i<lineNo; i++){
        fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*d \t %*lf \t %*lf \n", &Interim);
    
        if(Interim > ChiSqEval_kmax){
            chiSq_kmaxIndex = i -1;
        
            break;
        }
    }
    
    fclose(inputfile);
    
    printf("\n\nkmax limit for ChiSq: %d", chiSq_kmaxIndex);

    assignCovMat(mockNumber, chiSq_kmaxIndex, 2);
    
    // Be careful with 0 or 1 for the starting value !!
    for(k=0; k<mockNumber; k++){
            if(k<10)  sprintf(filepath, "%s/Data/Multipoles/zCube_xVel/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_00%d.dat", root_dir, k);
            else      sprintf(filepath, "%s/Data/Multipoles/zCube_xVel/Multipoles_zCube_xvel_clipThreshold_1.0e+03_fullCube_kbin_0.010_0%d.dat",  root_dir, k);
            
            // if(k<10)  sprintf(filepath, "%s/Data/Del2k/GaussCube/midK_Del2k_GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_kInterval_0.01_00%d.dat", root_dir, k);
            // else      sprintf(filepath, "%s/Data/Del2k/GaussCube/midK_Del2k_GaussCube_Realisation_clipThreshold_1.0e+03_fullCube_kInterval_0.01_0%d.dat",  root_dir, k);
            
            // if(k<10)  sprintf(filepath, "%s/Data/Del2k/GaussCube_BootStrap/midK_Del2k_GaussCube_BootStrap_clipThreshold_1.0e+03_fullCube_kInterval_0.01_00%d.dat", root_dir, k);
            // else      sprintf(filepath, "%s/Data/Del2k/GaussCube_BootStrap/midK_Del2k_GaussCube_BootStrap_clipThreshold_1.0e+03_fullCube_kInterval_0.01_0%d.dat",  root_dir, k);
            
            inputfile = fopen(filepath, "r");
            
            // printf("\n\nInput of mock number: %d", mockNumber);

            for(i=0; i<chiSq_kmaxIndex; i++){
	            fscanf(inputfile, "%lf \t %lf \t %lf \t %d \t %*lf \t %*lf \n", &kMultipoles[i], &Multipoles[k][0][i], &Multipoles[k][1][i], &ModeNumber[i]);
	       
                // printf("\n%e \t %e \t %e \t %d", kMultipoles[i], Multipoles[k][0][i], Multipoles[k][1][i], ModeNumber[i]);
            }

            // printf("\n\n");            

            fclose(inputfile);    
    }  
    
    for(k=0; k<chiSq_kmaxIndex; k++){
        for(j=0; j<hiMultipoleOrder; j++){
            // Be careful with 0 or 1 for the starting value !!
            for(i=0; i<mockNumber; i++){
                MeanMultipoles[j][k] += (1./mockNumber)*Multipoles[i][j][k];
            }
        }
    }
    
    printf("\n\nMean multipoles: \n");
    
    for(k=0; k<chiSq_kmaxIndex; k++)  printf("\n%e \t %e \t %e", kMultipoles[k], MeanMultipoles[0][k], MeanMultipoles[1][k]);
    
    for(k=0; k<chiSq_kmaxIndex; k++){
        // Be careful with 0 or 1 for the starting value !!
        for(i=0; i<mockNumber; i++){
            // Multipoles is a [CatalogNumber][hiMultipoleOrder][chiSq_kmaxIndex] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
            Multipoles[i][0][k] -= MeanMultipoles[0][k];
            Multipoles[i][1][k] -= MeanMultipoles[1][k];
        }
    }

    // Covariance is an N x N matrix, where N corresponds to hiMultipoleOrder*chiSq_kmaxIndex, here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
    
    for(j=0; j<chiSq_kmaxIndex; j++){    
        for(k=0; k<chiSq_kmaxIndex; k++){
            for(i=0; i<mockNumber; i++){
              // Multipoles is a [CatalogNumber][hiMultipoleOrder][chiSq_kmaxIndex] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
              
              // Mono-Mono elements.
              Covariance[j][k]                                                 +=  (1./mockNumber)*Multipoles[i][0][j]*Multipoles[i][0][k];   
              
              // Quad-Quad elements.
              Covariance[chiSq_kmaxIndex + j][chiSq_kmaxIndex + k]             +=  (1./mockNumber)*Multipoles[i][1][j]*Multipoles[i][1][k];   
                     
              // Mono-Quad.              
              Covariance[j][chiSq_kmaxIndex + k]                               +=  (1./mockNumber)*Multipoles[i][0][j]*Multipoles[i][1][k];   
                           
              // Quad-Mono.               
              Covariance[chiSq_kmaxIndex + j][k]                               +=  (1./mockNumber)*Multipoles[i][1][j]*Multipoles[i][0][k];   
            }
        }
    }
    
    return 0;
}
