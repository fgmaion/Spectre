int CovarianceMatrix(int mockNumber, int kBinNumb, int hiMultipoleOrder){
    // Multipoles is a [CatalogNumber][hiMultipoleOrder][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 

    assignCovMat(mockNumber, kBinNumb, hiMultipoleOrder);

    for(k=0; k<mockNumber; k++){
            sprintf(filepath, "%s/Data/Multipoles/Multipoles_zCube_clipThreshold_1.0e+03_subVol_%d_kbin_0.010_000.dat", root_dir, k);

            inputfile = fopen(filepath, "r");

            for(i=0; i<kBinNumb-1; i++){
	            fscanf(inputfile, "%lf \t %lf \t %lf \t %*d \t %*lf \t %*lf \n", &kMultipoles[i], &Multipoles[k][0][i], &Multipoles[k][1][i]);
            }

            fclose(inputfile);    
    }


    for(k=0; k<(kBinNumb-1); k++){
        for(j=0; j<hiMultipoleOrder; j++){
            for(i=0; i<mockNumber; i++){
                MeanMultipoles[j][k] += (1./mockNumber)*Multipoles[i][j][k];
            }
        }
    }
    
    
    for(k=0; k<kBinNumb-1;k++){
        for(i=0; i<mockNumber; i++){
            // Multipoles is a [CatalogNumber][hiMultipoleOrder][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
            Multipoles[i][0][k] -= MeanMultipoles[0][k];
            Multipoles[i][1][k] -= MeanMultipoles[1][k];
        }
    }
    
    // Covariance is an N x N matrix, where N corresponds to hiMultipoleOrder*(kBinNumb-1), here hiMultipoleOrder is due to Mono-Mono, Mono-Quad, Quad-Quad, etc... elements. Here hex-blah elements are 
    // ignored. 
    
    for(j=0; j<(kBinNumb-1); j++){
        k = j;
        
        // for(k=0; k<(kBinNumb-1); k++){
            for(i=0; i<mockNumber; i++){
              // Multipoles is a [CatalogNumber][hiMultipoleOrder][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
              
              // Mono-Mono elements.
              Covariance[j][k]                               +=  (1./mockNumber)*Multipoles[i][0][j]*Multipoles[i][0][k];   
              
              // Quad-Quad elements.
              Covariance[(kBinNumb-1) + j][(kBinNumb-1) + k] +=  (1./mockNumber)*Multipoles[i][1][j]*Multipoles[i][1][k];   
                     
              // Mono-Quad.              
              // Covariance[j][(kBinNumb-1) + k]                +=  (1./mockNumber)*Multipoles[i][0][j]*Multipoles[i][1][k];   
                           
              // Quad-Mono.               
              // Covariance[(kBinNumb-1) + j][k]                +=  (1./mockNumber)*Multipoles[i][1][j]*Multipoles[i][0][k];   
            }
        //}
    }
	  	     
	  	     
    sprintf(filepath, "%s/Data/Covariance/zCube_clipThreshold_1.0e+03_subVols_Covariance.dat", root_dir);

    output = fopen(filepath, "w"); 

    for(k=0; k<2*(kBinNumb-1); k++){
        for(j=0; j<2*(kBinNumb-1); j++){        
            fprintf(output, "%e \t", fabs(Covariance[j][k]));                  
        }
    
        fprintf(output, "\n");
    }

    fclose(output);
    
    return 0;
}
