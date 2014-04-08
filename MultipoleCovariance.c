int CovarianceMatrix(){
    double***     Multipoles;

    int mockNumber = 64;
    int kBinNumb   = 20;

    // Multipoles is a [CatalogNumber][3][kBinNumb] object, where [x][][] corresponds to mock number, [][x][] corresponds to x e [Mono, Quad, Hex], [][][x] mid or mean k bin. 
    Multipoles                                             = (double ***) malloc(mockNumber*sizeof(*Multipoles));
        
    for(j=0; j<mockNumber; j++) Multipoles[j]              = (double  **) malloc(3*sizeof(**Multipoles));
    
    for(j=0; j<mockNumber; j++){
        for(k=0; k<3; k++)      Multipoles[j][k]           = (double   *) malloc((kBinNumb-1)*sizeof(***Multipoles));
    }

    for(i=0; i<mockNumber; i++){
        for(j=0; j<3; j++){
            for(k=0; k<kBinNumb-1; k++){
                Multipoles[i][j][k] = 0.0;
            }
        }
    }
    
    kMultipoles                                            = (double  *) malloc((kBinNumb-1)*sizeof(*kMultipoles));

    MeanMultipoles                                         = (double **) malloc(3*sizeof(*MeanMultipoles));    
    
    for(j=0; j<3; j++)  MeanMultipoles[j]                  = (double  *) malloc((kBinNumb-1)*sizeof(**MeanMultipoles));
    
    for(j=0; j<3; j++){
        for(k=0; k<(kBinNumb-1); k++){
            MeanMultipoles[j][k] = 0.0;
        }
    }
    
    int counter = 0;
    
    for(k=0; k<8; k++){
       for(j=0; j<8; j++){
            sprintf(filepath, "%s/Data/ClippedPk/zSpace/2Dpk/Observed_Multipoles_Clipped_zPencilBeamCube_Jenkins1.0_xtrans_%d.00_ytrans_%d.00_kbin_0.040_000.dat", root_dir, 120*k, 120*j);

            inputfile = fopen(filepath, "r");

            for(i=0; i<kBinNumb-1; i++){
	            fscanf(inputfile, "%lf \t %lf \t %lf \t %*d \n", &kMultipoles[i], &Multipoles[counter][0][i], &Multipoles[counter][1][i]);
            }

            fclose(inputfile);    
            
            counter += 1; 
        }
    }

    for(k=0; k<(kBinNumb-1); k++){
        for(j=0; j<3; j++){

            for(i=0; i<mockNumber; i++){
                MeanMultipoles[j][k] += (1./mockNumber)*Multipoles[i][j][k];
            }
        }
        
        // printf("\n%e \t %e \t %e", MeanMultipoles[0][k], MeanMultipoles[1][k], MeanMultipoles[2][k]);
    }

    // Covariance is an N x N matrix, where N corresponds to 3*(kBinNumb-1), here 3 is due to Mono-Mono, Mono-Quad, Mono-Hex, Quad-Quad, etc... elements.     
    Covariance                                                = (double **) malloc(2*(kBinNumb-1)*sizeof(*Covariance));
    
    for(j=0; j<2*(kBinNumb-1); j++) Covariance[j]             = (double  *) malloc(2*(kBinNumb-1)*sizeof(**Covariance));

    for(k=0; k<2*(kBinNumb-1); k++){
        for(j=0; j<2*(kBinNumb-1); j++){
            Covariance[j][k] = 0.0;
        }
    }

    int MonoCount; 
    int QuadCount;
    
    for(MonoCount=0; MonoCount<2; MonoCount++){
        for(QuadCount=0; QuadCount<2; QuadCount++){
            for(j=0; j<(kBinNumb-1); j++){
                for(k=0; k<(kBinNumb-1); k++){
                    for(i=0; i<mockNumber; i++){
                           Covariance[MonoCount*(kBinNumb-1) + j][QuadCount*(kBinNumb-1) + k] +=  (1./mockNumber)*(Multipoles[i][MonoCount][j] - MeanMultipoles[MonoCount][j])*(Multipoles[i][QuadCount][k] - MeanMultipoles[QuadCount][k]);
                
                    }
                    
                    // printf("\n%d\t%d\t%e", MonoCount, QuadCount, Covariance[j][k]);
                }
            }
	    }	
	}
     
    sprintf(filepath, "%s/Data/Covariance/Clipped_HODCubeColumns_Covariance.dat", root_dir);

    output = fopen(filepath, "w"); 

    for(j=0; j<2*(kBinNumb-1); j++){
        for(k=0; k<2*(kBinNumb-1); k++){
            fprintf(output, "%e \t", fabs(Covariance[j][k]));                  
        }
    
        fprintf(output, "\n");
    }

    fclose(output);
    
    return 0;
}
