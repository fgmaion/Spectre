int AnisoGauss(double a, double b, double c){
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 0.0;
    
    double rx, ry, rz;
    
    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                rz   = CellSize*(k-n0/2);   
                ry   = CellSize*(j-n1/2);
                rx   = CellSize*(i-n2/2);
            
                Index = k*n1*n2 + j*n2 + i;
            
                Cell_SurveyLimitsMask[Index] = exp(-0.5*pow(rx/a, 2.))*exp(-0.5*pow(ry/b, 2.))*exp(-0.5*pow(rz/c, 2.));
            }
        }
    }
    
    return 0;
}


int FullCube(){
    // bool density = 1.0 everywhere for full cube.  Otherwise set to 0.0 initially. 
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 1.0;
    
    return 0;
}


int EmbeddedCube(int width){
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 0.0;

    for(k=width; k<n0-width; k++){
        for(j=width; j<n1-width; j++){
            for(i=width; i<n2-width; i++){
                Index = k*n1*n2 + j*n2 + i;
                
                Cell_SurveyLimitsMask[Index] = 1.0;
            }
        }
    }

    return 0;
}


int Spherical(float radius){
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 0.0;

    float r2 = 0.0;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                r2 = pow(CellSize*(k-n0/2), 2.) + pow(CellSize*(j-n1/2), 2.) + pow(CellSize*(i-n2/2), 2.);     
            
                Index = k*n1*n2 + j*n2 + i;
            
                if(r2<=radius*radius){
                    Cell_SurveyLimitsMask[Index] = 1.0;
                }
            }   
        }
    }

    return 0;
}


int Gaussian(float sigma){
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 0.0;

    float r2 = 0.0;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                r2 = pow(CellSize*(k-n0/2), 2.) + pow(CellSize*(j-n1/2), 2.) + pow(CellSize*(i-n2/2), 2.);     
            
                Index = k*n1*n2 + j*n2 + i;
            
                Cell_SurveyLimitsMask[Index] = exp(-0.5*r2/pow(sigma, 2.));
            }
        }
    }

    return 0;
}


int PencilBeamSurvey(int xlow, int xhi, int ylow, int yhi){
    for(j=0; j<n0*n1*n2; j++) Cell_SurveyLimitsMask[j] = 0.0;

    for(k=0; k<n0; k++){
        for(j=ylow; j<yhi; j++){
            for(i=xlow; i<xhi; i++){
                Index = k*n1*n2 + j*n2 + i;
                
                Cell_SurveyLimitsMask[Index] = 1.0;
            }
        }
    }

    return 0;
}


int VIPERS_mask(){
    sprintf(filepath, "%s/Data/SpectralDistortion/VIPERS_mask.dat", root_dir);
    
    inputfile     = fopen(filepath, "r");          

    for(j=0; j<n0*n1*n2; j++)  fscanf(inputfile, "%le", &Cell_SurveyLimitsMask[j]);

    fclose(inputfile);
    
    int nonEmptyCells = 0;
    
    for(j=0; j<n0*n1*n2; j++){
        if(fabs(Cell_SurveyLimitsMask[j]) > 0.0){  
            nonEmptyCells      +=   1;  
            Cell_VIPERSbools[j] = 1.0;
        }
    }
    
    printf("\n\nNon-empty cells in window fn.: %d", nonEmptyCells);   
    
    return 0;
}


int VIPERS_Binarymask(){
    sprintf(filepath, "%s/Data/SpectralDistortion/VIPERS_maskBinary.dat", root_dir);
    
    inputfile     = fopen(filepath, "r");          

    for(j=0; j<n0*n1*n2; j++)  fscanf(inputfile, "%le", &Cell_SurveyLimitsMask[j]);

    fclose(inputfile);
    
    int nonEmptyCells = 0;
    
    for(j=0; j<n0*n1*n2; j++){
        if(fabs(Cell_SurveyLimitsMask[j]) > 0.0){  
            nonEmptyCells      +=   1;  
            Cell_VIPERSbools[j] = 1.0;
        }
    }
    
    printf("\n\nNon-empty cells in window fn.: %d", nonEmptyCells);   
    
    return 0;
}
