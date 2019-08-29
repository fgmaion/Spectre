/*
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


int randoms_Sphere(double maxGals, double radius){
    printf("\nCreating sphere of randoms.");

    double x, y, z, r2;
    double GalNumber = 0.;

    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_randoms_500_sphere_%.2e.dat", root_dir, maxGals); 

    output = fopen(filepath, "w");    	
    
    while(GalNumber<maxGals){
        x  = 500.*(gsl_rng_uniform(gsl_ran_r) -0.5) + 250.;
        y  = 500.*(gsl_rng_uniform(gsl_ran_r) -0.5) + 250.;
        z  = 500.*(gsl_rng_uniform(gsl_ran_r) -0.5) + 250.;

        r2 = pow(x - 250., 2.) + pow(y - 250., 2.) + pow(z - 250., 2.);

        if(r2 <= radius*radius){
           fprintf(output, "%e \t %e \t %e \n", x, y, z);
            
           GalNumber += 1.0;
        }
   }
   
   fclose(output);
    
   return 0;
}


int randoms_anisoGauss(double maxGals, double a, double b, double c){
    printf("\nCreating anisotropic Gaussian distribution of randoms.");

    double x, y, z;
    double GalNumber = 0.;

    sprintf(filepath, "%s/Data/stacpolly/poissonSampled_randoms_500_anisoGauss_%.2e.dat", root_dir, maxGals); 

    output = fopen(filepath, "w");    	
    
    while(GalNumber<maxGals){
        x  = 250. + rand_gaussian(gsl_ran_r, a);
        y  = 250. + rand_gaussian(gsl_ran_r, b);
        z  = 250. + rand_gaussian(gsl_ran_r, c);
        
        if((x >= 0.) && (x <= 500.) && (y >= 0.) && (y <= 500.) && (z >= 0.) && (z <= 500.)){
            fprintf(output, "%e \t %e \t %e \n", x, y, z);
            
            GalNumber += 1.0;
        }
    }
   
    fclose(output);
    
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


int knownGRF_mask(){
    sprintf(filepath, "%s/Data/SpectralDistortion/GRF_mask_MonoAndQuad.dat", root_dir);
    
    inputfile     = fopen(filepath, "r");          

    for(j=0; j<n0*n1*n2; j++)  fscanf(inputfile, "%le", &surveyMask[j]);

    fclose(inputfile);

    return 0;
}
*/

int knownGRF_mask_smallCell(){
    sprintf(filepath, "%s/Data/SpectralDistortion/GRF_mask_MonoAndQuad_CellSize_2.00_papercheck.dat", root_dir);
    
    inputfile     = fopen(filepath, "r");          

    for(j=0; j<n0*n1*n2; j++)  fscanf(inputfile, "%le", &surveyMask[j]);

    fclose(inputfile);

    return 0;
}


int randoms_Cube(int maxGals){
    printf("\nCreating cube of randoms.");

    rand_number = maxGals;
    
    assign_randmemory(); 

    for(j=0; j<maxGals; j++){
      rand_x[j] = gsl_rng_uniform(gsl_ran_r)*1000.;
      rand_y[j] = gsl_rng_uniform(gsl_ran_r)*1000.;
      rand_z[j] = gsl_rng_uniform(gsl_ran_r)*1000.;
    
      // no redshift dependence to nbar. 
      rand_chi[j]     = 500.;
      rand_weight[j]  =   1.;
      rand_accept[j]  = true;
    }
    
    accepted_rand = rand_number;
    
    return 0;
}


int randoms_inenvironment(char filepath[]){
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    
    // Column 0: x co-ordinate                                                                                                     
    // Column 1: y co-ordinate                                                                                           
    // Column 2: z co-ordinate                                                                                             

    // Calculate number of in the catalogue (line number);

    ch          = 0;
    rand_number = 0;
    
    do{
        ch = fgetc(inputfile);     
           
        if(ch == '\n')
       	  rand_number += 1;
    } while(ch != EOF);

    printf("\n\nNumber of randoms in catalogue:  %d", rand_number);

    rewind(inputfile);
    
    assign_randmemory(); 

    for(j=0; j<rand_number; j++){
      fscanf(inputfile, "%lf \t %lf \t %lf \n", &rand_x[j], &rand_y[j], &rand_z[j]);
    
      // no redshift dependence to nbar. 
      rand_chi[j]     = 500.;
      rand_weight[j]  =   1.;
      rand_accept[j]  = true;
    }
    
    fclose(inputfile);

    accepted_rand = rand_number;
    
    return 0;
}
