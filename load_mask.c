int load_homogeneous_rands_window(int load, double sampling){
    // sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_parent_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);
    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim);
    
    inputfile   = fopen(filepath, "r");

    ch          = 0;
    rand_number = 0;

    do{
        ch = fgetc(inputfile);
        
        if(ch == '\n')
            rand_number += 1;
    } while(ch != EOF);
    
    // assumes catalogue is in a random order.
    lowerSampling_randomisedCatalogue(sampling);

    printf("\n\nsampling %.2e, %d randoms number", sampling, rand_number);

    if(load == 1){
        rewind(inputfile);

        assign_randmemory();

        for(j=0; j<rand_number; j++)   fscanf(inputfile, "%le \t %le \t %le \t %le \t %le \t %le", &rand_ra[j], &rand_dec[j], &rand_chi[j], &rand_x[j], &rand_y[j], &rand_z[j]);
    }

    fclose(inputfile);
    
    printf("\n\nrandoms, celestial co-ordinates.");
    
    printf("\n%e \t %e", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\n%e \t %e", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\n%e \t %e", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    StefanoReflection(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
    
    StefanoRotated(rand_number, CentreRA, CentreDec, rand_x, rand_y, rand_z);
    
    printf("\n\nStefano basis, randoms co-ordinates.");
    
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    // assign memory for grid representation of mask, Cell_SurveyLimitsMask. 
    prep_mask();    
    
    for(j=0; j<rand_number; j++){
        boxlabel = boxCoordinates(rand_x, rand_y, rand_z, j);
    
        surveyMask[boxlabel] = 1.0;
    }
    
    // free if not pair counting window. 
    // free(rand_x);
    // free(rand_y);
    // free(rand_z);
    
    return 0;
}


int lowerSampling_randomisedCatalogue(double sampling){
    rand_number = (int) ceil(rand_number*sampling);

    return 0;
}


int add_fkp(){
    double xx, yy, zz, chi, nbar;

    for(k=0; k<n0; k++){
        for(j=0; j<n1; j++){
            for(i=0; i<n2; i++){
                Index = k*n1*n2 + j*n2 + i;

                xx    = (i+0.5)*CellSize;
                yy    = (j+0.5)*CellSize;
                zz    = (k+0.5)*CellSize;
                
                chi   = invert_StefanoBasis(CentreRA, CentreDec, &xx, &yy, &zz);
                    
                nbar  = interp_nz(chi);
                    
                surveyMask[Index] *= nbar*fkpPk/(1. + nbar*fkpPk);
            }
        }
    }

    return 0;
}


int calc_DigitalAmplitudeCorrection(){  
    fkpWeightedVolume  = 0.0;
    fkpSqWeightsVolume = 0.0;
    
    // eqn. 16.120, pg 524., P(k) overestimated due to larger density of states. 
    for(j=0; j<n0*n1*n2; j++) fkpWeightedVolume  += surveyMask[j]; 

    for(j=0; j<n0*n1*n2; j++) fkpSqWeightsVolume += pow(surveyMask[j], 2.);
    
    fkpWeightedVolume                            *= CellVolume;
    fkpSqWeightsVolume                           *= CellVolume;
    
    printf("\n\nDigital estimate");
    
    printf("\nFKP     weighted volume:     %e    [Total Volume]", fkpWeightedVolume/TotalVolume);
    printf("\nFKP Sq. weighted volume:     %e    [Total Volume]",   fkpSqWeightsVolume/TotalVolume);
    
    return 0;
}


int calc_AnalogAmplitudeCorrection(){  
    // assuming unit fkp weights, currently. 
    fkpWeightedVolume  = sqdegs2steradians(W1area)*(pow(hiChi, 3.) - pow(loChi, 3.))/3.;
    
    // binary mask. 
    fkpSqWeightsVolume = fkpWeightedVolume;
    
    printf("\n\nAnalog estimate.");
        
    printf("\nFKP     weighted volume:     %e    ", fkpWeightedVolume);
    printf("\nFKP Sq. weighted volume:     %e    ", fkpSqWeightsVolume);
        
    return 0;
}
