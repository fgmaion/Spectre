int compare_randESR_PkCorrections(){
    double pk, WindowFunc;

    double rand_shot = 0.0, gal_shot = 0.0, rand_ESR_shot = 0.0;

    double* rand_iESR;

    // inverse ESR for randoms. 
    rand_iESR = malloc(rand_number*sizeof(double));


    sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_xyz_0.6_0.9_Nagoya_v6_Samhain_GranettESR.cat", root_dir, fieldFlag);

    inputfile = fopen(filepath, "r");

    for(j=0; j<rand_number; j++)  fscanf(inputfile, "%lf \t %*lf", &rand_iESR[j]);
    
    fclose(inputfile);


    // shot noise from randoms cat. calculation.    
    for(j=0; j<rand_number; j++){
        if(rand_accept[j]  == true)  rand_shot       += pow(rand_weight[j], 2.);
    } 
                    
    rand_shot     *= alpha*alpha;

    for(j=0; j<rand_number; j++){
      if(rand_accept[j]    == true)    rand_ESR_shot += pow(rand_weight[j], 2.)*rand_iESR[j];
    }

    rand_ESR_shot *= alpha;

    // shot noise from galaxies cat. calculation, including angular sampling. sum over galaxies, so no alpha scaling. 
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true)  gal_shot  += pow(fkp_galweight[j], 2.)/sampling[j];
    }

    printf("\n\nShot noise contributions: Randoms %.4lf, Galaxies %.4lf: ESR weighted randoms: %.4lf", rand_shot, gal_shot, rand_ESR_shot);
    
    return 0;
}

