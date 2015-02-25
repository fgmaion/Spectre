#define NCHAR 40
#define NR_END 1
#define FREE_ARG char*
#define deg2rad M_PI/180

int ct_lines(char *file){
  FILE *fic;
  
  int n;
  
  fic = fopen(file,"r");
  
  if(fic){
    n=0;
    
    while(!feof(fic)){
      fscanf(fic,"%*[^\n]\n");
      n++;
    }
    
    fclose(fic);
    
    return n;
  } 
  
  else{
    fprintf(stderr,"Can't read %s file!\n",file);
    
    return 0;
  }
}


int getStrings(char *line, char *strings, char *delimit, size_t *N){
  int i,j,begin,length;
  
  if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
  
  i = 0;
  j = 0;
  while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
    begin = i;
    while((line[i] == *delimit || line[i] == '\t') && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
    begin = i;
    while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
    length = i - begin;
    if(length > 0){
      strncpy(strings+NCHAR*j,&line[begin],length);
      strcpy(strings+NCHAR*j+length,"\0");
      j++;
    }
  }
  
  (*N) = j;
  
  if(*N > 0) return 1;
  else return 0;
}


int load_catSpecmask(int mock){
    // HOD mocks with spectroscopic mask, not angular density dependent spoc. 
    if(mock<10)        sprintf(mock_specmask, "%s/mocks_W1_v8.0_500/mock_00%d_specmask.dat", vipersHOD_dir, mock);
    else if(mock<100)  sprintf(mock_specmask, "%s/mocks_W1_v8.0_500/mock_0%d_specmask.dat",  vipersHOD_dir, mock);
    else               sprintf(mock_specmask, "%s/mocks_W1_v8.0_500/mock_%d_specmask.dat",   vipersHOD_dir, mock);
  
    // number of objects in catalogue.
    ch         = 0;
    nx         = 0;
    
    inputfile = fopen(mock_specmask, "r");
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  nx += 1;
    } while(ch != EOF);

    rewind(inputfile);

    rax  = realloc(rax, nx*sizeof(double));
    decx = realloc(decx, nx*sizeof(double));

    printf("\n\n%d objects in spectroscopic mask catalogue.", nx);
  
    for(j=0; j<nx; j++)  fscanf(inputfile, "%*d \t %le \t %le \t %*le \t %*le \t %*le \n", &rax[j], &decx[j]);

    // for(j=0; j<10; j++)  printf("\n%le \t %le", rax[j], decx[j]);

    fclose(inputfile);

    return 0;
}


int load_catSpec(int mock){
    // HOD mocks with spectroscopic mask, not angular density dependent spoc. 
    if(mock<10)        sprintf(mock_spec, "%s/mocks_W1_v8.0_500/mock_00%d_spec.dat", vipersHOD_dir, mock);
    else if(mock<100)  sprintf(mock_spec, "%s/mocks_W1_v8.0_500/mock_0%d_spec.dat",  vipersHOD_dir, mock);
    else               sprintf(mock_spec, "%s/mocks_W1_v8.0_500/mock_%d_spec.dat",   vipersHOD_dir, mock);
  
    // number of objects in catalogue.
    ch         = 0;
    n          = 0;
    
    inputfile = fopen(mock_spec, "r");
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  n   += 1;
    } while(ch != EOF);

    rewind(inputfile);
    
    printf("\n%d objects in spoc catalogue.", n);

    id  = realloc(id,  n*sizeof(double));
    ra  = realloc(ra,  n*sizeof(double));
    dec = realloc(dec, n*sizeof(double));
  
    for(j=0; j<n; j++)  fscanf(inputfile, "%d \t %le \t %le \t %*le \t %*le \t %*le \n", &id[j], &ra[j], &dec[j]);

    // for(j=0; j<10; j++) printf("\n%.6lf \t %.6lf", ra[j], dec[j]);

    fclose(inputfile);

    return 0;
}


int create_polygonMask(int mock){
    // 100 arc seconds by 60. 50x30
    double half_ra  = 30./60./60.; // degrees
    double half_dec = 50./60./60.; // degrees 
  
    // no restriction on order of vertices in polygon format. 
  
    if(mock<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_00%d_W1_spec_polygons.dat", root_dir, mock);
    else if(mock<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_0%d_W1_spec_polygons.dat",  root_dir, mock);
    else               sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_%d_W1_spec_polygons.dat",   root_dir, mock);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<n; i++){
      // rectangle. 
      fprintf(output, "polygon(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f) # id %d \n",   ra[i]-half_ra,  dec[i]+half_dec, 
                                                                                              ra[i]+half_ra,  dec[i]+half_dec,
                                                                                              ra[i]+half_ra,  dec[i]-half_dec,
                                                                                              ra[i]-half_ra,  dec[i]-half_dec,
                                                                                              id[i]
                                                                                             );
    }                                                        
    
    fclose(output);
    
    return 0;
}


int load_polygonMask(int mock){
    if(mock<10)        sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_00%d_W1_spec_polygons.dat", root_dir, mock);
    else if(mock<100)  sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_0%d_W1_spec_polygons.dat",  root_dir, mock);
    else               sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_500_polygons/mock_%d_W1_spec_polygons.dat",   root_dir, mock);
  
    inputfile     = fopen(maskfile, "r");  
  
    ch         = 0;
    Npolygons  = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Npolygons += 1;
    } while(ch != EOF);

    fclose(inputfile);

    polymask = init_poly(maskfile, outmask);

    return 0; 
}


int write_mockweights(int mock){
    if(mock<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_spocWeights/mock_00%d_W1_spec_weights.dat", root_dir, mock);
    else if(mock<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_spocWeights/mock_0%d_W1_spec_weights.dat",  root_dir, mock);
    else               sprintf(filepath, "%s/Venice/mocks_W1_v8.0_500_spocWeights/mock_%d_W1_spec_weights.dat",   root_dir, mock);
    
    output = fopen(filepath, "w");
  
    for(i=0; i<n; i++){  
        if((spoc_counts[i]/spoc_weights[i] < 0.) || (spoc_counts[i]/spoc_weights[i] > 1.)){ 
            printf("\n\nError in spec. weight calculation.");
            
            printf("\n\n%d \t %e \t %e", i, spoc_counts[i], spoc_weights[i]);
            
            return 1;
        }
    
        fprintf(output, "%d \t %lf \t %lf \t %lf \n", id[i], ra[i], dec[i], spoc_counts[i]/spoc_weights[i]);
    }

    fclose(output);

    return 0;
}


int polymaskCalc(){    
  // annoying, but relatively slow, memory accretion. sub sample mocks for each run. 
  for(ii=150; ii<307; ii++){  
      printf("\n\nCalculating SPOC weights for mock %d", ii);
  
      load_catSpecmask(ii);
      
      load_catSpec(ii);
      
      // Create polygon mask, rectangles surrounding spec galaxies. 
      // create_polygonMask(ii);
      
      // Polygon mask, rectangles surrounding spec galaxies. no restriction on order of vertices in polygon format. 
      load_polygonMask(ii);
    
      // load_Nagoya_v4();
      
      // weight for each galaxy in spec catalogue. 
      spoc_weights     = realloc(spoc_weights, n*sizeof(double));
      spoc_counts      = realloc(spoc_counts,  n*sizeof(double));  
      
      for(j=0; j<n; j++){  
        spoc_weights[j] = 0.0;
        spoc_counts[j]  = 0.0;
      }
      
      for(i=0; i<n; i++){
        // count n objects in spoc catalogue.
        addgal_overlapPoly_counts(polymask, outmask, ra[i]*deg2rad, dec[i]*deg2rad, spoc_counts);
    
        // only valid for non-overlapping polygons, returns on first polygon, ignoring the rest. 
        // poly_id = idPoly(polymask, outmask, ra[i]*deg2rad, dec[i]*deg2rad);
        // 
        // if(poly_id != -1)   spoc_counts[poly_id] += 1.;
      }
      
      for(i=0; i<nx; i++){
        // count nx objects in spectroscopic mask catalogue. nx > n.
        addgal_overlapPoly_counts(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad, spoc_weights);
    
        // invalid for overlapping polygons, returns on first polygon - ignoring the rest. 
        // poly_id = idPoly(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad);
      }
      
      write_mockweights(ii);
      
      // free memory for polygons catalogue.
      for(i=0; i<Npolygons; i++){
        free(polysAll[i].x);
        free(polysAll[i].y);
    
        free(polysAll[i].xmin);
        free(polysAll[i].xmax);
      }
        
      free_NodeP(polymask);
  }
  
  return 0;
}
