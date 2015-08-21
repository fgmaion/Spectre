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


int load_WX_Photometry(){
    //** FOR FULL INFO SEE README FILE IN ROOTDIR/W1_SPECTRO_V7_0.
    // Venice GenData on W1_PHOT_newFlag_photMask_selected.txt, double formatting required, be careful ra and dec ordering and col numbering. 
    // File    : W1_PHOT_newFlag_photMask_Nagoya_v6_selected.txt
    // Date    : 2015-07-10 14:01
    // Table   : W1_PHOT
    // Columns : 6
    // Rows    : xxx
    // SQL     : SELECT num, alpha, delta, newFlag, photoMask, spectroMask F
    //           ROM W1_PHOT
    // num	alpha	delta	newFlag	photoMask	spectroMask
    // CHAR	DOUBLE	DOUBLE	INT	TINYINT	TINYINT
    // 	deg	deg	

    // COMMENTS:
    // Selected by NEWFLAG == 1 (meeting i_AB < 22.5 selection). 
    // PHOTMASK == 1, inside Samhain.
    // Ignore SPECTROMASK, specifies Nagoya v4 selection. 
    // Nagoya V6 selected.

    // Photometry catalogue, W1_PHOT.txt in W1_Spectro_V5_0, corresponding to W1_SPECTRO_V7_0.txt
    sprintf(filepath, "%s/W1_Spectro_V7_0/W%d_PHOT_newFlag_photMask_Nagoya_v6_selected.txt",   root_dir, fieldFlag);
  
    // number of objects in catalogue.
    ch         = 0;
    nx         = 0;
    
    inputfile = fopen(filepath, "r");
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  nx += 1;
    } while(ch != EOF);

    rewind(inputfile);

    idx  = realloc(idx,  nx*sizeof(int));
    rax  = realloc(rax,  nx*sizeof(double));
    decx = realloc(decx, nx*sizeof(double));
    zx   = realloc(zx,   nx*sizeof(double));
    
    // printf("\n\n%d objects in spectroscopic mask catalogue.", nx);
  
    for(j=0; j<nx; j++)  fscanf(inputfile, "%d \t %le \t %le \t %*le \t %*le \t %*le \n", &idx[j], &rax[j], &decx[j]);

    printf("\n\nVIPERS: W%d photometry loaded.\n", fieldFlag);

    for(j=0; j<10; j++)  printf("%d \t %.4lf \t %.4lf \n", idx[j], rax[j], decx[j]);

    fclose(inputfile);

    return 0;
}


int load_catSpecmask(int mock){
    // HOD mocks with spectroscopic mask, not angular density dependent spoc. 
  if(mock<10)        sprintf(mock_specmask, "%s/mocks_W%d_v8.0_500/mock_00%d_specmask_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, mock);
  else if(mock<100)  sprintf(mock_specmask, "%s/mocks_W%d_v8.0_500/mock_0%d_specmask_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, mock);
  else               sprintf(mock_specmask, "%s/mocks_W%d_v8.0_500/mock_%d_specmask_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, mock);

    printf("\n\nLoading photometry: %s", mock_specmask);
  
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

    idx  = realloc(idx,  nx*sizeof(int));
    rax  = realloc(rax,  nx*sizeof(double));
    decx = realloc(decx, nx*sizeof(double));
    zx   = realloc(zx,   nx*sizeof(double));
    
    // printf("\n\n%d objects in spectroscopic mask catalogue.", nx);
  
    for(j=0; j<nx; j++)  fscanf(inputfile, "%d \t %le \t %le \t %le \t %*le \t %*le \n", &idx[j], &rax[j], &decx[j], &zx[j]);

    // for(j=0; j<10; j++)  printf("\n%.2lf \t %.2lf \t %.2lf", rax[j], decx[j], zx[j]);

    fclose(inputfile);

    return 0;
}


int load_catSpec(int mock){
    // HOD mocks with spectroscopic mask, not angular density dependent spoc. 
  if(mock<10)        sprintf(mock_spec, "%s/mocks_W%d_v8.0_500/mock_00%d_spec_Nagoya_v6_Samhain.dat", vipersHOD_dir, fieldFlag, mock);
  else if(mock<100)  sprintf(mock_spec, "%s/mocks_W%d_v8.0_500/mock_0%d_spec_Nagoya_v6_Samhain.dat",  vipersHOD_dir, fieldFlag, mock);
  else               sprintf(mock_spec, "%s/mocks_W%d_v8.0_500/mock_%d_spec_Nagoya_v6_Samhain.dat",   vipersHOD_dir, fieldFlag, mock);

  printf("\n\nLoading spectroscopy: %s", mock_spec);

    /** Mock TSR and Mock spec catalogues **/
    // if(mock<10)        sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_00%d_W1_GalSubsample_bypolygon.dat", root_dir, mock);
    // else if(mock<100)  sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_0%d_W1_GalSubsample_bypolygon.dat",  root_dir, mock);
    // else               sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_%d_W1_GalSubsample_bypolygon.dat",   root_dir, mock);
    
    // if(mock<10)        sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_00%d_W1_GalSubsample_nonpoisson_bypolygon.dat", root_dir, mock);
    // else if(mock<100)  sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_0%d_W1_GalSubsample_nonpoisson_bypolygon.dat",  root_dir, mock);
    // else               sprintf(mock_spec, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_%d_W1_GalSubsample_nonpoisson_bypolygon.dat",   root_dir, mock);
  
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

    Acceptanceflag = realloc(Acceptanceflag, n*sizeof(bool));
    
    printf("\n%d objects in spectroscopic mask catalogue.", nx);
  
    // for(j=0; j<n; j++)  fscanf(inputfile, "%d \t %le \t %le \t %*le \t %*le \t %*le \n", &id[j], &ra[j], &dec[j]);
    for(j=0; j<n; j++)  fscanf(inputfile, "%d \t %le \t %le \t %*le \t %*le \t %*le \n", &id[j], &ra[j], &dec[j]);

    for(j=0; j<n; j++)  Acceptanceflag[j] = true;
  
    // for(j=0; j<10; j++) printf("\n%.6lf \t %.6lf", ra[j], dec[j]);

    fclose(inputfile);

    return 0;
}


int load_catMockTSR(int mock){
    // HOD mocks with spectroscopic mask, not angular density dependent spoc. 
    if(mock<10)        sprintf(mock_specmask, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_00%d_W1_TSR.dat", root_dir, mock);
    else if(mock<100)  sprintf(mock_specmask, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_0%d_W1_TSR.dat",  root_dir, mock);
    else               sprintf(mock_specmask, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_%d_W1_TSR.dat",   root_dir, mock);
  
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

    idx  = realloc(idx,  nx*sizeof(int));
    rax  = realloc(rax,  nx*sizeof(double));
    decx = realloc(decx, nx*sizeof(double));
    zx   = realloc(zx,   nx*sizeof(double));
    
    printf("\n\n%d objects in spectroscopic mask catalogue.", nx);
  
    for(j=0; j<nx; j++)  fscanf(inputfile, "%le \t %le \t %le \n", &rax[j], &decx[j], &zx[j]);
    
    // for(j=0; j<10; j++)  printf("\n%le \t %le", rax[j], decx[j]);

    fclose(inputfile);

    return 0;
}


int create_polygonMask(int mock){
    // 100 arc seconds by 60. 50x30
    double half_ra  = 30./60./60.; // degrees
    double half_dec = 50./60./60.; // degrees 
  
    // double half_ra  = 0.025; // degrees
    // double half_dec =  0.05; // degrees 
  
    // no restriction on order of vertices in polygon format. 
    
    // analysis on **mock** catalogues. 
    if(data_mock_flag == 0){
    // if(mock<10)       sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_00%d_W1_mockSpec_nonpoisson_polygons.dat", root_dir, mock);
    // else if(mock<100) sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_0%d_W1_mockSpec_nonpoisson_polygons.dat",  root_dir, mock);
    // else              sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_%d_W1_mockSpec_nonpoisson_polygons.dat",   root_dir, mock);
    
    // if(mock<10)       sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_00%d_W1_spec_polygons.dat", root_dir, mock);
    // else if(mock<100) sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_0%d_W1_spec_polygons.dat",  root_dir, mock);
    // else              sprintf(filepath, "%s/Venice/mocks_W1_v8.0_polygons/mock_%d_W1_spec_polygons.dat",   root_dir, mock);
    
      if(mock<10)      sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_00%d_W%d_spec_polygons.dat", root_dir, fieldFlag, mock, fieldFlag);
      else if(mock<100)sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_0%d_W%d_spec_polygons.dat",  root_dir, fieldFlag, mock, fieldFlag);
      else             sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_%d_W%d_spec_polygons.dat",   root_dir, fieldFlag, mock, fieldFlag);
    }

    printf("\n\nWriting %s", filepath);
    
    // analysis on W1_SPECTRO_V7_0 catalogues. 
    if(data_mock_flag == 1){
        // simple variable assignment. 
        n = Vipers_Num;
     
        sprintf(filepath, "%s/W1_Spectro_V7_0/polygons_W%d_Spectro_V7_0.dat",   root_dir, fieldFlag);
    }
    
    output = fopen(filepath, "w");
    
    // n, id, ra, dec labels spec galaxies. 
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
    // analysis on **mock** catalogues. 
    if(data_mock_flag == 0){
        // if(mock<10)        sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_00%d_W1_spec_polygons.dat", root_dir, mock);
        // else if(mock<100)  sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_0%d_W1_spec_polygons.dat",  root_dir, mock);
        // else               sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_%d_W1_spec_polygons.dat",   root_dir, mock);
  
        // if(mock<10)        sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_00%d_W1_mockTSR_polygons.dat", root_dir, mock);
        // else if(mock<100)  sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_0%d_W1_mockTSR_polygons.dat",  root_dir, mock);
        // else               sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_%d_W1_mockTSR_polygons.dat",   root_dir, mock);
  
        // if(mock<10)        sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_00%d_W1_mockSpec_nonpoisson_polygons.dat", root_dir, mock);
        // else if(mock<100)  sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_0%d_W1_mockSpec_nonpoisson_polygons.dat",  root_dir, mock);
        // else               sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/mock_%d_W1_mockSpec_nonpoisson_polygons.dat",   root_dir, mock);
  
        // if(mock<10)        sprintf(maskfile, "%s/Venice/mocks_W1_v8.0_polygons/one_poly.dat", root_dir, mock);
  
      if(mock<10)        sprintf(maskfile, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_00%d_W%d_spec_polygons.dat", root_dir, fieldFlag, mock, fieldFlag);
      else if(mock<100)  sprintf(maskfile, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_0%d_W%d_spec_polygons.dat",  root_dir, fieldFlag, mock, fieldFlag);
      else               sprintf(maskfile, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/polygons/mock_%d_W%d_spec_polygons.dat",   root_dir, fieldFlag, mock, fieldFlag);
    }
    
    if(data_mock_flag == 1){
        // Analysis on W1_SPECTRO_V7_0
        sprintf(maskfile, "%s/W1_Spectro_V7_0/polygons_W%d_Spectro_V7_0.dat",   root_dir, fieldFlag);
    
        printf("\n");
    }

    printf("\n\nReading file: %s", maskfile);
    
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


int load_Nagoya_v6(){
    sprintf(maskfile, "%s/Venice/venice_3.8.3/allmasks_nagoya/regions_nagoya_6/regions_nagoya_6.0_W1.reg", root_dir);
  
    Npolygons = -99999;
  
    polymask = init_poly(maskfile, outmask);

    return 0; 
}


int localTSR_calc_Data(){    
      printf("\n\nCalculating local TSR weights for W%d_Nagoya_V7_0 data.", fieldFlag);
    
      DataInput();
      
      // load Photometry. i_AB<22.5, photmask and Nagoya v6. selection met. 
      load_WX_Photometry();
      
      // Acceptance criteria to be applied: 
      //   i) Redshift cut
      //  ii) zFlag cut, 2 to 9 inclusive. 
      // iii) photoMask == 1 
      assignAcceptance_WX_SPECTRO_V7();
      
      // Create polygon mask, rectangles surrounding spec galaxies. 
      create_polygonMask(0);
      
      // Polygon mask, rectangles surrounding spec galaxies. no restriction on order of vertices in polygon format. 
      load_polygonMask(0);
      
      // weight for each galaxy in spectroscopic W1_Nagoya_V7_0 catalogue. 
      spoc_weights_spec    = realloc(spoc_weights_spec,   n*sizeof(double));
      spoc_weights_photo   = realloc(spoc_weights_photo,  n*sizeof(double));
      
      for(j=0; j<n; j++){  
         spoc_weights_spec[j]  = 0.0;
        spoc_weights_photo[j]  = 0.0;
      }
      
      for(i=0; i<n; i++){
        if(Acceptanceflag[i] == true){
            // count n objects in spec catalogue.
            addgal_overlapPoly_counts(polymask, outmask, ra[i]*deg2rad, dec[i]*deg2rad, spoc_weights_spec);
        }
      }
      
      printf("\n\n");
      
      for(i=0; i<nx; i++){
        // All galaxies in input photometric catalogue meet selection criteria, i_AB<22.5, photomask and Nagoya v6.
        // count nx objects in spectroscopic mask catalogue. nx > n.
        addgal_overlapPoly_counts(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad, spoc_weights_photo);
      }
         
      write_localTSR(ii);
      
      return 0;
}


int write_localTSR(int mock){
    // analysis on **mock** catalogues. 
    if(data_mock_flag == 0){
        // if(mock<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_mockTSR_weights.dat", root_dir, mock);
        // else if(mock<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_mockTSR_weights.dat",  root_dir, mock);
        // else               sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_mockTSR_weights.dat",   root_dir, mock);
    
        // if(mock<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_mockSpec_nonpoisson_weights.dat", root_dir, mock);
        // else if(mock<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_mockSpec_nonpoisson_weights.dat",  root_dir, mock);
        // else               sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_mockSpec_nonpoisson_weights.dat",   root_dir, mock);
    
        // if(mock<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_00%d_W1_spec_weights.dat", root_dir, mock);
        // else if(mock<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_0%d_W1_spec_weights.dat",  root_dir, mock);
        // else               sprintf(filepath, "%s/Venice/mocks_W1_v8.0_spocWeights/mock_%d_W1_spec_weights.dat",   root_dir, mock);
    
      if(mock<10)        sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_00%d_W%d_spec_weights.dat", root_dir, fieldFlag, mock, fieldFlag);
      else if(mock<100)  sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_0%d_W%d_spec_weights.dat",  root_dir, fieldFlag, mock, fieldFlag);
      else               sprintf(filepath, "%s/W%d_Nagoya_v6_mocks_work/spec_weights/mock_%d_W%d_spec_weights.dat",   root_dir, fieldFlag, mock, fieldFlag);    
    }
    
    if(data_mock_flag == 1){
        // Analysis on W1_SPECTRO_V7_0
        sprintf(filepath, "%s/W1_Spectro_V7_0/spec_weights_W%d_Spectro_V7_0.dat",   root_dir, fieldFlag);    
    }

    printf("\n\nWriting spec weights : %s", filepath);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<n; i++){  
        if(Acceptanceflag[i] == false){
            fprintf(output, "%d \t %lf \n", id[i], 0.);
        }
    
        else{
            if((spoc_weights_spec[i] == 0.) || (spoc_weights_photo[i] == 0.)){ 
                printf("\nAnomaly in spec. weight calculation: %d \t %e \t %e", id[i], spoc_weights_spec[i], spoc_weights_photo[i]);
                
                // return 1;
            }
            
            else{
                // spec/specmask 
                fprintf(output, "%d \t %lf \n", id[i], spoc_weights_spec[i]/spoc_weights_photo[i]);    
            }
        }
    }
    
    fclose(output);
    
    return 0;
}


int localTSR_Calc(){    
  // annoying, but "relatively" slow, memory accretion. sub sample mocks for each run.
  
  int start = 489; 
  
  // for(ii=start; ii<start+50; ii++){  
  for(ii=start; ii<start+100; ii++){  
      printf("\n\nCalculating SPOC weights for mock %d", ii);
  
      load_catSpecmask(ii);
      
      load_catSpec(ii);
      
      // Create polygon mask, rectangles surrounding spec galaxies. 
      create_polygonMask(ii);
      
      // Polygon mask, rectangles surrounding spec galaxies. no restriction on order of vertices in polygon format. 
      load_polygonMask(ii);
      
      // weight for each galaxy in spec catalogue. 
      spoc_weights_spec      = realloc(spoc_weights_spec,   n*sizeof(double));
      spoc_weights_photo     = realloc(spoc_weights_photo,  n*sizeof(double));

      // spoc_counts      = realloc(spoc_counts,   n*sizeof(int));  
      
      for(j=0; j<n; j++){  
         spoc_weights_spec[j]  = 0.0;
        spoc_weights_photo[j]  = 0.0;
        // spoc_counts[j]   = 0;
      }
      
      for(i=0; i<n; i++){
        // count n objects in spec catalogue.
        addgal_overlapPoly_counts(polymask, outmask, ra[i]*deg2rad, dec[i]*deg2rad, spoc_weights_spec);
    
        // only valid for non-overlapping polygons, returns on first polygon, ignoring the rest. 
        // poly_id = idPoly(polymask, outmask, ra[i]*deg2rad, dec[i]*deg2rad);
        // 
        // if(poly_id != -1)   spoc_counts[poly_id] += 1.;
      }
      
      printf("\n\n");
      
      for(i=0; i<nx; i++){
        // count nx objects in spectroscopic mask catalogue. nx > n.
        addgal_overlapPoly_counts(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad, spoc_weights_photo);
    
        // invalid for overlapping polygons, returns on first polygon - ignoring the rest. 
        // poly_id = idPoly(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad);
      }
         
      write_localTSR(ii);
      /*
      // free memory for polygons catalogue.
      for(i=0; i<Npolygons; i++){
        free(polysAll[i].x);
        free(polysAll[i].y);
    
        free(polysAll[i].xmin);
        free(polysAll[i].xmax);
      }
        
      free_NodeP(polymask);
      */
  }
  
  return 0;
}


int sampleGals_byPoly(){
    int    empty_polys = 0, poly_id, minCounts = 9999999;

    double meanCounts;

    for(i=0; i<nx; i++){
      // count nx objects in spectroscopic mask catalogue. nx > n.
      addgal_overlapPoly_counts(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad, spoc_weights_spec);
      
      // returns polygon id in which the galaxy resides. 
      poly_id = idPoly(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad);
      
      // list of ids of each galaxy in each polygon.  
      quad_gal_ids[poly_id][spoc_counts[poly_id]] = i;
      
      // number of galaxies in polygon poly_id.
      spoc_counts[poly_id] += 1;
      
      if(spoc_counts[poly_id] > 900)  printf("\n\nError: insufficient memory assigned for a list of galaxies in each polygon, %d.", spoc_counts[poly_id]);
    }
    
    // original ordering. 
    // for(j=0; j<spoc_counts[0]; j++)  printf("\n%d", quad_gal_ids[0][j]);
    
    for(j=0; j<Npolygons; j++){
        // find polygon with minimum number of counts, polygon must be occupied. 
        if(spoc_counts[j]<minCounts)  if(spoc_counts[j] > 0) minCounts = spoc_counts[j];
    }
    
    // some polygons in Nagoya v4. mask are empty.
    for(i=0; i<Npolygons; i++)  if(spoc_counts[i] == 0) empty_polys += 1;
      
      
    printf("\n\nNumber of non-empty polygons: %d", Npolygons - empty_polys);
    
    printf("\n\nfewest occupied polygon contains %d gals", minCounts);
    
    
    int rand, keep;
    
    // global scope int. 
    NN = minCounts;
    
    // for(j=0; j<spoc_counts[0]; j++)  printf("\n%d", quad_gal_ids[0][j]);
    
    // The idea here is to take a list of galaxies associated to a given quadrant, and 
    // subsample said list until there are minCounts galaxies in the polygons.  Then
    // assign the appropriate weight. 
    
    // mock sampling. 
    for(i=0; i<Npolygons; i++){
        if(spoc_counts[i] > 0){    
            // polygon must be occupied. 
            for(j=0; j<minCounts; j++){
                // randomly pick a galaxy from the list associated to the polygon. 
                // once a galaxy has made the subsample, ignore it.  
                rand = j + gsl_rng_uniform_int(gsl_ran_r, spoc_counts[i] - j);
        
                keep = quad_gal_ids[i][j];
            
                // swap this galaxy to the front of the list, the first minCounts galaxies will be 
                // kept for the subsample catalogue.
                quad_gal_ids[i][j]    = quad_gal_ids[i][rand];
        
                quad_gal_ids[i][rand] = keep;
            }
        }
    }

    // for(j=0; j<spoc_counts[0]; j++)  printf("\n%d", quad_gal_ids[0][j]);
    
    return 0;
}


int write_TSR_subsample(){
    int gal_id;

    // if(ii<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_00%d_W1_GalFullsample_bypolygon.dat", root_dir, ii);
    // else if(ii<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_0%d_W1_GalFullsample_bypolygon.dat",  root_dir, ii);
    // else             sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_%d_W1_GalFullsample_bypolygon.dat",   root_dir, ii);
    
    if(ii<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_00%d_W1_GalSubsample_bypolygon.dat", root_dir, ii);
    else if(ii<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_0%d_W1_GalSubsample_bypolygon.dat",  root_dir, ii);
    else             sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockTSR/mock_%d_W1_GalSubsample_bypolygon.dat",   root_dir, ii);
    
    output = fopen(filepath, "w");
          
    for(i=0; i<Npolygons; i++){
      if(spoc_counts[i]>0){      
        for(j=0; j<NN; j++){  
          // poly_id = idPoly(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad);
          gal_id   =  quad_gal_ids[i][j];
                                                                        
          // associated weight is minCounts/originalCounts.                                                          
          fprintf(output, "%d \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", idx[gal_id], rax[gal_id], decx[gal_id], zx[gal_id], 1.*NN/spoc_counts[i]);
        }
        
        /*
        // full catalogue. associated weight is 0. 
        for(j=NN; j<spoc_counts[i]; j++){            
            gal_id   =  quad_gal_ids[i][j];
        
            // weight 0., did not make subsample catalogue.
            fprintf(output, "%e \t %e \t %e \t %e \n", rax[gal_id], decx[gal_id], zx[gal_id], 0.);
        }*/
      }
    }

    fclose(output);

    return 0;
}


int randoms_mockTSR_weights(){
    double weight;

    int   poly_id;

    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_mockTSR_xyz_%.1f_%.1f.cat", lo_zlim, hi_zlim); 
    
    output = fopen(filepath, "w");
    
    for(i=0; i<rand_number; i++){      
      // returns polygon id in which the galaxy resides. 
      poly_id = idPoly(polymask, outmask, rand_ra[i]*deg2rad, rand_dec[i]*deg2rad);
    
      weight  = (double) NN/spoc_counts[poly_id];
      
      // if(weight > 0.0) printf("\n%e", weight);
      
      fprintf(output, "%le \t %le \t %le \t %le \t %le \t %le \t %le \n", rand_ra[i], rand_dec[i], rand_chi[i], rand_x[i], rand_y[i], rand_z[i], weight);
    }

    fclose(output);

    return 0;
}


int assign_spocweights_memory(){
  // weight for each galaxy in spec catalogue. 
  spoc_weights_spec  = realloc(spoc_weights_spec, Npolygons*sizeof(double));
  spoc_counts        = realloc(spoc_counts,  Npolygons*sizeof(int));  
      
  for(j=0; j<Npolygons; j++){  
    spoc_weights_spec[j] = 0.0;
    spoc_counts[j]  =   0;
  }

  // list of those galaxies which occupy given polygon. 
  quad_gal_ids  = realloc(quad_gal_ids, Npolygons*sizeof(int*));
  
  // limited to have 900 galaxies per quadrant. 
  for(j=0; j<Npolygons; j++)  quad_gal_ids[j] = malloc(20*sizeof(int));

  for(i=0; i<Npolygons; i++){
    for(j=0; j<20; j++){
        quad_gal_ids[i][j] = -99;
    }
  }

  return 0;
}


int mask_area(){
  // calc. of area must be done independently of any other calc. in spec_weights.c, 
  // polygon memory is altered in place. 
  
  // Cannot be used for Nagoya v6 and Samhain due to overlapping polygons / polygons in Samhain
  // outside Nagoya v6.
  load_Nagoya_v6();
  
  // Nagoya v4; cannnot be read from line number in mask file as each line does
  // not correspond to a polygon. 
  Npolygons = 764;
  
  double total_area = 0.0;
  
  // total area subtended by the mask.
  for(jj=0; jj<Npolygons; jj++)  total_area += clockwise_vertices(jj);
  
  printf("\n\nTotal W1 area: %.3f sq. degs.", total_area);
  
  return 0;
}

// ****************************************** //
// toy TSR and SSR catalogue/weight creation.
// ****************************************** //

int quadrant_sampling(){    
  // annoying, but relatively slow, memory accretion. sub sample mocks for each run.
  load_Nagoya_v6();
  
  // Nagoya v4; cannnot be read from line number in mask file as each line does
  // not correspond to a polygon. 
  
  Npolygons = 764;
  
  double delta_min;
   
  for(ii=1; ii<307; ii++){  
      printf("\n\nSub sampling galaxies to homogeneous density, on a per polygon basis. mock %d", ii);
  
      assign_spocweights_memory();
  
      load_catSpecmask(ii);
      
      // load_catMockTSR(ii);
      
      sampleGals_byPoly();
      
      write_TSR_subsample();    
  
      // randoms_mockTSR_weights();
  }

  return 0;
}


int write_mockSpec_subsample(){
    if(ii<10)        sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_00%d_W1_GalSubsample_nonpoisson_bypolygon.dat", root_dir, ii);
    else if(ii<100)  sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_0%d_W1_GalSubsample_nonpoisson_bypolygon.dat",  root_dir, ii);
    else             sprintf(filepath, "%s/Venice/mocks_W1_v8.0_Nagoya_v4_mockSpec/mock_%d_W1_GalSubsample_nonpoisson_bypolygon.dat",   root_dir, ii);
    
    int keep;
    int gal_id;
    int poly_id;
    int retained;
    
    double         weight;
    double retention_prob = 1.0;
    double draw;
    
    output = fopen(filepath, "w");
            
    // mock sampling. 
    for(i=0; i<Npolygons; i++){
        if(spoc_counts[i] == 1){    
            gal_id = quad_gal_ids[i][0];
            
            // 'retention_prob' % sampling. 
            draw = gsl_rng_uniform(gsl_ran_r);
            
            weight = retention_prob;
            
            if(draw > 1. - retention_prob)  fprintf(output, "%d \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", idx[gal_id], rax[gal_id], decx[gal_id], zx[gal_id], weight);  
        }
        /*        
        // Poisson sampling (almost).
        if(spoc_counts[i] > 1){    
            retained = gsl_ran_poisson(gsl_ran_r, 1.);
            
            if(retained > spoc_counts[i]) retained = spoc_counts[i];
            
            // polygon must be occupied.             
            for(j=0; j<retained; j++){
                keep     = j + gsl_rng_uniform_int(gsl_ran_r, spoc_counts[i] - j);
            
                // choose a 'retained' number of gals from the spoc_counts[i] gals in the polygon.  
                gal_id   = quad_gal_ids[i][keep];
                
                weight   = retained/spoc_counts[i];
            
                fprintf(output, "%d \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", idx[gal_id], rax[gal_id], decx[gal_id], zx[gal_id], weight);  
    
                quad_gal_ids[i][keep] = quad_gal_ids[i][j];
    
                quad_gal_ids[i][j]    = gal_id;
            }
         }*/
        
        // Absolute cull
        if(spoc_counts[i]>1){    
            gal_id   = quad_gal_ids[i][gsl_rng_uniform_int(gsl_ran_r, spoc_counts[i])];
                
            weight   = 1./spoc_counts[i];
            
            fprintf(output, "%d \t %.6lf \t %.6lf \t %.6lf \t %.6lf \n", idx[gal_id], rax[gal_id], decx[gal_id], zx[gal_id], weight);  
        }
    }

    fclose(output);

    return 0;
}


int create_gridoverlay(){
    // 100 arc seconds by 60. 50x30
    // double half_ra  = 1.5/60./60.; // degrees
    // double half_dec = 54./60./60.; // degrees 
    
    double half_ra  = 0.0009; // degrees
    double half_dec; // degrees 
  
    half_dec = half_ra*40.;
  
    double ra, dec;
    
    // no restriction on order of vertices in polygon format. 
  
    int rows, cols;
  
    rows = (int) ceil((UpperRAlimit  - LowerRAlimit)/(2.*half_ra));
    
    cols = (int) ceil((UpperDecLimit - LowerDecLimit)/(2.*half_dec));
    
  
    sprintf(filepath, "%s/Venice/parent_gridoverlay_polygons.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(i=0; i<cols; i++){
        for(j=0; j<rows; j++){
            ra  = LowerRAlimit  + (2.*j + 1.)*half_ra;
            
            dec = LowerDecLimit + (2.*i + 1.)*half_dec;
        
            fprintf(output, "polygon(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f) \n",   ra-half_ra,  dec+half_dec, 
                                                                                            ra+half_ra,  dec+half_dec,
                                                                                            ra+half_ra,  dec-half_dec,
                                                                                            ra-half_ra,  dec-half_dec);
        }
    }                                                        
    
    fclose(output);
    
    return 0;
}


int load_gridoverlay(){
    sprintf(maskfile, "%s/Venice/parent_gridoverlay_polygons.dat", root_dir);
  
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


int toySpoc(){
    // division of parent catalogue into slits, 60x100 arc seconds, sampling performed on a per slit basis. 
    // if number of gals in cell is 0 or 1 sampling is f_Q
    // 2 |->1, 3->1, etc. implying a sampling of 1/2, 1/3. 

    int poly_id;
    
    create_gridoverlay();

    load_gridoverlay();
    
    for(ii=300; ii<307; ii++){   
      printf("\nmock %d", ii);
    
      assign_spocweights_memory();
  
      load_catSpecmask(ii);
      
      for(i=0; i<nx; i++){
        // count nx objects in spectroscopic mask catalogue. nx > n.
        addgal_overlapPoly_counts(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad, spoc_weights_spec);
        
        // returns polygon id in which the galaxy resides. 
        poly_id = idPoly(polymask, outmask, rax[i]*deg2rad, decx[i]*deg2rad);
        
        // list of ids of each galaxy in each polygon.  
        if(poly_id != -1)  quad_gal_ids[poly_id][spoc_counts[poly_id]] = i;
        
        // number of galaxies in polygon poly_id.
        spoc_counts[poly_id] += 1;
      }
      
      write_mockSpec_subsample();
    }
    
    return 0;
}
