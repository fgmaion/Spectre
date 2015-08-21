double splint_VIPERS_maskMono(double r){
    if(r<0.1)   return 1.;
    
    if(r>900.)  return 0.;
    
    if(r<hiRes_hihiRes_join){
        double Interim;
    
        splint(VIPERS_maskr_hihi, VIPERS_maskMono_hihi, VIPERS_maskMono2D_hihi, VIPERS_mask_lineNo_hihi, r, &Interim);
    
        // amplitude determined by powerlaw_regression at r=1. -> A parameter. 
        return Interim;
    }
    
    else if(r<loRes_highRes_join){
        double Interim;
    
        splint(VIPERS_maskr_hi, VIPERS_maskMono_hi, VIPERS_maskMono2D_hi, VIPERS_mask_lineNo_hi, r, &Interim);
    
        // amplitude determined by powerlaw_regression at r=1. -> A parameter. 
        return Interim;
    }
    
    else{
        double Interim;
    
        splint(VIPERS_maskr_lo, VIPERS_maskMono_lo, VIPERS_maskMono2D_lo, VIPERS_mask_lineNo_lo, r, &Interim);
    
        // amplitude determined by powerlaw_regression at r=1. -> A parameter. 
        return Interim;
    }
}


double splint_VIPERS_maskQuad(double r){
    if(r<0.5) return 0.;
    
    if(r>900.0) return 0.;
    
    if(r<hiRes_hihiRes_join){
        double Interim;
        
        splint(VIPERS_maskr_hihi, VIPERS_maskQuad_hihi, VIPERS_maskQuad2D_hihi, VIPERS_mask_lineNo_hihi, r, &Interim);
    
        return Interim;
    }
    
    else if(r<loRes_highRes_join){
        double Interim;
        
        splint(VIPERS_maskr_hi, VIPERS_maskQuad_hi, VIPERS_maskQuad2D_hi, VIPERS_mask_lineNo_hi, r, &Interim);
    
        return Interim;
    }
    
    else{
        double Interim;
    
        splint(VIPERS_maskr_lo, VIPERS_maskQuad_lo, VIPERS_maskQuad2D_lo, VIPERS_mask_lineNo_lo, r, &Interim);
    
        return Interim;
    }
}


double splint_VIPERS_maskHex(double r){
    if(r<0.5) return 0.;
    
    if(r>900.0) return 0.;
    
    if(r<hiRes_hihiRes_join){
        double Interim;
        
        splint(VIPERS_maskr_hihi, VIPERS_maskHex_hihi, VIPERS_maskHex2D_hihi, VIPERS_mask_lineNo_hihi, r, &Interim);
    
        return Interim;
    }
    
    else if(r<loRes_highRes_join){
        double Interim;
        
        splint(VIPERS_maskr_hi, VIPERS_maskHex_hi, VIPERS_maskHex2D_hi, VIPERS_mask_lineNo_hi, r, &Interim);
    
        return Interim;
    }
    
    else{
        double Interim;
    
        splint(VIPERS_maskr_lo, VIPERS_maskHex_lo, VIPERS_maskHex2D_lo, VIPERS_mask_lineNo_lo, r, &Interim);
    
        return Interim;
    }
}


double splint_VIPERS_maskMultipoles(double r, int transformOrder){
    switch(transformOrder){
        case 0:
            return  splint_VIPERS_maskMono(r);  
        case 2:
            return  splint_VIPERS_maskQuad(r);
        case 4:
            return  splint_VIPERS_maskHex(r);
    }
}


int print_windowCorrfn(){
    char windowxi_filepath[200];
  
    sprintf(windowxi_filepath, "%s_mixedRes_hex.dat", filepath);

    double r;
    
    output = fopen(windowxi_filepath, "w");
    
    for(j=1; j<20000; j++){  
        r = j*0.025;
    
        if(r<900.)  fprintf(output, "%e \t %e \t %e \t %e \n", r,  splint_VIPERS_maskMono(r), splint_VIPERS_maskQuad(r), splint_VIPERS_maskHex(r));
    }

    fclose(output);

    return 0;
}


int prep_VIPERS_maskMultipoles(){
    //** Mask autocorrelation, super high resolution: very large number density of pairs counting for r ~ 1. 

    //** Nagoya v4. **//
    // sprintf(filepath, "%s/Data/500s/maskmultipoles_W1_Nagoya_v4_xi_incnbar_%.1f_%.1f", root_dir, lo_zlim, hi_zlim);    
    
    //** Nagoya v5. **//
    // sprintf(filepath, "%s/W1_Spectro_V5_0/maskmultipoles_W1_Nagoya_v5_xi_incnbar_%.1f_%.1f", root_dir, lo_zlim, hi_zlim);    
    // sprintf(filepath, "%s/W1_Spectro_V5_0/maskmultipoles_W1_Nagoya_v5_Samhain_incnbar_xi_%.1f_%.1f", root_dir, lo_zlim, hi_zlim);
    
    //** Nagoya v6. **//
    // sprintf(filepath, "%s/W1_Nagoya_v6_mocks_work/randoms/paircounts/maskmultipoles_W1_Nagoya_v6_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1f_%.1f", root_dir, lo_zlim, hi_zlim);
    
    //** Nagoya v7. **//
    sprintf(filepath, "%s/W1_Spectro_V7_0/maskmultipoles_W1_Nagoya_v7_Samhain_incmock_specweight_nbar_fkpweighted_8000.00_xi_%.1lf_%.1lf", root_dir, lo_zlim, hi_zlim);

    char   loRes_filepath[200];
    char   hiRes_filepath[200];
    char hihiRes_filepath[200];

    sprintf(  loRes_filepath, "%s_loRes_hex.dat",   filepath);
    sprintf(  hiRes_filepath, "%s_hiRes_hex.dat",   filepath);
    sprintf(hihiRes_filepath, "%s_hihiRes_hex.dat", filepath);
    
    inputfile = fopen(hihiRes_filepath, "r");
    
    VIPERS_mask_lineNo_hihi = 0;
    ch                      = 0;

    do{
        ch = fgetc(inputfile);
        if(ch == '\n')  VIPERS_mask_lineNo_hihi += 1;
    } while(ch != EOF);

    rewind(inputfile);

    VIPERS_maskr_hihi        = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    
    VIPERS_maskMono_hihi     = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    VIPERS_maskQuad_hihi     = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    VIPERS_maskHex_hihi      = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    
    VIPERS_maskMono2D_hihi   = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    VIPERS_maskQuad2D_hihi   = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    VIPERS_maskHex2D_hihi    = malloc(VIPERS_mask_lineNo_hihi*sizeof(double));
    
    for(j=0; j<VIPERS_mask_lineNo_hihi; j++)  fscanf(inputfile, "%le \t %le \t %le \t %le \n", &VIPERS_maskr_hihi[j], &VIPERS_maskMono_hihi[j], &VIPERS_maskQuad_hihi[j], &VIPERS_maskHex_hihi[j]);
    
    for(j=0; j<VIPERS_mask_lineNo_hihi; j++){
        VIPERS_maskMono_hihi[j] /= pow(VIPERS_maskr_hihi[j], 3.); 
        VIPERS_maskQuad_hihi[j] /= pow(VIPERS_maskr_hihi[j], 3.); 
        VIPERS_maskHex_hihi[j]  /= pow(VIPERS_maskr_hihi[j], 3.);
    }
    
    flatSlope_amp(VIPERS_mask_lineNo_hihi, 0.2, 0.3, 1.0, VIPERS_maskr_hihi, VIPERS_maskMono_hihi, &mask_monopolenorm_hihi);

    for(j=0; j<VIPERS_mask_lineNo_hihi; j++){
        VIPERS_maskMono_hihi[j] /= mask_monopolenorm_hihi;
        VIPERS_maskQuad_hihi[j] /= mask_monopolenorm_hihi;
        VIPERS_maskHex_hihi[j]  /= mask_monopolenorm_hihi;
    }

    fclose(inputfile);
    
    spline(VIPERS_maskr_hihi, VIPERS_maskMono_hihi, VIPERS_mask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_maskMono2D_hihi);
    spline(VIPERS_maskr_hihi, VIPERS_maskQuad_hihi, VIPERS_mask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_maskQuad2D_hihi);
    spline(VIPERS_maskr_hihi, VIPERS_maskHex_hihi,  VIPERS_mask_lineNo_hihi, 1.0e31, 1.0e31, VIPERS_maskHex2D_hihi);
    
    // set the join value. 
    hiRes_hihiRes_join = VIPERS_maskr_hihi[VIPERS_mask_lineNo_hihi - 1];
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // high resolution

    inputfile = fopen(hiRes_filepath, "r");
    
    VIPERS_mask_lineNo_hi = 0;
    ch                    = 0;

    do{
        ch = fgetc(inputfile);
        if(ch == '\n')  VIPERS_mask_lineNo_hi += 1;
    } while(ch != EOF);

    rewind(inputfile);

    VIPERS_maskr_hi        = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    
    VIPERS_maskMono_hi     = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    VIPERS_maskQuad_hi     = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    VIPERS_maskHex_hi      = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    
    VIPERS_maskMono2D_hi   = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    VIPERS_maskQuad2D_hi   = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    VIPERS_maskHex2D_hi    = malloc(VIPERS_mask_lineNo_hi*sizeof(double));
    
    for(j=0; j<VIPERS_mask_lineNo_hi; j++)  fscanf(inputfile, "%le \t %le \t %le \t %le \n", &VIPERS_maskr_hi[j], &VIPERS_maskMono_hi[j], &VIPERS_maskQuad_hi[j], &VIPERS_maskHex_hi[j]);
    
    // calc. hi res amplitude factor. pre-Nagoya v4, normalisation range: 0.9<k<1.0
    // previously 0.1 to 0.2 for 0.7<z<1.1 Nagoya v4 mask. 
    // powerlaw_regression(VIPERS_mask_lineNo_hi, 1.0, 2.0, 1.0, VIPERS_maskr_hi, VIPERS_maskMono_hi, &mask_monopolenorm_hi, &dummy);
    
    for(j=0; j<VIPERS_mask_lineNo_hi; j++){
        VIPERS_maskMono_hi[j] /= pow(VIPERS_maskr_hi[j], 3.);
        VIPERS_maskQuad_hi[j] /= pow(VIPERS_maskr_hi[j], 3.);
        VIPERS_maskHex_hi[j]  /= pow(VIPERS_maskr_hi[j], 3.);
    }

    // calculate amplitude factor for hi res. 
    lowRes_amplitudeCalc(VIPERS_mask_lineNo_hihi, VIPERS_mask_lineNo_hi, 1., 1.5, VIPERS_maskr_hihi, VIPERS_maskr_hi, VIPERS_maskMono_hihi, VIPERS_maskMono_hi, &mask_monopolenorm_hi);

    for(j=0; j<VIPERS_mask_lineNo_hi; j++){
        VIPERS_maskMono_hi[j] *= mask_monopolenorm_hi;
        VIPERS_maskQuad_hi[j] *= mask_monopolenorm_hi;
        VIPERS_maskHex_hi[j]  *= mask_monopolenorm_hi;
    }

    fclose(inputfile);
    
    spline(VIPERS_maskr_hi, VIPERS_maskMono_hi, VIPERS_mask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_maskMono2D_hi);
    spline(VIPERS_maskr_hi, VIPERS_maskQuad_hi, VIPERS_mask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_maskQuad2D_hi);
    spline(VIPERS_maskr_hi, VIPERS_maskHex_hi,  VIPERS_mask_lineNo_hi, 1.0e31, 1.0e31, VIPERS_maskHex2D_hi);
    
    // set the join value. 
    loRes_highRes_join = VIPERS_maskr_hi[VIPERS_mask_lineNo_hi - 1];
    
    printf("\n\nhi  res. norm: %e",   mask_monopolenorm_hi);
    printf("\njoin: %e", loRes_highRes_join);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // lower resolution on larger scales.

    inputfile = fopen(loRes_filepath, "r");
    
    VIPERS_mask_lineNo_lo = 0;
    ch                    = 0;

    do{
        ch = fgetc(inputfile);
        if(ch == '\n')  VIPERS_mask_lineNo_lo += 1;
    } while(ch != EOF);

    rewind(inputfile);

    VIPERS_maskr_lo        = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    
    VIPERS_maskMono_lo     = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    VIPERS_maskQuad_lo     = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    VIPERS_maskHex_lo      = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    
    VIPERS_maskMono2D_lo   = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    VIPERS_maskQuad2D_lo   = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    VIPERS_maskHex2D_lo    = malloc(VIPERS_mask_lineNo_lo*sizeof(double));
    
    for(j=0; j<VIPERS_mask_lineNo_lo; j++)  fscanf(inputfile, "%le \t %le \t %le \t %le \n", &VIPERS_maskr_lo[j], &VIPERS_maskMono_lo[j], &VIPERS_maskQuad_lo[j], &VIPERS_maskHex_lo[j]);
    
    for(j=0; j<VIPERS_mask_lineNo_lo; j++){
        VIPERS_maskMono_lo[j] /= pow(VIPERS_maskr_lo[j], 3.);
        VIPERS_maskQuad_lo[j] /= pow(VIPERS_maskr_lo[j], 3.);
        VIPERS_maskHex_lo[j]  /= pow(VIPERS_maskr_lo[j], 3.);
    }

    fclose(inputfile);
    
    // calculate amplitude factor for lo res. 
    lowRes_amplitudeCalc(VIPERS_mask_lineNo_hi, VIPERS_mask_lineNo_lo, 8., 9., VIPERS_maskr_hi, VIPERS_maskr_lo, VIPERS_maskMono_hi, VIPERS_maskMono_lo, &mask_monopolenorm_lo);
    
    printf("\nlow res. norm: %e", mask_monopolenorm_lo);
    
    for(j=0; j<VIPERS_mask_lineNo_lo; j++){
        VIPERS_maskMono_lo[j] *= mask_monopolenorm_lo;
        VIPERS_maskQuad_lo[j] *= mask_monopolenorm_lo;
        VIPERS_maskHex_lo[j]  *= mask_monopolenorm_lo;
    }
    
    spline(VIPERS_maskr_lo, VIPERS_maskMono_lo, VIPERS_mask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_maskMono2D_lo);
    spline(VIPERS_maskr_lo, VIPERS_maskQuad_lo, VIPERS_mask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_maskQuad2D_lo);
    spline(VIPERS_maskr_lo, VIPERS_maskHex_lo,  VIPERS_mask_lineNo_lo, 1.0e31, 1.0e31, VIPERS_maskHex2D_lo);
    
    // for integral constraint.
    // prepVIPERS_kSpaceMultipole();
    
    print_windowCorrfn();
    
    return 0;
}


double splint_VIPERS_kSpaceMono(double k){
    if(k<0.001)   return 1.;
    
    if(k>1.000)   return 0.;
    
    else{
        double Interim;
    
        splint(VIPERS_k, VIPERS_kMono, VIPERS_kMono2D,  VIPERS_kSpace_multipoles_lineNo, k, &Interim);
    
        return Interim;
    }
}


double splint_VIPERS_kSpaceQuad(double k){
    if(k<0.001)   return 1.;
    
    if(r>1.000)   return 0.;
    
    else{
        double Interim;
    
        splint(VIPERS_k, VIPERS_kQuad, VIPERS_kQuad2D,  VIPERS_kSpace_multipoles_lineNo, k, &Interim);
    
        return Interim;
    }
}


int prepVIPERS_kSpaceMultipole(){
    sprintf(filepath, "%s/Data/500s/500s_kmaskMultipoles.dat", root_dir);

    inputfile = fopen(filepath, "r");
    
    VIPERS_kSpace_multipoles_lineNo = 0;
    ch                              = 0;

    do{
        ch = fgetc(inputfile);
        if(ch == '\n')  VIPERS_kSpace_multipoles_lineNo += 1;
    } while(ch != EOF);

    rewind(inputfile);

    VIPERS_k               = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
    
    VIPERS_kMono           = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
    VIPERS_kMono2D         = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
    
    VIPERS_kQuad           = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
    VIPERS_kQuad2D         = malloc(VIPERS_kSpace_multipoles_lineNo*sizeof(double));
    
    for(j=0; j<VIPERS_kSpace_multipoles_lineNo; j++)  fscanf(inputfile, "%le \t %le \t %le \n", &VIPERS_k[j], &VIPERS_kMono[j], &VIPERS_kQuad[j]);
    
    fclose(inputfile);

    spline(VIPERS_k, VIPERS_kMono, VIPERS_kSpace_multipoles_lineNo, 1.0e31, 1.0e31, VIPERS_kMono2D);
    spline(VIPERS_k, VIPERS_kQuad, VIPERS_kSpace_multipoles_lineNo, 1.0e31, 1.0e31, VIPERS_kQuad2D);
    
    return 0;
}


int lowRes_amplitudeCalc(int Nhi, int Nlo, double rmin, double rmax, double rhi[], double rlo[], double Dhi[], double Dlo[], double* norm){
    // Least squares fit for amplitude factor of low resolution, between lo and high res in the range rmin to rmax.
    double sum_DiFi       = 0.0;
    double sum_FiFi       = 0.0;
        
    int index_hi, index_lo, ceiling;
    
    for(j=0; j<Nlo; j++){  
        if(rlo[j]>rmin){
            index_lo = j;
            
            break;
        } 
    }
    
    // printf("\n!! %e \t %e", rmin, rlo[index_lo]);
    
    for(j=0; j<Nhi; j++){  
        if(rhi[j]>rmin){
            index_hi = j;
            
            break;
        } 
    }
    
    // printf("\n!! %e \t %e", rmin, rhi[index_lo]); 
        
    for(j=0; j<Nhi; j++){  
        if(rhi[j]>rmax){
            ceiling = j - 1;
            
            break;
        } 
    }
        
    for(j=0; j<Nhi; j++){    
        if((j+index_hi) < ceiling){
            // printf("\n%e \t %e \t %e \t %e", rhi[index_hi + j], rlo[index_lo + j], Dhi[index_hi + j], Dlo[index_lo + j]);   
            
            sum_DiFi     += Dlo[index_lo + j]*Dhi[index_hi + j];
            
            sum_FiFi     += Dlo[index_lo + j]*Dlo[index_lo + j];       
        }
    }
        
    printf("\n\n%e \t %e", sum_DiFi, sum_FiFi);   
        
    *norm = sum_DiFi/sum_FiFi;

    return 0;
}


int flatSlope_amp(int N, double rmin, double rmax, double sign, double ri[], double Di[], double* norm){
    int goodN             =   0;
    double sum_Di         = 0.0;

    for(i=0; i<N; i++){                
        if((ri[i]>= rmin) && (ri[i] <= rmax) && (sign*Di[i] > 0.0)){
            goodN           +=     1;
        
            sum_Di          += Di[i];
        }
    }

    // amp for flat line through points is just the mean. 
    *norm = sum_Di/goodN;

    return 0;
}


int powerlaw_regression(int N, double rmin, double rmax, double sign, double ri[], double Di[], double* norm, double* index){
    // W    = A r^(n + 3.).
    // ln W = A' + n*log(r).
    
    // A'   = ln(A)
    
    double detA           = 0.0;
    
    // parameters.
    double Adash          = 0.0;
    double A              = 0.0;
    double n              = 0.0;
    
    double sum_lnri       = 0.0;
    double sum_lnri2      = 0.0;
    double sum_lnDi       = 0.0;
    double sum_lnDi_lnri  = 0.0;  
    
    double b1             = 0.0;
    double b2             = 0.0;
    
    int goodN             =   0;

    for(i=0; i<N; i++){                
        if((ri[i]>= rmin) && (ri[i] <= rmax) && (sign*Di[i] > 0.0)){
            goodN           += 1;
            
            // printf("\n%e \t %e", ri[i], sign*Di[i]);
            
            sum_lnri        += log(ri[i]);
            sum_lnri2       += log(ri[i])*log(ri[i]);
        
            sum_lnDi        += log(sign*Di[i]);
            sum_lnDi_lnri   += log(sign*Di[i])*log(ri[i]);
        }
    }
    
    // printf("\n\n%d good N", goodN);
    
    // For a matrix A, paramater vector (A', n)^T, and vector B
    // Required to invert AP  = B. 2x2 matrix inversion.
         	
    // A reads as (a b) for a = N, b = sum log(r_i), c = sum log(r_i), d = sum log(r_i)*log(r_i).
    //            (c d)        

    // and B    = (b1, b2)^T for b1 = sum[log(Di) - 3log(ri)], b2 = sum[ln(Di)*ln(ri) - 3ln(ri)^2]
     
    // det = ad - bc.
    detA   = goodN*sum_lnri2 - sum_lnri*sum_lnri;
      
    // (A', n)^T = (A^-1)B  = (1./detA)*( d*b1 - b*b2) 
    //                                  (-c*b1 + a*b2)
    
    b1    = sum_lnDi      - 3.*sum_lnri;
    b2    = sum_lnDi_lnri - 3.*sum_lnri2;
     
    Adash = pow(detA, -1.)*(sum_lnri2*b1 - sum_lnri*b2); 
    
    A     = exp(Adash);
    
    n     = pow(detA, -1.)*(-sum_lnri*b1 + goodN*b2); 
    
    printf("\n\nPower law regression: %e \t %e", sign*A, n);
    
    *norm = sign*A;

    *index = n;

    return 0;
}
