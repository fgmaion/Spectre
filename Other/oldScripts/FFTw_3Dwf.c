int wfPkCalc(){    
    double*       W2_vecr;
    double*       W2_veck;
    
    double* FFTW2_vecr_re;
    double* FFTW2_vecr_im;

    double*        Corrfn;

    W2_veck            = malloc(n0*n1*n2*sizeof(double));
    W2_vecr            = malloc(n0*n1*n2*sizeof(double));
    
    FFTW2_vecr_re      = malloc(n0*n1*n2*sizeof(double));
    FFTW2_vecr_im      = malloc(n0*n1*n2*sizeof(double));
    
    PkCube             = malloc(n0*n1*n2*sizeof(*PkCube));
    Corrfn             = malloc(n0*n1*n2*sizeof(*Corrfn));
    
    in                 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n0*n1*n2);
    out                = (fftw_complex*) fftw_malloc(n0*n1*n2*sizeof(fftw_complex)); 
    
    p                  = fftw_plan_dft_3d(n0, n1, n2,  in, out, FFTW_FORWARD,  FFTW_ESTIMATE);
    iplan              = fftw_plan_dft_3d(n0, n1, n2, out,  in, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    printf("\n\nBeginning window function calculation.");
    
    // The true density field is multiplied by a mask, surveyMask [W(x)].
    for(j=0; j<n0*n1*n2; j++) in[j][0] = surveyMask[j];
    for(j=0; j<n0*n1*n2; j++) in[j][1] = 0.0;
    
    printf("\nPerforming FFT.");
    
    fftw_execute(p);
    
    printf("\nWindow function FFT complete.");
    
    // window. 
    for(j=0; j<n0*n1*n2; j++) W2_veck[j] = pow(n0*n1*n2, -2.0)*(pow(out[j][0], 2.) + pow(out[j][1], 2.));
    
    // 3D FFT of W^2(vec k), power spectrum window.
    for(j=0; j<n0*n1*n2; j++)  out[j][0]  = W2_veck[j];
    for(j=0; j<n0*n1*n2; j++)  out[j][1]  = 0.0;
    
    fftw_execute(iplan);
    
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_re[j] =  in[j][0]/in[0][0];
    for(j=0; j<n0*n1*n2; j++) FFTW2_vecr_im[j] =                0.0;  // correlation fn. is purely real, in[j][1] = 0.0 for all j. 
    
    for(j=0; j<n0*n1*n2; j++) W2_vecr[j]       = pow(FFTW2_vecr_re[j], 2.) + pow(FFTW2_vecr_im[j], 2.);
    
    int iii;
    
    for(iii=0; iii<1; iii++){
        fsigma8       =   0.6;  
    
        bsigma8       =  1.15;
    
        velDispersion =    3.;
    
        formPkCube();
        
        // in -> out
        for(j=0; j<n0*n1*n2; j++) out[j][0] = PkCube[j];
        for(j=0; j<n0*n1*n2; j++) out[j][1] =       0.0;
        
        printf("\nPerforming FFT.");
    
        // plan -> iplan
        fftw_execute(iplan);
    
        for(j=0; j<n0*n1*n2; j++) Corrfn[j] = in[j][0];
        
        // Carry out the convolution by inverse FFT of the correlation fn. x FFT of W^2(k).
        for(j=0; j<n0*n1*n2; j++){ 
            // Cannot use kMonopole values due to aliased measurement. 
            in[j][0]        =  Corrfn[j]*FFTW2_vecr_re[j];
            in[j][1]        =                         0.0; // Corrfn[j]*FFTW2_vecr_im[j];
        }
    
        fftw_execute(p);
        
        
        printf("\nWindow fn. convolution complete.");

        prep_pkRegression(-2., log10(modkMax), kbin_no);
        
        polarPk_modeCount = 0;
    
        for(k=0; k<n0; k++){
            for(j=0; j<n1; j++){
                for(i=0; i<n2; i++){
                    k_x = kIntervalx*i;
                    k_y = kIntervaly*j;
                    k_z = kIntervalz*k;

                    if(k_x>xNyquistWaveNumber)  k_x   -= n2*kIntervalx;
                    if(k_y>yNyquistWaveNumber)  k_y   -= n1*kIntervaly;
                    if(k_z>zNyquistWaveNumber)  k_z   -= n0*kIntervalz;

                    kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                    kmodulus                           = pow(kSq, 0.5);
                
                    mu                                 = k_z/kmodulus;
                    if(kmodulus < 0.000001)       mu   = 0.0; 
                
                    Index                              = k*n1*n2 + j*n2 + i;
                
                    // Binning calculation following inverse transform. 
                    polar_pk[polarPk_modeCount][0]     = kmodulus;
		            polar_pk[polarPk_modeCount][1]     = fabs(mu);
		        
		            // Convolved P(k), amplitude corrected for the effects of the convolution. 
		            polar_pk[polarPk_modeCount][2]     = pow(n0*n1*n2, -1.)*out[Index][0];
		            
		            polarPk_modeCount                 += 1;
                }
            }
        }
        
        sprintf(filepath, "%s/Data/maskedRSD_draftwork/kaiserLorentz_Convergence_convolvedpk_256.dat", root_dir);

        MultipoleCalc(kbin_no, mean_modk, kMonopole, kQuadrupole, polar_pk, polarPk_modeCount, filepath, 0.0, 1.0, 1);
    }
    
    return 0;
}
