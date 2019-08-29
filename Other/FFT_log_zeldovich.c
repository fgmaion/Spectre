int zeldovich_memory(int FFTlogRes, double beta, double velDispersion){
    //Can bias simply by changing q, positive or negative. input func is multiplie by suitable (kr)^q below. 
    eta_para    = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), -0.0, 0 + 0.5);

    eta_perp    = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.), -2.0, 1 + 0.5);
    
    meanInfall  = FFTLog_init(FFTlogRes, pow(10., -10.), pow(10., 14.),  0.0, 1 + 0.5);    


    FFTLog_setInput(eta_para, beta, velDispersion);

    FFTLog_setInput(eta_perp, beta, velDispersion);
    
    FFTLog_setInput(meanInfall, beta, velDispersion);
    
    return 0;
}


int zeldovich_eta_para(FFTLog_config* fc, FFTLog_config* perp){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    for(i=0; i<fc->N; i++)  fc->input[i][0]    = pow(fc->krvals[i][0], -0.5)*pow(2.*pi, -3./2.)*fc->pk[i][0];
    for(i=0; i<fc->N; i++)  fc->input[i][1]    = 0.0;
   
    // power law bias may be required.                                                                                                                        
    for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][0], -fc->q);            
                                                                     
    FFTLog(fc, fc->forwardplan, fc->backwardplan);

    // and debias.                                                                                                                                           
    for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][1], -fc->q);  

    for(i=0; i<fc->N;   i++) fc->xi[i][0]      = fc->output[i][0]*pow(fc->krvals[i][1], -3./2.);
    
    // -2.*eta_perp.
    for(i=0; i<fc->N;   i++) fc->xi[i][0]     -= 2.*perp->xi[i][0];
    
    return 0;
}


int zeldovich_eta_perp(FFTLog_config* fc){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    // for(i=0; i<fc->N; i++)  fc->input[i][0] = pow(fc->krvals[i][0], 2.)*exp(-0.5*fc->krvals[i][0]*fc->krvals[i][0]); 
    for(i=0; i<fc->N; i++)  fc->input[i][0]    = fc->pk[i][0]*pow(2.*pi*fc->krvals[i][0], -3./2.);
    for(i=0; i<fc->N; i++)  fc->input[i][1]    = 0.0;
    
    // power law bias maybe required.                                                                                                                        
    for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][0], -fc->q);                                                                           
                                                                                                                                                             
    FFTLog(fc, fc->forwardplan, fc->backwardplan);

    // and debias.                                                                                                                                           
    for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][1], -fc->q);  

    // for(i=0; i<fc->N;   i++) fc->xi[i][0]   = fc->output[i][0]*pow(fc->krvals[i][1], -1.);
    for(i=0; i<fc->N;   i++) fc->xi[i][0]      = fc->output[i][0]*pow(fc->krvals[i][1], -5./2.);
    
    return 0;
}


int zeldovich_meanInfall(FFTLog_config* fc){
    // Last argument of FFTLog_init is the order of the Hankel transform. 
    double transformOrder;
    
    transformOrder  =  fc->mu - 0.5;
    
    for(i=0; i<fc->N; i++)  fc->input[i][0] = -1.*pow(fc->krvals[i][0], 0.5)*pow(2.*pi, -3./2.)*fc->pk[i][0];
    for(i=0; i<fc->N; i++)  fc->input[i][1] = 0.0;
    
    // power law bias maybe required.                                                                                                                        
    // for(i=0; i<fc->N; i++)   fc->input[i][0]  *= pow(fc->krvals[i][0], -fc->q);                                                                           
                                                                                                                                                              
    FFTLog(fc, fc->forwardplan, fc->backwardplan);

    // and debias.                                                                                                                                           
    // for(i=0; i<fc->N; i++)  fc->output[i][0]  *= pow(fc->krvals[i][1], -fc->q);  

    for(i=0; i<fc->N;   i++) fc->xi[i][0]   = fc->output[i][0]*pow(fc->krvals[i][1], -3./2.);
    
    return 0;
}


int sigma_eta(double* sigma2eta){
    *sigma2eta  = qromb(pt2Pk, 0.0001, 10.);
    
    *sigma2eta *= pow(6.*pi*pi, -1.);
    
    printf("\n\nSigma^2_eta: %e\n\n", *sigma2eta);
    
    return 0;
}



int cubature_trialfunc(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    double sigma = *((double *) fdata); // we can pass σ via fdata argument
    
    double sum   = 0.0;
    
    for(i=0; i<ndim; i++) sum += x[i]*x[i];
    
    // compute the output value: note that fdim should == 1 from below
    fval[0] = exp(-sigma*sum);
    
    return 0; // success
}


int cubature_zeldovich_GaussianIntegral(unsigned ndim, const double *q, void *fdata, unsigned fdim, double *fval){
    double detA     = 0.0;
    double modq     = 0.0;
    double x2       = 0.0;
    double exponent = 0.0;
    
    double sigma2perp;
    double sigma2para;
    
    double invA_vecx[3];
    
    for(i=0; i<3; i++) invA_vecx[i] = 0.0;
    
    for(i=0; i<ndim; i++) modq += q[i]*q[i];

    modq = sqrt(modq);

    sigma2perp = splint_sigma2perp(modq);
    sigma2para = splint_sigma2para(modq);

    // Here q[i] represents the Lagrangian position, over which the integral is taken. and vec x = vec r - vec q.   
    for(i=0; i<ndim; i++)  zel_vecx[i] = zel_vecr[i] - q[i];
    
    for(i=0; i<ndim; i++)  x2         += zel_vecx[i]*zel_vecx[i];
    
    // exponent  = -0.5*pow(sigma2perp, -1.)*x2;
    // exponent += -0.5*(sigma2perp - sigma2para)*pow(sigma2perp*sigma2para, -1.)*pow(zel_vecx[0] + zel_vecx[1] + zel_vecx[2], 2.); 
    
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            invA[i][j]           = 0.0;
            
            if(i==j) invA[i][j] += pow(sigma2perp, -1.);
            
            invA[i][j]          += (sigma2perp - sigma2para)*pow(sigma2perp*sigma2para, -1.)*pow(modq, -2.)*q[i]*q[j];
        }
    }
    
    for(i=0; i<3; i++){  
        for(j=0; j<3; j++){
            invA_vecx[i] += invA[i][j]*zel_vecx[j];
        }
    }
    
    for(i=0; i<3; i++)  exponent += zel_vecx[i]*invA_vecx[i];
    
    exponent *= -0.5;
    
    // Actually (2.pi)^(3./2.) |A|^0.5
    detA      = sigma2perp*sqrt(sigma2para)*pow(2.*pi, 3./2.);
    
    // compute the output value: note that fdim should == 1 from below
    fval[0]   = exp(exponent)/detA;

    return 0; // success
}


int calc_zeldovichxi(){
    // For an odd numbered resolution, transform and its inverse are identical except q -> -q.
    FFTlogRes        = 4097;
    
    zeldovich_memory(FFTlogRes, beta, velDispersion);
    
    
    zeldovich_eta_perp(eta_perp);
    
    zeldovich_eta_para(eta_para, eta_perp);
    
    zeldovich_meanInfall(meanInfall);
    
    
    sigma_eta(&sigma2eta);
    
        
    print_linearpk();
    
    for(j=0; j<FFTlogRes; j++){
        if(eta_perp->krvals[j][1] > 0.1){
            startint = j;
            
            break;
        }
    }
    
    for(j=0; j<FFTlogRes; j++){
        if(eta_perp->krvals[j][1] > 1000.){
            endint = j;
            
            break;
        }
    }
    
    sigma2rvals    = malloc((endint - startint)*sizeof(double));
    
    sigma2perp     = malloc((endint - startint)*sizeof(double));
    sigma2para     = malloc((endint - startint)*sizeof(double));
    
    sigma2perp_2d  = malloc((endint - startint)*sizeof(double));
    sigma2para_2d  = malloc((endint - startint)*sizeof(double));
    
    for(j=startint; j<endint; j++){
        sigma2rvals[j - startint] = eta_perp->krvals[j][1];
    
        sigma2perp[j - startint]  = 2.*(sigma2eta - eta_perp->xi[j][0]);
        sigma2para[j - startint]  = 2.*(sigma2eta - eta_para->xi[j][0]);
        
        // printf("\n%e \t %e \t %e", sigma2rvals[j - startint], sigma2perp[j - startint], sigma2para[j - startint]);
    }
    
    print_zeldovichDispersions();
    
    spline(sigma2rvals, sigma2perp, (endint-startint), 1.0e31, 1.0e31, sigma2perp_2d);
    spline(sigma2rvals, sigma2para, (endint-startint), 1.0e31, 1.0e31, sigma2para_2d);
    
    
    invA = malloc(3*sizeof(double*));
    
    for(j=0; j<3; j++) invA[j] = malloc(3*sizeof(double));
    
    for(i=0; i<3; i++){
        for(j=0; j<3; j++){
            invA[i][j] = 0.0;
        }
    }
    
    zel_vecr[0]     =   0.;
    zel_vecr[1]     =   0.;
    
    zel_scale       =  40.0;
    
    sprintf(filepath, "%s/Data/zeldovich_pk/r_r2xi.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(jj=40; jj>=0; jj--){
        zel_vecr[2]     =                   jj*5.;
    
        zel_qmin[0]     = zel_vecr[0] - zel_scale;
        zel_qmin[1]     = zel_vecr[1] - zel_scale;
        zel_qmin[2]     = zel_vecr[2] - zel_scale;
        
        zel_qmax[0]     = zel_vecr[0] + zel_scale;
        zel_qmax[1]     = zel_vecr[1] + zel_scale;
        zel_qmax[2]     = zel_vecr[2] + zel_scale;
    
        zel_r           = pow(zel_vecr[0]*zel_vecr[0] + zel_vecr[1]*zel_vecr[1] + zel_vecr[2]*zel_vecr[2], 0.5);
    
        hcubature(1, cubature_zeldovich_GaussianIntegral, &dummy, 3, zel_qmin, zel_qmax, 0, 0, 1e-4, ERROR_INDIVIDUAL, &zel_xi, &zel_err);
     
        printf("\n%.2e \t %.2e", zel_r, zel_xi - 1.);       
    
        fprintf(output, "%e \t %e\n", zel_r, pow(zel_r, 2.)*(zel_xi -1.));
    }
    
    fclose(output);
    
    return 0;
}


double splint_sigma2perp(double r){
    double Interim;
    
    if(r<0.1)  return 0.;
    
    if(r>999.){
        splint(sigma2rvals, sigma2perp, sigma2perp_2d, (endint-startint), 990., &Interim);
    
        return Interim;
    }
    
    else{
        splint(sigma2rvals, sigma2perp, sigma2perp_2d, (endint-startint), r, &Interim);
    
        return Interim;
    }
}


double splint_sigma2para(double r){
    double Interim;
    
    if(r<0.1)  return 0.;
    
    if(r>999.){
        splint(sigma2rvals, sigma2para, sigma2para_2d, (endint-startint), 990., &Interim);
    
        return Interim;
    }
    
    else{
        splint(sigma2rvals, sigma2para, sigma2para_2d, (endint-startint), r, &Interim);
    
        return Interim;
    }
}


int print_zeldovichDispersions(){
    sprintf(filepath, "%s/Data/zeldovich_pk/q_sigma2perp_sigma2para_meanInfall.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<1200; j++){  
        fprintf(output, "%e \t %e \t %e \n", 1.*j, splint_sigma2perp(1.*j), splint_sigma2para(1.*j));
    }
    
    fclose(output);
    
    return 0;
} 


int print_linearpk(){
    sprintf(filepath, "%s/Data/zeldovich_pk/k_linear_hodpk.dat", root_dir);
    
    output = fopen(filepath, "w");
    
    for(j=0; j<eta_perp->N; j++)  fprintf(output, "%e \t %e \t %e \n", eta_perp->krvals[j][0], linearPk_Gaussian(eta_perp->krvals[j][0]), HODPk_Gaussian(eta_perp->krvals[j][0]));

    fclose(output);

    return 0;
} 
