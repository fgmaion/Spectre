int assign_randmemory(){
    rand_ra       = (double *) realloc(rand_ra,  rand_number*sizeof(*rand_ra));
    rand_dec      = (double *) realloc(rand_dec, rand_number*sizeof(*rand_dec));
    rand_chi      = (double *) realloc(rand_chi, rand_number*sizeof(*rand_chi));

    rand_x        = (double *) realloc(rand_x,   rand_number*sizeof(*rand_x));
    rand_y        = (double *) realloc(rand_y,   rand_number*sizeof(*rand_y));
    rand_z        = (double *) realloc(rand_z,   rand_number*sizeof(*rand_z));
    
    rand_weight   = (double *) realloc(rand_weight, rand_number*sizeof(*rand_weight));
    
    // rand_accept   = (bool *)   realloc(rand_accept, rand_number*sizeof(*rand_accept));
    
    return 0;
}


int randoms_maskGen(){   
    rand_number = 600000;
                  
    assign_randmemory();
    
    for(j=0; j<rand_number; j++){
         rand_ra[j] =  LowerRAlimit  + (UpperRAlimit  -  LowerRAlimit)*gsl_rng_uniform(gsl_ran_r);
         
        // uniform in sin(dec). 
        rand_dec[j] = sin(LowerDecLimit*pi/180.) + (sin(UpperDecLimit*pi/180.) - sin(LowerDecLimit*pi/180.))*gsl_rng_uniform(gsl_ran_r);
        
        rand_dec[j] = (180./pi)*asin(rand_dec[j]);
    }
    
    for(j=0; j<rand_number; j++){  
      // Index     = (int) gsl_rng_uniform_int(gsl_ran_r, redshiftsNumber);
      // rand_z[j] = redshifts[Index]; 

      Interim      = pow(loChi, 3.) + (pow(hiChi, 3.) - pow(loChi, 3.))*gsl_rng_uniform(gsl_ran_r);

      rand_chi[j]  = pow(Interim, 1./3.);
    }
    
    for(j=0; j<rand_number; j++){  
      rand_ra[j]    *= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   *= (pi/180.0);                                 // Converted to radians.
    
      rand_x[j]      = rand_chi[j]*cos(rand_dec[j])*cos(rand_ra[j]);        
      rand_y[j]      = rand_chi[j]*cos(rand_dec[j])*sin(rand_ra[j]);
      rand_z[j]      = rand_chi[j]*sin(rand_dec[j]);
      
      rand_ra[j]    /= (pi/180.0);                                 // Converted to radians.
      rand_dec[j]   /= (pi/180.0);                                 // Converted to radians.
    }
    
    // double draw, max_nz;
    
    // assumes n bar monotonically decreasing from loChi to hiChi.
    // max_nz = interp_nz(loChi); // mocks
    
    // smoothed, upweighted nbar of W1_Spectro_V7_0.
    // max_nz = interp_nz(1656.);
    
    sprintf(filepath, "%s/W1_Spectro_V7_2/randoms/randoms_W%d_parent_xyz_%.1lf_%.1lf_%d.cat", root_dir, fieldFlag, lo_zlim, hi_zlim, rand_number);

    output = fopen(filepath, "w");
    
    for(j=0; j<rand_number; j++){  
        // draw   = gsl_rng_uniform(gsl_ran_r);
    
        // nbar dependence. volume dependence (at cnst. number density) taken care of above.
        // if(draw < interp_nz(rand_chi[j])/max_nz)  
        fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \n", rand_ra[j], rand_dec[j], rand_chi[j], rand_x[j], rand_y[j], rand_z[j]);
    }
    
    fclose(output);
    
    printf("\n\nrandoms co-ordinates.");
    
    printf("\n%e \t %e", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\n%e \t %e", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\n%e \t %e", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    return 0;
}


int randoms_maskGen_GriddingError(){   
    double xx, yy, zz, xerr, yerr, zerr; 
    
    double scale;
    
    scale  = AxisLimsArray[1][0]/fft_size;
    
    
    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/500s/randoms_W1_500s_Nagoya_v4_gridded_%.2f_xyz_%.1f_%.1f.cat", scale, lo_zlim, hi_zlim);

    output = fopen(filepath, "w");
    
    for(j=0; j<rand_number; j++){  
        xerr = scale*(gsl_rng_uniform(gsl_ran_r) - 0.5);
        yerr = scale*(gsl_rng_uniform(gsl_ran_r) - 0.5);
        zerr = scale*(gsl_rng_uniform(gsl_ran_r) - 0.5);
    
        fprintf(output, "%e \t %e \t %e \n", rand_x[j] + xerr, rand_y[j] + yerr, rand_z[j] + zerr);
    }
    
    fclose(output);
    
    return 0;
}
