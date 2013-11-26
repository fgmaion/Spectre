int sphereCentre(){
    const gsl_rng_type* T;
    gsl_rng*            r;
    
    gsl_rng_env_setup();
    
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    
    // Generate random x, y, z co-ordinates for the sphere centre.  Use periodicity of box to roll the co-ordinates
    // such that the sphere lies in the centre of the new co-ordinate system. 
    
    xcentre       = 1000.*gsl_rng_uniform(r);
    ycentre       = 1000.*gsl_rng_uniform(r);
    zcentre       = 1000.*gsl_rng_uniform(r);
    
    printf("\n\n%e \t %e \t %e", xcentre, ycentre, zcentre);
    
    xroll         = -1.*(xcentre - 0.5*(AxisLimsArray[1][0] - AxisLimsArray[0][0]));
    yroll         = -1.*(ycentre - 0.5*(AxisLimsArray[1][1] - AxisLimsArray[0][1]));
    zroll         = -1.*(zcentre - 0.5*(AxisLimsArray[1][2] - AxisLimsArray[0][2]));
    
    gsl_rng_free(r);

    return 0;
}


int rollcube(float xCoor[], float yCoor[], float zCoor[], int galNumber){
    sphereCentre();
    
    for(j=0; j<galNumber; j++){
        xCoor[j] += xroll;
    
        if(xCoor[j] > 1000.) xCoor[j] -= 1000.;
        if(xCoor[j] <    0.) xCoor[j] += 1000.;
    
        yCoor[j] += yroll;
        
        if(yCoor[j] > 1000.) yCoor[j] -= 1000.;
        if(yCoor[j] <    0.) yCoor[j] += 1000.;
    
        zCoor[j] += zroll;
        
        if(zCoor[j] > 1000.) zCoor[j] -= 1000.;
        if(zCoor[j] <    0.) zCoor[j] += 1000.;
    
    }

    return 0;
}


double signum(double abscissa){
    if(abscissa >= 0.0){  
        return  1.0;
    }
    
    else{
        return -1.0;
    }
}

