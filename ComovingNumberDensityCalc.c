int ComovingNumberDensityCalc(){
    prepNumberDensityCalc();
    
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<zBinNumber; i++){
	        if((redshiftSlices[i] < zUtilized[j]) && (zUtilized[j] < redshiftSlices[i+1])){
                NumberAtRedshift[i]   += 1;
                
                MeanSliceRedshift[i]  += zUtilized[j];
            
                SelectedGalNumber_nz  += 1;
            }
        }
    }

    printf("\nNumber of selected galaxies for n(z) calc:  %d", SelectedGalNumber_nz);

    for(j=1; j<zBinNumber; j++){
        ComovingVolumeAtZ[j]      = VIPERS_SolidAngle*pow(3., -1.)*(pow(interp_comovingDistance(redshiftSlices[j+1]), 3.) - pow(interp_comovingDistance(redshiftSlices[j]), 3.));
        ComovingNumberDensity[j]  = (float) NumberAtRedshift[j]/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j] /= NumberAtRedshift[j];
        }
        
        else{
            MeanSliceRedshift[j]  = 0.5*(redshiftSlices[j+1] + redshiftSlices[j]); 
        }
    }
    
    for(j=1; j<zBinNumber+1; j++) ChiSlices[j] = (float) interp_comovingDistance(MeanSliceRedshift[j]);
    
    spline(ChiSlices, ComovingNumberDensity, zBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);

    sprintf(filepath, "%s/Data/nz/HODMocks_nz_%.2f.dat", root_dir, zBinWidth);
    output = fopen(filepath, "w");
    
    for(j=1; j<zBinNumber; j++){
        fprintf(output, "%e \t %e \t %d \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j], ComovingVolumeAtZ[j], ComovingNumberDensity[j], sdltNz(MeanSliceRedshift[j]), (float) NumberAtRedshift[j]/steradians2sqdegs(VIPERS_SolidAngle));
    }

    fclose(output);

    return 0;
}


int prepNumberDensityCalc(){
  zBinNumber               =  (int)      ceil((1.4 - 0.0)/zBinWidth);

  // Numerical recipes indexing for spline. arrays begin at a[1].                                                           \
                                                                                                                               
  redshiftSlices           =  (double *)  realloc(redshiftSlices,          (zBinNumber+1)*sizeof(*redshiftSlices));
  ChiSlices                =  (float  *)  realloc(ChiSlices,               (zBinNumber+1)*sizeof(*ChiSlices));
  NumberAtRedshift         =  (int    *)  realloc(NumberAtRedshift,         zBinNumber*sizeof(*NumberAtRedshift));
  ComovingNumberDensity    =  (float  *)  realloc(ComovingNumberDensity,    zBinNumber*sizeof(*ComovingNumberDensity));
  ComovingVolumeAtZ        =  (double *)  realloc(ComovingVolumeAtZ,        zBinNumber*sizeof(*ComovingVolumeAtZ));
  MeanSliceRedshift        =  (double *)  realloc(MeanSliceRedshift,        zBinNumber*sizeof(*MeanSliceRedshift));
    
  //Second derivatives.   
  ComovingNumberDensity2d  =  (float *)  realloc(ComovingNumberDensity2d,  zBinNumber*sizeof(*ComovingNumberDensity2d));

  for(j=1; j<zBinNumber+1; j++)          redshiftSlices[j]   = 0.0 + (j-1)*zBinWidth;

  for(j=0; j<zBinNumber;   j++){         
      NumberAtRedshift[j]       = 0;
      ComovingNumberDensity[j]  = 0.0;  
      ComovingVolumeAtZ[j]      = 0.0;
      MeanSliceRedshift[j]      = 0.0;
  }
  
  SelectedGalNumber_nz          = 0;

  return 0;
}


double interp_nz(double z){
    float fz = (float) z;
    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d, zBinNumber-1, fz, &fInterim);
    
    return (double) fInterim;
}


double CSR(double z){
    double  b = 17.465;
    double zt =  0.424;

    return 0.5*(1. - gsl_sf_erf(b*(zt - z)));
}


double sdltNz(double z){
    double alpha = 8.603; 
    double  beta = 1.448;
    double    z0 = 0.191;
    double     A = 3.103;
    
    // Seemingly a 40% correction needs to be applied, eqn 2. from sdlt correlation fn. paper.
    return (100./40.)*A*pow(z/z0, alpha)*exp(-1.*pow(z/z0, beta))*CSR(z);
}
