int Gaussianfilter(float array[], float xarray[], int len, float sigma){
  int newlen;
  int kernelSize;
  
  float Sum;
  float Norm = 0.0;
  float Interval;

  Interval   = xarray[10] - xarray[9];

  kernelSize = (len - 1)/8;

  newlen     = len + 2*kernelSize;

  NewArray   =  (float *)  realloc(NewArray, newlen*sizeof(*NewArray));

  for(j=0; j<newlen; j++) NewArray[j] = 0.0;

  // padded array.                                                                                                                                          
  for(j=kernelSize; j<len+kernelSize;j++)  NewArray[j] = array[j-kernelSize];

  for(i=-kernelSize; i<kernelSize+1; i++)  Norm += exp(-0.5*pow(((float) i)*Interval/sigma, 2.0));

  for(j=kernelSize; j<len+kernelSize; j++){
    Sum    = 0.0;

    for(i=-kernelSize; i<kernelSize+1; i++){
      Sum +=  pow(Norm, -1.)*exp(-0.5*pow(((float) i)*Interval/sigma, 2.0))*NewArray[j+i];
    }

    array[j-kernelSize] = Sum;
  }

  return 0;
}



double CSR(double z){
    double  b = 17.465;
    double zt =  0.424;

    return 0.5*(1. - gsl_sf_erf(b*(zt - z)));
}


double sdltNz(double z){
    // Units of sq. degs
    double alpha = 8.603; 
    double  beta = 1.448;
    double    z0 = 0.191;
    double     A = 3.103;
    
    // Seemingly a 40% correction needs to be applied, eqn 2. from sdlt correlation fn. paper.
    return (100./40.)*A*pow(z/z0, alpha)*exp(-1.*pow(z/z0, beta))*CSR(z);
}


double sdltNz_minChi2(double z){
    // Units of sq. degs.
    // double alpha =  8.386; 
    // double  beta =  1.465;
    // double    z0 =  0.202;
    // double     A =  126.9; 40% correction?
  
    double alpha =  -5.59693061;
    double beta  =  -2.3637079;
    double z0    =   0.879776911;
    double  A    =   3.93191208*pow(10., 5.);
  
    return A*pow(z/z0, alpha)*exp(-1.*pow(z/z0, beta))*CSR(z);
}


double minChi2_nz(double Chi){
    double z      = interp_inverseComovingDistance(Chi);
    
    // Assuming n(z) cnst. over the size of a cell. dr -> dz. 
    double dz     = HubbleCnst(z)*CellSize/(2.9979*pow(10., 5.));

    double dVol_z = pow(Chi, 2.)*4.*pi*CellSize;
       
    // Solid angle conversion, from one sq. deg. to VIPERS survey.
    return steradians2sqdegs(4.*pi)*sdltNz_minChi2(z)*dz/dVol_z;
}


double NChi_dChi(double Chi, double dChi){
     double z        = interp_inverseComovingDistance(Chi);
     double w1Area   = 10.837; 
     double w1w4Area = 21.47;
    
    return sdltNz_minChi2(z)*(w1w4Area/w1Area)*HubbleCnst(z)/(2.9979*pow(10., 5.))*dChi;
}


int ComovingNumberDensityCalc(){
    float TotalW1W4area = 21.47; // sq. degs.

    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    prepNumberDensityCalc();
    
    for(j=1; j<chiBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(0.4) + chiBinWidth*(j-1);

    // W1 catalogue.
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i]<interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
                NumberAtRedshift[i]   += 1;
                
                MeanSliceRedshift[i]  += zUtilized[j];
            
                SelectedGalNumber_nz  += 1;
            }
        }
    }

    // Now the W4 field. 
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CatalogueInput(filepath);

    // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
    zUtilized       =   &zcos[0];

    assignAcceptance();

    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1])){
	            NumberAtRedshift[i]   += 1;

	            MeanSliceRedshift[i]  += zUtilized[j];

	            SelectedGalNumber_nz  += 1;
	        }
        }
    }

    for(j=1; j<chiBinNumber; j++) filteredNumberAtRedshift[j] = (float) NumberAtRedshift[j];
    Gaussianfilter(filteredNumberAtRedshift, ChiSlices, chiBinNumber, nzSigma);
    
    for(j=1; j<chiBinNumber; j++) filtered_divfuncln_Atz[j] = (float) log((float) NumberAtRedshift[j]/NChi_dChi(ChiSlices[j], chiBinWidth));
    Gaussianfilter(filtered_divfuncln_Atz, ChiSlices, chiBinNumber, nzSigma);
    
    for(j=1; j<chiBinNumber; j++) filtered_divfuncln_Atz[j] = (float) NChi_dChi(ChiSlices[j], chiBinWidth)*exp(filtered_divfuncln_Atz[j]);

    for(j=1; j<chiBinNumber; j++){      // W1 + W4
        ComovingVolumeAtZ[j]            = sqdegs2steradians(TotalW1W4area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));
        ComovingNumberDensity[j]        = (float) (filteredNumberAtRedshift[j])/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j]     /= NumberAtRedshift[j];
        }
    }
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/nz/HODMocks_00%d_Nz_chiSliced_%.1f_W1andW4_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, nzSigma);
    else              sprintf(filepath, "%s/Data/nz/HODMocks_0%d_Nz_chiSliced_%.1f_W1andW4_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, nzSigma); 
    
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++) fprintf(output, "%e \t %e \t %d \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j], filteredNumberAtRedshift[j], ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth), filtered_divfuncln_Atz[j]);
    
    fclose(output);
    
    return 0;
}


int splineGaussfilteredW1_W4_nz(){
  prepNumberDensityCalc();

  if(loopCount<10)  sprintf(filepath, "%s/Data/nz/HODMocks_00%d_Nz_chiSliced_%.1f_W1andW4_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, nzSigma);
  else              sprintf(filepath, "%s/Data/nz/HODMocks_0%d_Nz_chiSliced_%.1f_W1andW4_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, nzSigma);

  inputfile     = fopen(filepath, "r");

  for(j=1; j<chiBinNumber; j++){
    // assignment suppression character, %*e.                                                                              
    fscanf(inputfile, "%le \t %e \t %d \t %e \t %e \t %*e \t %*e \n", &MeanSliceRedshift[j], &ChiSlices[j], &NumberAtRedshift[j], &filteredNumberAtRedshift[j], &ComovingNumberDensity[j]);
  
    printf("\n %le \t %e \t %d \t %e \t %e", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j], filteredNumberAtRedshift[j], ComovingNumberDensity[j]);
  }
  
  fclose(inputfile);

  spline(ChiSlices, ComovingNumberDensity, chiBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);

  return 0;
}


int prepNumberDensityCalc(){
  chiBinNumber             =  (int)       ceil((interp_comovingDistance(1.2) - interp_comovingDistance(0.4))/chiBinWidth);   

  for(j=0; j<8; j++){
    if(chiBinNumber%8 != 0){
      chiBinNumber    += 1;
    } 
  
    else{
      break;
    }
  }

  printf("\n\nChi bin width:  %f", chiBinWidth);
  printf("\nChi bin number: %d", chiBinNumber);

  // Numerical recipes indexing for spline. arrays begin at a[1].                                                      
                                                                                                                             
  redshiftSlices           =  (double *)  realloc(redshiftSlices,          (chiBinNumber+1)*sizeof(*redshiftSlices));
  ChiSlices                =  (float  *)  realloc(ChiSlices,               (chiBinNumber+1)*sizeof(*ChiSlices));
  NumberAtRedshift         =  (int    *)  realloc(NumberAtRedshift,         chiBinNumber*sizeof(*NumberAtRedshift));
  ComovingNumberDensity    =  (float  *)  realloc(ComovingNumberDensity,    chiBinNumber*sizeof(*ComovingNumberDensity));
  ComovingVolumeAtZ        =  (double *)  realloc(ComovingVolumeAtZ,        chiBinNumber*sizeof(*ComovingVolumeAtZ));
  MeanSliceRedshift        =  (double *)  realloc(MeanSliceRedshift,        chiBinNumber*sizeof(*MeanSliceRedshift));

  filteredNumberAtRedshift =  (float  *)  realloc(filteredNumberAtRedshift, chiBinNumber*sizeof(*filteredNumberAtRedshift));
  filtered_divfuncln_Atz   =  (float  *)  realloc(filtered_divfuncln_Atz,   chiBinNumber*sizeof(*filtered_divfuncln_Atz));

  //Second derivatives.   
  ComovingNumberDensity2d  =  (float *)   realloc(ComovingNumberDensity2d,  chiBinNumber*sizeof(*ComovingNumberDensity2d));

  for(j=0; j<chiBinNumber;   j++){         
      NumberAtRedshift[j]         = 0;
      ComovingNumberDensity[j]    = 0.0;  
      ComovingVolumeAtZ[j]        = 0.0;
      MeanSliceRedshift[j]        = 0.0;
      filteredNumberAtRedshift[j] = 0.0; 
      filtered_divfuncln_Atz[j]   = 0.0;
  }
  
  SelectedGalNumber_nz            = 0;

  return 0;
}



double interp_nz(double Chi){
    float fchi = (float) Chi;    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d,  chiBinNumber-1, fchi, &fInterim);
    
    return (double) fInterim;
}
