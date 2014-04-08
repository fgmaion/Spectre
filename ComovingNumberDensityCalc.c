double CSR(double z){
    double  b = 17.465;
    double zt =  0.424;

    return 0.5*(1. - gsl_sf_erf(b*(zt - z)));
}


double sdltNz(double z){
    // Units of sq. degs. Params as declared in Sylvain's clustering paper.
    double alpha = 8.603; 
    double  beta = 1.448;
    double    z0 = 0.191;
    double     A = 3.103;
    
    // Seemingly a 40% correction needs to be applied, eqn 2. from sdlt correlation fn. paper.
    return (100./40.)*A*pow(z/z0, alpha)*exp(-1.*pow(z/z0, beta))*CSR(z);
}


double sdltNz_minChi2(double z){
    // Units of sq. degs. Params as best fit to N(z).
  
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


float NChi_dChi(float Chi, float dChi){
    double doubleChi = (double) Chi;
    double z         = interp_inverseComovingDistance(doubleChi);
    
    // Ensure normalisation over Chi Slices is one, to then multipy by Observed number of galaxies norm.
    // Any previous normalisation, e.g. TotalW1W4area, is irrelevant.
    
    return (float) TotalObservedGalaxies*sdltNz_minChi2(z)*(TotalW1W4area/W1area)*HubbleCnst(z)/(2.9979*pow(10., 5.))*dChi/(2.263029e+05);
}


float identitytransform(float xVal, int  j){
    return xVal;
}


float ln_factor_fittedNz(float xVal, int m){
    return (float) log(xVal/NChi_dChi(ChiSlices[m], chiBinWidth));
}


float factor_fittedNz_exp(float xVal, int m){
    return NChi_dChi(ChiSlices[m], chiBinWidth)*exp(xVal);
}


int Gaussianfilter(float array[], float xarray[], int len, float sigma){
    int   newlen;
    int   kernelSize;
  
    float Sum;
    float Norm = 0.0;
    float Interval;

    Interval   = xarray[10] - xarray[9];

    kernelSize = (len - 1)/8;

    newlen     = len + 2*kernelSize;

    NewArray   =  (float *)  realloc(NewArray, newlen*sizeof(*NewArray));

    for(j=0; j<newlen; j++) NewArray[j] = 0.0;

    // padded array.                                                                                                                                          
    for(j=kernelSize + 1; j<len+kernelSize;j++)  NewArray[j] = array[j-kernelSize];

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


int NormedGaussianfilter(float (*transformArray)(float, int), float (*inversetransform)(float, int), float array[], float xarray[], int len, float sigma, double externNorm, float output[]){

    float Interim[len];
  
    double InterimNorm = 0.0;
    
    for(j=1; j<len; j++)  Interim[j]   = (*transformArray)(array[j], j);
    
    Gaussianfilter(Interim, xarray, len, sigma);
    
    for(j=1; j<len; j++)  Interim[j]   = (*inversetransform)(Interim[j], j);
      
    for(j=1; j<len; j++)  InterimNorm += (double) Interim[j];
    
    for(j=1; j<len; j++)  output[j]    = (externNorm/InterimNorm)*Interim[j];

    return 0;
}


int ComovingNumberDensityCalc(){    
    TotalObservedGalaxies    =  0.0;
    dimmestAcceptedMagnitude = -99.;

    // Equal intervals in comoving distance, for both W1 and W4 fields.  
    prepNumberDensityCalc();
    
    for(j=1; j<chiBinNumber+1; j++) ChiSlices[j] = interp_comovingDistance(0.4) + chiBinWidth*(j-1);
    
    // W1 catalogue.
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i]<interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1]) && (M_B[j]<absMagCut)){
                NumberAtRedshift[i]   += 1.;
                
                MeanSliceRedshift[i]  += zUtilized[j];
                
                if(M_B[j] > dimmestAcceptedMagnitude){
                    dimmestAcceptedMagnitude = M_B[j];
                }
            }
        }
    }

    printf("\nDimmest magnitude accepted in W1:  %le", dimmestAcceptedMagnitude);
    
    // Now the W4 field. 
    if(loopCount<10)  sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_00%d_ALLINFO.cat", vipersHOD_dir, loopCount);
    else              sprintf(filepath, "%s/mocks_W4_v1.2/mock_W4_0%d_ALLINFO.cat",  vipersHOD_dir, loopCount);

    CatalogueInput(filepath);
    
    // Choice of redshift from zcos, zpec, zphot, zobs.                                                                                   
    zUtilized       =   &zcos[0];

    assignAcceptance();
    
    // Reinitialise.
    dimmestAcceptedMagnitude =  -99.0;
    
    for(j=0; j<Vipers_Num; j++){
        for(i=1; i<chiBinNumber; i++){
	        if((ChiSlices[i] < interp_comovingDistance(zUtilized[j])) && (interp_comovingDistance(zUtilized[j]) < ChiSlices[i+1]) && (M_B[j]<absMagCut)){
	            NumberAtRedshift[i]   += 1.;

	            MeanSliceRedshift[i]  += zUtilized[j];
	            
	            if(M_B[j] > dimmestAcceptedMagnitude){
                    dimmestAcceptedMagnitude = M_B[j];
                }
	        }
        }
    }
    
    printf("\nDimmest magnitude accepted in W4:  %le", dimmestAcceptedMagnitude);

    for(j=1; j<chiBinNumber; j++)  TotalObservedGalaxies       += (double) NumberAtRedshift[j];
    
    NormedGaussianfilter(&identitytransform,  &identitytransform,   NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies, filteredNumberAtRedshift);
    
    NormedGaussianfilter(&ln_factor_fittedNz, &factor_fittedNz_exp, NumberAtRedshift, ChiSlices, chiBinNumber, nzSigma, TotalObservedGalaxies, filtered_divfuncln_Atz);

    for(j=1; j<chiBinNumber; j++){      // W1 + W4
        ComovingVolumeAtZ[j]            = sqdegs2steradians(TotalW1W4area)*pow(3., -1.)*(pow(ChiSlices[j+1], 3.) - pow(ChiSlices[j], 3.));
        
        ComovingNumberDensity[j]        = (float) (filteredNumberAtRedshift[j])/ComovingVolumeAtZ[j];
        
        if(NumberAtRedshift[j] > 0){
            MeanSliceRedshift[j]       /= NumberAtRedshift[j];
        }
    }
    
    if(loopCount<10)  sprintf(filepath, "%s/Data/nz/HODMocks_00%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma);
    else              sprintf(filepath, "%s/Data/nz/HODMocks_0%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma); 
    
    output = fopen(filepath, "w");

    for(j=1; j<chiBinNumber; j++) fprintf(output, "%e \t %e \t %e \t %e \t %e \t %e \t %e\n", MeanSliceRedshift[j], ChiSlices[j], NumberAtRedshift[j]/TotalW1W4area, filteredNumberAtRedshift[j]/TotalW1W4area, ComovingNumberDensity[j], NChi_dChi(ChiSlices[j], chiBinWidth)/TotalW1W4area, filtered_divfuncln_Atz[j]/TotalW1W4area);
    
    fclose(output);
    
    return 0;
}


int splineGaussfilteredW1_W4_nz(){
  prepNumberDensityCalc();

  if(loopCount<10)  sprintf(filepath, "%s/Data/nz/HODMocks_00%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma);
  else              sprintf(filepath, "%s/Data/nz/HODMocks_0%d_SingleMockNz_chiSliced_%.1f_W1andW4_VolLim_%.2f_Gaussfiltered_%.1f.dat", root_dir, loopCount, chiBinWidth, absMagCut, nzSigma);

  printf("\nReading n(z) file for NGP calc.\n%s", filepath);

  inputfile     = fopen(filepath, "r");

  for(j=1; j<chiBinNumber; j++){
    // assignment suppression character, %*e.                                                                              
    fscanf(inputfile, "%le \t %e \t %e \t %e \t %e \t %*e \t %*e \n", &MeanSliceRedshift[j], &ChiSlices[j], &NumberAtRedshift[j], &filteredNumberAtRedshift[j], &ComovingNumberDensity[j]);
  }
  
  fclose(inputfile);

  spline(ChiSlices, ComovingNumberDensity, chiBinNumber-1, 1.0e31, 1.0e31, ComovingNumberDensity2d);

  return 0;
}


double interp_nz(double Chi){
    float fchi = (float) Chi;    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d,  chiBinNumber-1, fchi, &fInterim);
    
    return (double) fInterim;
}


int splineStefano_nz(){
    ChiSlices                =  (float  *)  realloc(ChiSlices,                101*sizeof(*ChiSlices));    
    ComovingNumberDensity    =  (float  *)  realloc(ComovingNumberDensity,    101*sizeof(*ComovingNumberDensity));
    ComovingNumberDensity2d  =  (float  *)  realloc(ComovingNumberDensity2d,  101*sizeof(*ComovingNumberDensity2d));
    
    sprintf(filepath, "%s/Stefano/nz.txt", root_dir);
    
    inputfile = fopen(filepath, "r");
    
    for(j=0; j<101; j++){
        ChiSlices[j]               = 0.0;
        ComovingNumberDensity[j]   = 0.0;
        ComovingNumberDensity2d[j] = 0.0;
    }  
    
    for(j=1; j<101; j++)  fscanf(inputfile, "%f \t %*f \t %*f \t %f \n", &ChiSlices[j], &ComovingNumberDensity[j]);
    
    fclose(inputfile);
    
    // for(j=1; j<101; j++)  printf("%f \t %f \n", ChiSlices[j], ComovingNumberDensity[j]);
    
    spline(ChiSlices, ComovingNumberDensity, 100, 1.0e31, 1.0e31, ComovingNumberDensity2d);
    
    return 0;
}


double interp_Stefano_nz(double Chi){
    float fchi = (float) Chi;    
    float fInterim;
    
    splint(ChiSlices, ComovingNumberDensity, ComovingNumberDensity2d,  100, fchi, &fInterim);
    
    return (double) fInterim;
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
  NumberAtRedshift         =  (float  *)  realloc(NumberAtRedshift,         chiBinNumber*sizeof(*NumberAtRedshift));
  ComovingNumberDensity    =  (float  *)  realloc(ComovingNumberDensity,    chiBinNumber*sizeof(*ComovingNumberDensity));
  ComovingVolumeAtZ        =  (double *)  realloc(ComovingVolumeAtZ,        chiBinNumber*sizeof(*ComovingVolumeAtZ));
  MeanSliceRedshift        =  (double *)  realloc(MeanSliceRedshift,        chiBinNumber*sizeof(*MeanSliceRedshift));

  filteredNumberAtRedshift =  (float  *)  realloc(filteredNumberAtRedshift, chiBinNumber*sizeof(*filteredNumberAtRedshift));
  filtered_divfuncln_Atz   =  (float  *)  realloc(filtered_divfuncln_Atz,   chiBinNumber*sizeof(*filtered_divfuncln_Atz));

  //Second derivatives.   
  ComovingNumberDensity2d  =  (float *)   realloc(ComovingNumberDensity2d,  chiBinNumber*sizeof(*ComovingNumberDensity2d));

  for(j=0; j<chiBinNumber;   j++){      
      NumberAtRedshift[j]         = 0.0;
      ComovingNumberDensity[j]    = 0.0;  
      ComovingVolumeAtZ[j]        = 0.0;
      MeanSliceRedshift[j]        = 0.0;
      filteredNumberAtRedshift[j] = 0.0; 
      filtered_divfuncln_Atz[j]   = 0.0;
  }

  return 0;
}
