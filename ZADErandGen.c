int randGenerate(){
  char         pointing[192][10];
  float        alpha[192], dec[192]; 
  int          Q1[192], Q2[192], Q3[192], Q4[192];
  float        PointingArea =  (7.0 + 7.0 + 2.0)*(2.4 + 8.0 + 8.0)/(60.0*60.0);   // Pointing area in sq. degrees, fig 6. Guzzo et al 2013. 
  
  IntervalChi3       = MaxChi3 - MinChi3;

  sprintf(filepath, "/disk1/mjw/VIPERS/pointing_locations.txt");
  inputfile = fopen(filepath, "r");  

  if(inputfile == NULL){
    printf("Error opening %s\n", filepath);  
    return 1;
  }

  // W1 Pointings, of which there are 192 in the W1 field. 
  // Ptg        Alpha        Delta  Q1  Q2  Q3  Q4
  //                                 1   0  -1   1
  //  1 complete
  // -1 failed
  //  0 yet to be completed

  for(j=0;j<192; j++){ 
    fscanf(inputfile, "%s \t %f \t %f \t %d \t %d \t %d \t %d \n", pointing[j], &alpha[j], &dec[j], &Q1[j], &Q2[j], &Q3[j], &Q4[j]); 

    if(Q1[j] == 1) NuQuadrants += 1;
    if(Q2[j] == 1) NuQuadrants += 1;
    if(Q3[j] == 1) NuQuadrants += 1;
    if(Q4[j] == 1) NuQuadrants += 1;
  }

  fclose(inputfile);

  NuRandoms = NuQuadrants*(RandPerPoint/4);

  rand_chi  = malloc(NuRandoms*sizeof(*rand_chi));
  rand_dec  = malloc(NuRandoms*sizeof(*rand_dec));
  rand_rA   = malloc(NuRandoms*sizeof(*rand_rA));

  for(j=0;j<192; j++){         // 192: Number of Pointings in W1, in the ** v3 data release **, see pointings_3.0.txt.
    if(Q1[j]== 1){             // Q1: Top left Quadrant of pointing. 
        Max_SinDec      =  sin((dec[j] + (1.2 + 8.0)*(1./60.0))*pi/180.0); // Figure 6, Guzzo 2013
        Min_SinDec      =  sin((dec[j] - (0.0 + 0.0)*(1./60.0))*pi/180.0);  
        Max_rA          =  (alpha[j]+(1.0+7.0)/60.0)*pi/180.0;             // fig. 6, Guzzo et al.  2013,  
        Min_rA          =  (alpha[j]-(0.0+0.0)/60.0)*pi/180.0;             // fig. 6, Guzzo et al.  2013,
    	lineNum         =  pointRandGen(lineNum, Min_SinDec, Max_SinDec, Min_rA, Max_rA);
    }

    if(Q2[j]==1){              // Q2: Top right Quadrant of pointing. 
        Max_SinDec      =  sin((dec[j] + (1.2 + 8.0)*(1./60.0))*pi/180.0);  // Figure 6, Guzzo 2013
        Min_SinDec      =  sin((dec[j] - (0.0 + 0.0)*(1./60.0))*pi/180.0);  
        Max_rA          =  (alpha[j]+(0.0+0.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,  
        Min_rA          =  (alpha[j]-(1.0+7.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,
        lineNum         =  pointRandGen(lineNum, Min_SinDec, Max_SinDec, Min_rA, Max_rA);  
    }

    if(Q3[j]==1){              // Q3: Bottom right Quadrant of pointing.
        Max_SinDec      =  sin((dec[j] + (0.0 + 0.0)*(1./60.0))*pi/180.0);  // Figure 6, Guzzo 2013
        Min_SinDec      =  sin((dec[j] - (1.2 + 8.0)*(1./60.0))*pi/180.0);  
        Max_rA          =  (alpha[j]+(0.0+0.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,  
        Min_rA          =  (alpha[j]-(1.0+7.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,
	lineNum             =  pointRandGen(lineNum, Min_SinDec, Max_SinDec, Min_rA, Max_rA);
    }

    if(Q4[j]==1){              // Q4: Bottom left Quadrant of pointing.
        Max_SinDec      =  sin((dec[j] + (0.0 + 0.0)*(1./60.0) )*pi/180.0); // Figure 6, Guzzo 2013
        Min_SinDec      =  sin((dec[j] - (1.2 + 8.0)*(1./60.0) )*pi/180.0);  
        Max_rA          =  (alpha[j]+(1.0+7.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,  
        Min_rA          =  (alpha[j]-(0.0+0.0)/60.0)*pi/180.0;              // fig. 6, Guzzo et al.  2013,
	lineNum             =  pointRandGen(lineNum, Min_SinDec, Max_SinDec, Min_rA, Max_rA);
    }
  } 
  return 0;
}


int pointRandGen(int lineNum, float Min_SinDec, float Max_SinDec, float Min_rA, float Max_rA){
    float Interval_SinDec =  Max_SinDec - Min_SinDec;
    float Interval_rA     =  Max_rA     - Min_rA;

    // dVol = r^2 sin(theta) dr d(theta) d(phi)                                                                                                                                                   
    //      = d(r^3) d(cos theta) d(phi)                                                                                                                                                          
    //      = d(r^3) d(cos pi/2 - dec) d(ra)                                                                                                                                                      
    //      = d(r^3) d(sin dec) d(ra)                                                                                                                                                             

    // Uniform random deviate in r^3, sin dec and right ascension.                                                                                                                                

    for(i=0;i<RandPerPoint/4; i++){
        rand_chi[lineNum] = pow(MinChi3 + IntervalChi3*randGen(&randCall_chi3), 1.0/3.0);            // previously a uniform deviate in r^3, h^-1 Mpc.                                                     
        rand_dec[lineNum] = (180.0/pi)*asin(Min_SinDec + Interval_SinDec*randGen(&randCall_sinDec)); // previously a uniform deviate in sin dec, degrees.                                                  
        rand_rA[lineNum]  = (180.0/pi)*(Min_rA + Interval_rA*randGen(&randCall_rA));                 // uniform deviate in right ascension, degrees                                                        
	    lineNum          += 1;
    }
    return lineNum;
}
