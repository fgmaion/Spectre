int randGenerate(){ 
  // Selection variables. 
  MinChi3            = pow(interp_comovingDistance(redshiftLowLimit), 3.);                        // h^-1 Mpc, see fig. 14, Guzzo et al.  2013
  MaxChi3            = pow(interp_comovingDistance(redshiftHiLimit),  3.);                        // Redshift limited sample.

  IntervalChi3       = MaxChi3 - MinChi3;
  NuRandoms          = 150000000;

  rand_chi           = malloc(NuRandoms*sizeof(*rand_chi));
  rand_dec           = malloc(NuRandoms*sizeof(*rand_dec));
  rand_rA            = malloc(NuRandoms*sizeof(*rand_rA));

  Max_SinDec         =  sin(-4.18*pi/180.0);
  Min_SinDec         =  sin(-5.40*pi/180.0);  

  Max_rA             =    38.6*pi/180.0;  
  Min_rA             =    30.0*pi/180.0;  

  lineNum            =  pointRandGen(lineNum, Min_SinDec, Max_SinDec, Min_rA, Max_rA);
  printf("\nRandom generation complete.");
  
return 0;
}


int pointRandGen(int lineNum, float Min_SinDec, float Max_SinDec, float Min_rA, float Max_rA){
    float Interval_SinDec =  Max_SinDec - Min_SinDec;
    float Interval_rA     =  Max_rA     - Min_rA;
    
    // dVol = r^2 sin(theta) dr d(theta)d(phi)                                                                                                             
    //      = d(r^3) d(cos theta) d(phi)                                                                                                                                              
    //      = d(r^3) d(cos pi/2 - dec) d(ra)                                                                    
    //      = d(r^3) d(sin dec) d(ra)                                                                                                                        

    // Uniform random deviate in r^3, sin dec and right ascension.                                                                                        

    for(i=0;i<NuRandoms; i++){
        rand_chi[lineNum] = pow(MinChi3 + IntervalChi3*randGen(&randCall_chi3), 1.0/3.0);            // previously a uniform deviate in r^3, h^-1 Mpc.                              
        rand_dec[lineNum] = (180.0/pi)*asin(Min_SinDec + Interval_SinDec*randGen(&randCall_sinDec)); // previously a uniform deviate in sin dec, degrees.               
        rand_rA[lineNum]  = (180.0/pi)*(Min_rA + Interval_rA*randGen(&randCall_rA));                 // uniform deviate in right ascension, degrees                                     
	    lineNum          += 1;
    }
    
    return lineNum;
}
