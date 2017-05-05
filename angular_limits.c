int set_angularlimits(int dataormock, int fieldFlag){
  if(dataormock == 0){ // mock
    W1area                    =     10.692;                 // Nagoya v6 spectroscopic mask (& Samhain).
    W4area                    =      5.155;                 // Dec cut at -5.97 in the mocks;   

    if(fieldFlag == 1){
      LowerRAlimit            =     30.175;                 
      UpperRAlimit            =     38.797;
      CentreRA                =     34.492;                 

      LowerDecLimit           =     -5.970;
      UpperDecLimit           =     -4.171;
      CentreDec               =     -5.091;                
    }

    else if (fieldFlag ==4){
      LowerRAlimit            =    330.046;                 
      UpperRAlimit            =    335.389;
      CentreRA                =    332.638;

      LowerDecLimit           =      0.862;
      UpperDecLimit           =     2.3696;
      CentreDec               =      1.583;
    }

    else{
      return 1;
    }
  }

  else if(dataormock == 1){ // data
    W1area                    =      10.763;                
    W4area                    =       5.155;                

    if(fieldFlag == 1){
      LowerRAlimit            =     30.1893;
      UpperRAlimit            =     38.8022;
      CentreRA                =     34.3213;

      LowerDecLimit           =     -5.9801;
      UpperDecLimit           =     -4.1715;
      CentreDec               =     -5.1188;
    }

    else if(fieldFlag == 4){
      LowerRAlimit            =    330.0452;
      UpperRAlimit            =    335.3890;
      CentreRA                =    332.7049;

      LowerDecLimit           =      0.8621;
      UpperDecLimit           =     2.36950;
      CentreDec               =      1.5549;
    }

    else{
      return 1;
    }
  }

  else{
    return 1;
  }

  TotalW1W4area               = W1area + W4area;            // Required for <n(z)> calculation.
  
  return  0;
}
