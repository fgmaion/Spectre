int randNGP(){
  for(j=0; j<n0*n1*n2; j++) booldensity[j]  = 0.0;

  sprintf(filepath, "%s/Data/BoolDensity/NGP_HODMocksV3_ByCellNumber.dat", root_dir);        
                                                            
  inputfile = fopen(filepath, "rb");                                                                                                                                  
  
  if(inputfile == NULL){                                                                                                                                             
    printf("\nCalculating Randoms NGP.");

    randGenerate();
    
    RandCoorCalc();
    
    redshiftSpaceRotation(34.5, -4.79, rand_x, rand_y, rand_z, NuRandoms, 180.0 + 94.79);
    printf("\nRotation of -z axis to VIPERS los for randoms complete.");

    // Jenkins run. 
    JenkinsFold(rand_x, NuRandoms, 0);
    JenkinsFold(rand_y, NuRandoms, 1);
    JenkinsFold(rand_z, NuRandoms, 2);

    for(j=0; j<NuRandoms; j++){
        if((rand_rshift[j] > redshiftLowLimit) && (rand_rshift[j] < redshiftHiLimit)){
            // Volume limited sample between redshift 0.7 and 0.9 within the boundaries of the survey.                                                                                                             
            xlabel                      = (int) floor((rand_x[j] - AxisLimsArray[0][0])/CellSize);
            ylabel                      = (int) floor((rand_y[j] - AxisLimsArray[0][1])/CellSize);
            zlabel                      = (int) floor((rand_z[j] - AxisLimsArray[0][2])/CellSize);

            boxlabel                    =       xlabel + n2*ylabel + n2*n1*zlabel;

            if(booldensity[boxlabel] < 0.001){
                booldensity[boxlabel]   = 1;
            }
        }
    }       
  
  output = fopen(filepath, "wb");
  fwrite(booldensity, sizeof(booldensity[0]), n0*n1*n2, output);
  fclose(output);
  
  printf("\nFreeing memory assigned for randoms.");
  freeRand(); 
  }

  else{
        printf("\nReading in randoms bool density.");
        fread(booldensity, sizeof(booldensity[0]), n0*n1*n2, inputfile);    
        fclose(inputfile);                                                             
  }

  return 0;
}


