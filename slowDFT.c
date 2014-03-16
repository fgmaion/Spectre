double slowDFT(double kx, double ky, double kz){
    double Real  = 0.0;
    double Imag  = 0.0;
    double theta = 0.0;
    double Pk    = 0.0;

    for(jj=0; jj<Vipers_Num; jj++){
       theta = kx*xCoor[jj] + ky*yCoor[jj] + kz*zCoor[jj];         
                    
       Real += cos(theta);
       Imag += sin(theta);
    }

    Real /= Vipers_Num;
    Imag /= Vipers_Num;
    
    Pk = pow(Real, 2.) + pow(Imag, 2.);
 
    return Pk;
}


int slowDFTcalc(){
    int rowNumber = 0;
    
    for(k=0; k<n0; k=k+10){
        for(j=0; j<n1; j=j+10){
            for(i=0; i<n2; i=i+10){
                k_x = kIntervalx*i;
                k_y = kIntervaly*j;
                k_z = kIntervalz*k;

                if(k_x>NyquistWaveNumber)  k_x    -= n2*kIntervalx;
                if(k_y>NyquistWaveNumber)  k_y    -= n1*kIntervaly;
                if(k_z>NyquistWaveNumber)  k_z    -= n0*kIntervalz;

                Index                              = k*n1*n2 + j*n2 + i;

                kSq                                = pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.);
                
                kmodulus                           = pow(pow(k_x, 2.) + pow(k_y, 2.) + pow(k_z, 2.), 0.5);
                
                mu                                 = k_x/kmodulus;
                if(kmodulus < 0.000001)       mu   = 0.0;      

                PkArray[rowNumber][0]              = kmodulus;
                
                PkArray[rowNumber][1]              = slowDFT(k_x, k_y, k_z);
            
                rowNumber                         += 1;
            }
        }
    }
                
    PkBinningCalc(rowNumber, PkArray);

    sprintf(filepath, "/disk1/mjw/HOD_MockRun/Data/Del2k/slowDFT_HOD_001.dat");

    output = fopen(filepath, "w");
    
    for(j=0; j<kBinNumb-1; j++) fprintf(output, "%e \t %e \t %e \t %d \t %e \t %e\n", meanKBin[j], del2[j], TotalVolume*binnedPk[j], modesPerBin[j], linearErrors[j], (*pt2Pk)(meanKBin[j]));

    fclose(output);
    
    return 0;
}   
