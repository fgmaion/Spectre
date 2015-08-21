int boxCoordinates(double xCoor[], double yCoor[], double zCoor[], int rowNumber){
    xlabel                  = (int) floor((xCoor[rowNumber] - AxisLimsArray[0][2])/xCellSize);    
    ylabel                  = (int) floor((yCoor[rowNumber] - AxisLimsArray[0][1])/yCellSize);
    zlabel                  = (int) floor((zCoor[rowNumber] - AxisLimsArray[0][0])/zCellSize);
    
    boxlabel                = (int)                         xlabel + n2*ylabel + n2*n1*zlabel;

    return boxlabel;
}


int calc_overdensity(){    
    // clean overdensity before galaxy assignment. 
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][0] = 0.0;
    for(j=0; j<n0*n1*n2; j++)  overdensity[j][1] = 0.0;
    
    // cic assign the galaxies to overdensity, weighted by sampling, fkp weights and clipping weight.
    for(j=0; j<Vipers_Num; j++){
      if(Acceptanceflag[j] == true){  
	    cic_assign(1, xCoor[j], yCoor[j],    zCoor[j], (1./sampling[j])*fkp_galweight[j]*clip_galweight[j]);
      
	    // printf("\n%.2lf \t %.2lf \t %.2lf", 1./sampling[j], fkp_galweight[j], clip_galweight[j]);
      }
    }
       
    for(j=0; j<n0*n1*n2; j++)  surveyMask[j] = 0.0;
        
    for(j=0; j<rand_number; j++){
        // assign randoms to surveyMask weighted by fkp.
        if(rand_accept[j] == true)     cic_assign(0, rand_x[j], rand_y[j], rand_z[j], rand_weight[j]);    
    }
    
    // free if not pair counting window. 
    // free(rand_x);
    // free(rand_y);
    // free(rand_z); 
    
    return 0;
}
