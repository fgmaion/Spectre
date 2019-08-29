#define FOLDFACTOR 2.0

int fold_volume(){
    // Jenkin's run to beat aliasing of P(k) near the Nyquist wavenumber. 
    AxisLimsArray[1][0] /= FOLDFACTOR;                     
    AxisLimsArray[1][1] /= FOLDFACTOR;                     
    AxisLimsArray[1][2] /= FOLDFACTOR;                     
         
    return 0;
}


int fold_randoms(){
    // Jenkins run to beat aliasing. 
    JenkinsFold(rand_x, rand_number, 0);  
    JenkinsFold(rand_y, rand_number, 1);
    JenkinsFold(rand_z, rand_number, 2);
    
    printf("\n\nJenkins folded, Stefano basis, randoms co-ordinates.");
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_x, rand_number), arrayMax(rand_x, rand_number));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_y, rand_number), arrayMax(rand_y, rand_number));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", arrayMin(rand_z, rand_number), arrayMax(rand_z, rand_number));
    
    return 0;
}


int fold_galaxies(){
    // Jenkins run to beat aliasing. 
    JenkinsFold(xCoor, Vipers_Num, 2);
    JenkinsFold(yCoor, Vipers_Num, 1);
    JenkinsFold(zCoor, Vipers_Num, 0);
          
    printf("\n\nJenkins folded, Stefano basis, galaxy co-ordinates.");
    printf("\nx: %.1lf \t %.1lf h^-1 Mpc", AcceptedMin(xCoor, Acceptanceflag, Vipers_Num), AcceptedMax(xCoor, Acceptanceflag, Vipers_Num));
    printf("\ny: %.1lf \t %.1lf h^-1 Mpc", AcceptedMin(yCoor, Acceptanceflag, Vipers_Num), AcceptedMax(yCoor, Acceptanceflag, Vipers_Num));
    printf("\nz: %.1lf \t %.1lf h^-1 Mpc", AcceptedMin(zCoor, Acceptanceflag, Vipers_Num), AcceptedMax(zCoor, Acceptanceflag, Vipers_Num));
        
    return 0;
}


int JenkinsFold(double original[], int lenArray, int axis){
  for(j=0; j<lenArray; j++)  original[j] = fmod(original[j], AxisLimsArray[1][axis]/FOLDFACTOR));               

  return 0;
}


int JenkinsFoldTest(double original[], int lenArray, int axis){
    double Interim;

    for(j=0; j<10; j++){
        Interim = original[j];
    
        original[j] = fmod(original[j], (AxisLimsArray[1][axis] - AxisLimsArray[0][axis]));               
    
        printf("\n%e \t %e", Interim, original[j]);
    }

    return 0;
}
