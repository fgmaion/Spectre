int wedgeMockCreate(double ya, double za, double yb, double zb, double yc, double zc, double yd, double zd, double xlo, double xhi){
/*
    Continous slice in x dir. 
    Wedge shape in y-z, plane. 
    
      d                  c
      \-----------------/
       \               /   
     l1 \             / l2
         \           / 
          \---------/
          a         b
    
    l1: y = m1*x + c1 etc.
    
    ^ z dir     -> y dir
    |

*/

    double m1, m2;
    double c1, c2;
    
    m1 = (zd - za)/(yd - ya);
    m2 = (zc - zb)/(yc - yb);
    
    c1 = (za*yd - ya*zd)/(yd - ya);
    c2 = (zb*yc - yb*zc)/(yc - yb);
    
    sprintf(filepath, "%s/Data/HODCube/zcube_gal_-20.0.dat", root_dir);
    
    printf("\n\nOpening catalogue: %s", filepath);
    
    inputfile     = fopen(filepath, "r");  
    if(inputfile == NULL){  
        printf("Error opening %s\n", filepath); 
        return 1;
    }

    // Column  0: x coordinate.                                                                                                         
    // Column  1: y coordinate.                                                                                         
    // Column  2: z coordinate.  
                                                                                                  
    ch         = 0;
    Vipers_Num = 0;
    
    do{
        ch = fgetc(inputfile);        
        if(ch == '\n')
       	  Vipers_Num += 1;
    } while (ch != EOF);

    printf("\n\nNumber of galaxies in catalogue:  %d", Vipers_Num);

    rewind(inputfile);
    
    double       xi[5481349];
    double       yi[5481349];
    double       zi[5481349];
    double inc_flag[5481349];
    
    for(j=0; j<Vipers_Num; j++){
        fscanf(inputfile, "%lf \t %lf \t %lf \n", &xi[j], &yi[j], &zi[j]);
        
        inc_flag[j] = 0.0;
    }

    printf("\nInput of catalogue for wedge-mock creation successful.");
    
    int realisationNumber;
    int jjj, iii;
    
    for(jjj=0; jjj<2; jjj++){
        for(iii=0; iii<20; iii++){
            realisationNumber = jjj*20 + iii;
            
            printf("\n%d", realisationNumber);
            
            sprintf(filepath, "%s/Data/HODCube/WedgeMocks/zMock_%d.dat", root_dir, realisationNumber);
    
            output = fopen(filepath, "w");
        
            for(j=0; j<Vipers_Num; j++){
                if((zi[j] > za) && (zi[j] < zd) && (zi[j] >= m1*yi[j] + c1) && (zi[j] >= m2*yi[j] + c2) && (xi[j]>xlo) && (xi[j]<xhi)){
                    fprintf(output, "%le \t %le \t %le \n", xi[j], yi[j], zi[j]);
                
                    if(inc_flag[j] > 0.0){
                        printf("\nError: Overlapping catalogs.");
                    }
                
                    else{
                        inc_flag[j] = 1.0;                        
                    }
                }
            }
    
            fclose(output);
    
            rollxyz(50., 0., 0., xi, yi, zi, Vipers_Num);
        }
     
        rollxyz(0., 500., 0., xi, yi, zi, Vipers_Num);  
    }
    
    return 0;
}
