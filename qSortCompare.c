

int FirstColumnCompare(const void *pa, const void *pb){
    // pointer to const. int      // This is a derefernce (indirection) operator. 
    const double *a                = *(const double **) pa;
    const double *b                = *(const double **) pb;

    if(a[0] == b[0]){
        if(a[1] - b[1] < 0)   return -1;
        else                  return  1;
    }
    
    else{ 
        if(a[0] - b[0] < 0)   return -1;
        
        else                  return  1;
    }
}


int SecondColumnCompare(const void *pa, const void *pb){
    // pointer to const. int      // This is a derefernce (indirection) operator. 
    const double *a                = *(const double **) pa;
    const double *b                = *(const double **) pb;

    if(a[1] == b[1]){
        if(a[2] - b[2] < 0)   return -1;
        else                  return  1;
    }
    
    else{ 
        if(a[1] - b[1] < 0)   return -1;
        else                  return  1;
    }
}


int TwoColumnCompareTest(){
    double** testArray;

    int length = 40;

    testArray = malloc(length*sizeof(double*));                         // rows

    for(j=0; j<length; j++){ 
        testArray[j]     = malloc(3*sizeof(double));
        testArray[j][0]  = (double) 1./(rand() % 7 + 1.);
        testArray[j][1]  = (double) 1./(rand() % 7 + 1.);
        testArray[j][2]  = (double) 1./(rand() % 7 + 1.);
    
        printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);
    }

    printf("\n\n");

    qsort(testArray, length, sizeof(testArray[0]), FirstColumnCompare);

    printf("\n\n");

    for(j=0; j<length; j++) printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);

    qsort(testArray, length, sizeof(testArray[0]), SecondColumnCompare);

    printf("\n\n");

    for(j=0; j<length; j++) printf("\n%.6f \t \t %.6f \t \t %.6f", testArray[j][0], testArray[j][1], testArray[j][2]);

    for(j=0; j<length; j++) free(testArray[j]);
    free(testArray);
    
    return 0;
}
