int assignAcceptance(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = true;
    
    for(j=0; j<Vipers_Num; j++){
        if((redshiftLowLimit > zUtilized[j]) || (zUtilized[j] > redshiftHiLimit) || (M_B[j] > absMagCut)){
            Acceptanceflag[j]  = false;
        }
    }

    return 0;
}


int assignAcceptanceCube(){
    for(j=0; j<Vipers_Num; j++)                     Acceptanceflag[j]  = true;

    return 0;
}
