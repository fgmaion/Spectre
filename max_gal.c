int max_array(int a[], int num_elements){
  int i, max=-32000;
  
  for(i=0; i<num_elements; i++){
      if(a[i]>max){
	max=a[i];
      }
  }

  return(max);
}


int max_gal(){
  int Ngal[154], max;

  for(j=1; j<154; j++){
    sprintf(filepath, "%s/mock_%03d_VAC_Nagoya_v6_Samhain.dat",  vipersHOD_dir, j);
    
    inputfile = fopen(filepath, "r");

    line_count(inputfile, &Ngal[j]);

    printf("\n%d \t %d", j, Ngal[j]);
    
    fclose(inputfile);
  }

  max = max_array(Ngal, 154);
  
  printf("\n\nMax gals. %d", max);
  
  return max;
}
