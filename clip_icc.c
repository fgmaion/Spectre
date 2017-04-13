#include "header_pk.h"
#include "max_gal.c"

int prep_CatalogueInput_500s(){
  // Maximum number of galaxies present in any mock of the collection (i.e. those for covariance estimate).
  if(strcmp(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W1") == 0)  max_gals = 61765;
  else if(strcmp(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W4") == 0)  max_gals = 30261;

  else{
    printf("\n\nClipping weights are not available.");

    exit(EXIT_FAILURE);
  };
  
  // Initialise to max. memory required to process all mocks.
  ra             =  (double *)  malloc(max_gals*sizeof(*ra));
  dec            =  (double *)  malloc(max_gals*sizeof(*dec));
  zobs           =  (double *)  malloc(max_gals*sizeof(*zobs));
  M_B            =  (double *)  malloc(max_gals*sizeof(*M_B));
  Acceptanceflag =  (bool   *)  malloc(max_gals*sizeof(*Acceptanceflag));
  rDist          =  (double *)  malloc(max_gals*sizeof(*rDist));
  xCoor          =  (double *)  malloc(max_gals*sizeof(*xCoor));
  yCoor          =  (double *)  malloc(max_gals*sizeof(*yCoor));
  zCoor          =  (double *)  malloc(max_gals*sizeof(*zCoor));
  sampling       =  (double *)  malloc(max_gals*sizeof(*sampling));
  fkp_galweight  =  (double *)  malloc(max_gals*sizeof(*fkp_galweight));
  clip_galweight =  (double *)  malloc(max_gals*sizeof(*clip_galweight));

  gal_z = &zobs[0];  // Choice of redshift from zcos, zpec, zphot, zobs.
  
  return 0;
}


int load_clippingweights(){
  int       line_no;
  char  format[200];

  // default ordering: double d0s[4] = {4., 6., 10., 1000.};
  if(data_mock_flag == 0)  sprintf(filepath, "%s/W1_Spectro_V7_4/mocks_v1.7/clip_weights/W%d/mock_%03d_z_%.1lf_%.1lf_%d.dat", root_dir, fieldFlag, loopCount, lo_zlim, hi_zlim, fft_size);
  if(data_mock_flag == 1)  sprintf(filepath, "%s/W1_Spectro_V7_4/data_v1.7/clip_weights/W%d/data_%.1lf_z_%.1lf_%d.dat",  root_dir, fieldFlag, lo_zlim, hi_zlim, fft_size);

  printf("\n\nFilepath: %s", filepath);
  
  inputfile = fopen(filepath, "r");
  
  line_count(inputfile, &line_no);

  for(j=0; j<line_no; j++){
    if(d0 ==    2)  fscanf(inputfile, "%lf \t %*lf \t %*lf \t %*lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 ==    4)  fscanf(inputfile, "%*lf \t %lf \t %*lf \t %*lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 ==    6)  fscanf(inputfile, "%*lf \t %*lf \t %lf \t %*lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 ==   10)  fscanf(inputfile, "%*lf \t %*lf \t %*lf \t %lf \t %*lf \t %*lf", &clip_galweight[j]);
    if(d0 == 1000)  fscanf(inputfile, "%*lf \t %*lf \t %*lf \t %*lf \t %lf \t %*lf", &clip_galweight[j]);
  }
  
  fclose(inputfile);
  
  // for(j=0; j<130; j++)  printf("\n%d \t %.6lf", j, clip_galweight[j]);

  return 0;
}


int clip_icc(void){
  sprintf(vipersHOD_dir, "/home/mjw/HOD_MockRun/W1_Spectro_V7_2/mocks_v1.7/W%d", fieldFlag);

  loopCount =  10;
  fft_size  = 512;
  
  prep_CatalogueInput_500s();

  load_clippingweights();
  
  return 0;
}
