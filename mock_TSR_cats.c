int CatalogueInput_mockTSR(char filepath[]){
  printf("\n\nOpening catalogue: %s", filepath);

  inputfile     = fopen(filepath, "r");

  if(inputfile == NULL){
    printf("\nError opening %s\n", filepath);
    return 1;
  }

  // Column  0: ra
  // Column  1: dec
  // Column  2: zobs
  // Column  3: sampling

  ch         = 0;
  Vipers_Num = 0;

  do{
    ch = fgetc(inputfile);
    if(ch == '\n')
      Vipers_Num += 1;
  } while (ch != EOF);

  rewind(inputfile);

  id             =  (int   *)   malloc(Vipers_Num*sizeof(*id));
  ra             =  (double *)  malloc(Vipers_Num*sizeof(*ra));
  dec            =  (double *)  malloc(Vipers_Num*sizeof(*dec));
  zobs           =  (double *)  malloc(Vipers_Num*sizeof(*zobs));

  // derived parameters.
  Acceptanceflag =  (bool  *)   malloc(Vipers_Num*sizeof(*Acceptanceflag));
  polarAngle     =  (double *)  malloc(Vipers_Num*sizeof(*polarAngle));
  rDist          =  (double *)  malloc(Vipers_Num*sizeof(*rDist));
  xCoor          =  (double *)  malloc(Vipers_Num*sizeof(*xCoor));
  yCoor          =  (double *)  malloc(Vipers_Num*sizeof(*yCoor));
  zCoor          =  (double *)  malloc(Vipers_Num*sizeof(*zCoor));
  sampling       =  (double *)  malloc(Vipers_Num*sizeof(*sampling));

  // redshift range 0.7<z<0.8, as traced by randoms. magnitude cut for known linear bias. volume limited to z=0.85
  for(j=0; j<Vipers_Num; j++){
    id[j] = j;

    fscanf(inputfile, "%d \t %le \t %le \t %le \t %le \n", &id[j], &ra[j], &dec[j], &zobs[j], &sampling[j]);
  }

  fclose(inputfile);

  printf("\nHOD 500s catalogue input successful.");

  return 0;
}

