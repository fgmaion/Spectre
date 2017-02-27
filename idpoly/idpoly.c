#include <string.h>

#define NCHAR 40
#define NR_END 1
#define FREE_ARG char*
#define deg2rad M_PI/180

/*
double *dvector(long nl, long nh)
{
  double *v;
  
  v=(double *) calloc((size_t)(nh-nl+1+NR_END),sizeof(double));
  if (!v) printf("allocation failure in dvector()\n");
  return v-nl+NR_END;
}*/

/*
void free_dvector(double *v, long nl, long nh)
{
  free((FREE_ARG) (v+nl-NR_END));
}*/

int ct_lines(char *file)
{
  FILE *fic;
  int n;
  
  fic=fopen(file,"r");
  if (fic) {
    n=0;
    while(!feof(fic)) {
      fscanf(fic,"%*[^\n]\n");
      n++;
    }
    fclose(fic);
    return n;
  } else {
    fprintf(stderr,"Can't read %s file!\n",file);
    return 0;
  }
}

int getStrings(char *line, char *strings, char *delimit, size_t *N)
{
  int i,j,begin,length;
  
  if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
  
  i = 0;
  j = 0;
  while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
    begin = i;
    while((line[i] == *delimit || line[i] == '\t') && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
    begin = i;
    while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
    length = i - begin;
    if(length > 0){
      strncpy(strings+NCHAR*j,&line[begin],length);
      strcpy(strings+NCHAR*j+length,"\0");
      j++;
    }
  }
  
  (*N) = j;
  
  if(*N > 0) return 1;
  else return 0;
}


int polymaskCalc(){
  FILE *fic;
  
  char file[BUFSIZ],maskfile[BUFSIZ],weightfile[BUFSIZ];
  
  double outmask[2];
  NodeP *polymask=NULL;
  
  double *ra,*dec,*wei;
  
  size_t i,n,nw;
  /*
  if(argc != 3 && argc != 4) {
      printf("usage: %s <catalogue> <mask_file> <weight_file>\n",argv[0]);
      exit(1);
  }*/

  sscanf(argv[1],"%s",file);
  sscanf(argv[2],"%s",maskfile);
  if (argc==4) sscanf(argv[3],"%s",weightfile);

  fic=fopen(file,"r");
  if (fic==0) {
    fprintf(stderr,"Can't open catalogue!\n");
    exit(1);
  } else fclose(fic);
  
  fic=fopen(maskfile,"r");
  if (fic==0) {
    fprintf(stderr,"Can't open mask!\n");
    exit(1);
  } else fclose(fic);

  if (argc==4) {
    fic=fopen(weightfile,"r");
    if (fic==0) {
      fprintf(stderr,"Can't open weight file!\n");
      exit(1);
    } else fclose(fic);
  }

  polymask=init_poly(maskfile,outmask);

  n=ct_lines(file);

  ra=dvector(1,n);
  dec=dvector(1,n);

  fic=fopen(file,"r");
  for (i=1;i<=n;i++) {
    fscanf(fic,"%lf %lf\n",&ra[i],&dec[i]);
    ra[i]*=deg2rad;
    dec[i]*=deg2rad;
  }
  fclose(fic);

  if (argc==4) {
    nw=ct_lines(weightfile);
    
    if (nw!=polymask->NpolysAll) {
      fprintf(stderr,"Different number of weights than polygons!\n");
      exit(1);
    }
    
    wei=dvector(1,nw);
    
    fic=fopen(weightfile,"r");
    for (i=1;i<=nw;i++) fscanf(fic,"%lf\n",&wei[i]);
    fclose(fic);
  }
  
  sprintf(file,"%s.out",file);

  fic=fopen(file,"w");
  for (i=1;i<=n;i++) {
    int id=idPoly(polymask,outmask,ra[i],dec[i]);
    if (argc==4) fprintf(fic,"%lf %lf %d %lf\n",ra[i],dec[i],id+1,id+1>0?wei[id+1]:0);
    else fprintf(fic,"%lf %lf %d\n",ra[i],dec[i],id+1);
  }
  fclose(fic);

  free_NodeP(polymask);

  // free_dvector(ra,1,n);
  // free_dvector(dec,1,n);
  // free_dvector(wei,1,nw);

  return 0;
}
