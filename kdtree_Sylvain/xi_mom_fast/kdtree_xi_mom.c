#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libkdtree.h"
#include "mfunct.h"

#define deg2rad M_PI/180
#define Dh 299792.458/H_o

double H_o, omega_m, omega_v;

/***! Distance and volume !***/
double invE(double z){
    return 1./(sqrt(omega_m*pow(1. + z, 3.) + omega_v));
}

double Dc(double z){
    int j;
    double xm, dx, s = 0.0;

    static double x[]={0.0, 0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
    static double w[]={0.0, 0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

    xm = 0.5*z;
    
    // s=0;
    for (j=1; j<=5; j++) {
        dx = xm*x[j];
        s += w[j]*(invE(xm+dx)+invE(xm-dx));
    }
    s*=xm;

    return Dh*s;
}

double Pl_2 (double x){
    // Second order Legendre polynomial
    return (3.0*x*x-1.0)*0.5;
}


double Pl_4 (double x){
    // Fourth order Legendre polynomial
    return (35.0*x*x*x*x-30.0*x*x+3.0)/8.0;
}

int xi0_comp(int nlogbins, int nlinbins, double linbinsz, double **xi_s_mu, double *xi){
  int i, j;
  
  double mu;

  for(i=1; i<=nlogbins; i++){
    xi[i]=0;
    
    for(j=1; j<=nlinbins; j++) {
      
      mu=linbinsz*0.5+(j-1)*linbinsz;
      
      if (xi_s_mu[i][j]<-1) xi_s_mu[i][j]=0;
      
      if (mu>=0 && mu<=1) xi[i]+=0.5*linbinsz*xi_s_mu[i][j]*2;
    }
  } 
  
  return 0;
}


int xi2_comp(int nlogbins,int nlinbins,double linbinsz,double **xi_s_mu,double *xi) 
{
  int i,j;
  double mu;

  for (i=1;i<=nlogbins;i++) {
    xi[i]=0;
    for (j=1;j<=nlinbins;j++) {
      mu=linbinsz*0.5+(j-1)*linbinsz;
      if (xi_s_mu[i][j]<-1) xi_s_mu[i][j]=0;
      if (mu>=0 && mu<=1) xi[i]+=2.5*linbinsz*xi_s_mu[i][j]*Pl_2(mu)*2;
    }
  }
  return 0;
}

int xi4_comp(int nlogbins,int nlinbins,double linbinsz,double **xi_s_mu,double *xi) 
{
  int i,j;
  double mu;

  for (i=1;i<=nlogbins;i++) {
    xi[i]=0;
    for (j=1;j<=nlinbins;j++) {
      mu=linbinsz*0.5+(j-1)*linbinsz;
      if (xi_s_mu[i][j]<-1) xi_s_mu[i][j]=0;
      if (mu>=0 && mu<=1) xi[i]+=4.5*linbinsz*xi_s_mu[i][j]*Pl_4(mu)*2;
    }
  }
  return 0;
}

void readCat(FILE *fileIn, int x0col, int x1col, int x2col, int wcol, double *ra, double *dec, double *z, double *w, size_t *N){
  size_t i,Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  rewind(fileIn);

  *N = 0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL) if(*(line) != '#') (*N)++;
  
  rewind(fileIn);

  i=0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
    if(getStrings(line,item," ",&Ncol)){
      ra[++i]  = getDoubleValue(item,x0col)*deg2rad;
      dec[i] = getDoubleValue(item,x1col)*deg2rad;
      z[i] = getDoubleValue(item,x2col);
      w[i] = getDoubleValue(item,wcol);
    }
  }
}

void landy_szalay(int nlogbins,int nlinbins,double efngg,double efnrr,double efngr,double **ggsmu,double **rrsmu,double **grsmu,double **xi_s_mu)
{
  int i,j;
  double norm1,norm2;
  
  norm1=efnrr/efngg;
  norm2=efnrr/efngr;

  for (i=1;i<=nlogbins;i++) {
    for (j=1;j<=nlinbins;j++) {
      if (rrsmu[i][j]==0) xi_s_mu[i][j]=0;
      else xi_s_mu[i][j]=norm1*ggsmu[i][j]/rrsmu[i][j]-2*norm2*grsmu[i][j]/rrsmu[i][j]+1;
    }
  }
}

//-----------------------------------------------------------------------

int main(int argc, char **argv)
{
  FILE *fileIn1,*fileIn2,*fileOut;
  int x0col1,x1col1,x2col1,x0col2,x1col2,x2col2,wcol1,wcol2;
  size_t i,j,N1,N2;

  double startlog,zerolog,logbinsz,nlogbins;
  double startlin,zerolin,linbinsz,nlinbins;
  
  double *ra,*dec,*z,*w;
  double *ra_r,*dec_r,*z_r,*w_r;

  double **gg,**rr,**gr,**xi,*xi0,*xi2,*xi4,ngg,ngr,nrr;

  // Initialization ------------------------------------------//

  for(i=0;i<argc;i++){
    if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc < 11){
      printf("Usage:  %s <file_gal> <id_x> <id_y> <id_z> <id_weight> <file_rnd> <id_x> <id_y> <id_z> <id_weight>\n",argv[0]);
      exit(1);
    }
  }

  // Catalogues-------------------------------------------------------------//

  N1=ct_lines(argv[1]);
  N2=ct_lines(argv[6]);

  ra  = dvector(1,N1);
  dec = dvector(1,N1);
  z   = dvector(1,N1);
  w   = dvector(1,N1);
  ra_r  = dvector(1,N2);
  dec_r = dvector(1,N2);
  z_r   = dvector(1,N2);
  w_r   = dvector(1,N2);

  fileIn1 = fopenAndCheck(argv[1],"r");   
  x0col1 = atoi(argv[2]);
  x1col1 = atoi(argv[3]);
  x2col1 = atoi(argv[4]);
  wcol1 = atoi(argv[5]);
  fileIn2 = fopenAndCheck(argv[6],"r"); 
  x0col2 = atoi(argv[7]);
  x1col2 = atoi(argv[8]);
  x2col2 = atoi(argv[9]);
  wcol2 = atoi(argv[10]);

  printf("Reading %s...\n",argv[1]);
  readCat(fileIn1,x0col1,x1col1,x2col1,wcol1,ra,dec,z,w,&N1);
  printf("Data  : %ld galaxies\n",N1);
  readCat(fileIn2,x0col2,x1col2,x2col2,wcol2,ra_r,dec_r,z_r,w_r,&N2);
  printf("Random: %ld galaxies\n",N2);

  fclose(fileIn1);
  fclose(fileIn2);

  // Cosmo   ---------------------------------------------------------------//

  H_o=70.;
  omega_m=0.27;
  omega_v=0.73;
  for (i=1;i<=N1;i++) z[i]=Dc(z[i]);
  for (i=1;i<=N2;i++) z_r[i]=Dc(z_r[i]);

  // Binning ---------------------------------------------------------------//

  zerolog=-1.1;
  nlogbins=370;
  logbinsz=0.01;
  // + half bin. 
  startlog=zerolog+logbinsz*0.5;

  zerolin=0;
  nlinbins=100;
  linbinsz=1.0/(double)nlinbins;
  startlin=zerolin+linbinsz*0.5;

  // Allocation ---------------------------------------------------------------//

  gg=dmatrix(1,nlogbins,1,nlinbins);
  gr=dmatrix(1,nlogbins,1,nlinbins);
  rr=dmatrix(1,nlogbins,1,nlinbins); 
  for(i=1;i<=nlogbins;i++) for(j=1;j<=nlinbins;j++) gg[i][j]=gr[i][j]=rr[i][j]=0;

  // Count pairs ------------------------------------------------//

  printf("Computing DD pairs...\n");
  //comp_kdt_pairs_loglin_smu_w(ra,dec,z,w,N1,ra,dec,z,w,N1,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,gg,&ngg);
  comp_kdt_pairs_loglin_smu(ra,dec,z,N1,ra,dec,z,N1,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,gg,&ngg);

  printf("Computing DR pairs...\n");
  //comp_kdt_pairs_loglin_smu_w(ra,dec,z,w,N1,ra_r,dec_r,z_r,w_r,N2,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,gr,&ngr);
  comp_kdt_pairs_loglin_smu(ra,dec,z,N1,ra_r,dec_r,z_r,N2,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,gr,&ngr);

  printf("Computing RR pairs...\n");
  
  fileIn1=fopen("RR.dat","r");
  if (fileIn1) {
    fscanf(fileIn1,"%lf",&nrr);
    for(i=1;i<=nlogbins;i++) for(j=1;j<=nlinbins;j++) fscanf(fileIn1,"%lf",&rr[i][j]);
    fclose(fileIn1);
  } else {
    //comp_kdt_pairs_loglin_smu_w(ra_r,dec_r,z_r,w_r,N2,ra_r,dec_r,z_r,w_r,N2,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,rr,&nrr);
    comp_kdt_pairs_loglin_smu(ra_r,dec_r,z_r,N2,ra_r,dec_r,z_r,N2,zerolog,nlogbins,logbinsz,zerolin,nlinbins,linbinsz,rr,&nrr);

    fileIn1=fopen("RR.dat","w");
    fprintf(fileIn1,"%.9lf\n",nrr);
    for(i=1;i<=nlogbins;i++) for(j=1;j<=nlinbins;j++) fprintf(fileIn1,"%lf ",rr[i][j]);
    fclose(fileIn1);
  }

  // Compute xi(s,mu) ----------------------------------------//

  printf("Computing Correlation functions\n");

  xi=dmatrix(1,nlogbins,1,nlinbins);
  xi0=dvector(1,nlogbins);
  xi2=dvector(1,nlogbins);
  xi4=dvector(1,nlogbins);

  landy_szalay(nlogbins,nlinbins,ngg,nrr,ngr,gg,rr,gr,xi);

  fileOut=fopen("xismu.dat","w");
  for(i=1;i<=nlogbins;i++) {
    for(j=1;j<=nlinbins;j++) fprintf(fileOut,"%lf %lf %lf\n",pow(10,startlog+(i-1)*logbinsz),startlin+(j-1)*linbinsz,xi[i][j]);
    fprintf(fileOut,"\n");
  }
  fclose(fileOut);
 
  xi0_comp(nlogbins,nlinbins,linbinsz,xi,xi0);
  xi2_comp(nlogbins,nlinbins,linbinsz,xi,xi2);
  xi4_comp(nlogbins,nlinbins,linbinsz,xi,xi4);

  fileOut=fopen("xil.dat","w");
  for(i=1;i<=nlogbins;i++) fprintf(fileOut,"%lf %lf %lf %lf\n",pow(10,startlog+(i-1)*logbinsz),xi0[i],xi2[i],xi4[i]);
  fclose(fileOut);

  // Frees memory --------------------------------------------//

  free_dmatrix(gg,1,nlogbins,1,nlinbins);
  free_dmatrix(gr,1,nlogbins,1,nlinbins);
  free_dmatrix(rr,1,nlogbins,1,nlinbins);
  free_dmatrix(xi,1,nlogbins,1,nlinbins);

  free_dvector(xi0,1,nlogbins);
  free_dvector(xi2,1,nlogbins);
  free_dvector(xi4,1,nlogbins);

  free_dvector(ra,1,N1);
  free_dvector(dec,1,N1);
  free_dvector(z,1,N1);
  free_dvector(w,1,N1);
  free_dvector(ra_r,1,N2);
  free_dvector(dec_r,1,N2);
  free_dvector(z_r,1,N2);
  free_dvector(w_r,1,N2);
  
  return 0;
}
