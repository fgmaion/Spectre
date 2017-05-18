int check_environmentvars(void){
  extern char **environ;

  int i = 0;

  while(environ[i]) {
    printf("%s\n", environ[i++]); // prints in form of "variable=value"
  }
  
  return 0;
}

int camb_call(int dononlinear, double redshift){
  int         status = 0;
  double camb_sigma8 = 0.;
  
  char   sys_command[200];

  char       set_Hub[200];
  char      set_Om_v[200];
  char    set_Om_cdm[200];
  char  set_redshift[200];
  char set_nonlinear[200];
  
  // set environment variables; strings are linked and therefore must be new for each instance. 
  sprintf(set_nonlinear, "DONONLINEAR=%d", dononlinear);
   putenv(set_nonlinear);
  
  sprintf(set_Om_cdm, "OM_CDM=%.3lf", Om_m - Om_b);
   putenv(set_Om_cdm);

  sprintf(set_Om_v, "OM_LMB=%.3lf", Om_v);
   putenv(set_Om_v);
 
  sprintf(set_Hub, "HUBBLE=%.2lf", 100.*h);
   putenv(set_Hub);

  sprintf(set_redshift, "REDSHIFT=%.2lf", redshift);
   putenv(set_redshift);
      
  // Issue CAMB call; uses sed to set var in params.ini and runs. 
  sprintf(sys_command, "/home/mjw/CAMB/camb_call.sh");
  
  status = system(sys_command);

  if(status==1){
    printf("\n\nCamb call failed.");

    return 1;
  }
          
  return 0;
}
