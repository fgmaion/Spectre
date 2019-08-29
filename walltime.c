int start_walltime(){
  start_time = time(NULL);

  return 0;
}


int walltime(char string[]){

  time_t now = time(NULL);

  printf("\n%s:  %lf seconds", string, difftime(now, start_time));
  
  return 0;
}
