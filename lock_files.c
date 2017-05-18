int delete_lockfile(){
  int    status = 0; 
  char   sys_command[200];

  sprintf(sys_command, "rm -r /home/mjw/IO_lock;");

  status = system(sys_command);

  if(status == 0){
    printf("\n\nIO lock dir. destroyed.");
  }

  else{
    printf("\n\nFailed to destroy IO lock dir.");
  }
  
  return 0;
}
