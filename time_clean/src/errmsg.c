#include "local_err.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*....................................................................*/
errType
init_local_err(void){
  errType err;

  err.status = 0;
  err.message[0] = '\0';

return err;
}

/*....................................................................*/
errType
write_local_err(int status, char *message){
  errType err;

  err.status = status;
  strncpy(err.message, message, ERR_STR_LEN);

return err;
}

/*....................................................................*/
void
error(int exitStatus, char *message){
  printf("ERROR: %s\n", message);
exit(exitStatus);
}

/*....................................................................*/
void
warning(char *message){
  printf("Warning: %s\n", message);
}

