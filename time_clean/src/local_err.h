/*
 *  local_err.h
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#ifndef LOCAL_ERR_H
#define LOCAL_ERR_H

#define ERR_STR_LEN	255

extern char errMessage[ERR_STR_LEN+1];

typedef struct  {
  int status;
  char message[ERR_STR_LEN+1];
} errType;

errType	init_local_err(void);
errType	write_local_err(int status, char *message);
void	error(int exitStatus, char *message);
void	warning(char *message);

#endif /* LOCAL_ERR_H */

