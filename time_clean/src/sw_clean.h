/*
 *  sw_clean.h
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#ifndef SW_CLEAN_H
#define SW_CLEAN_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <fftw3.h>
#include <fitsio.h>

#include "local_err.h"

#ifndef TRUE
#define TRUE  (_Bool)1
#endif

#ifndef FALSE
#define FALSE (_Bool)0
#endif

#define BUFFERSIZE	1024
#define MAX_LEN_FNAME	 256

typedef struct {
  int xiSize,yiSize,maxNumIter;
  unsigned short nBeams;
  double cleanLoopGain;
} cleanInfoType;

typedef struct {
  double *values;
  unsigned short nValues;
  int xi,yi;
} cleanCompType;

typedef struct {
  cleanCompType *comps;
  int nComps,availableNComps;
} cleanCompsType;

typedef struct {
  int xiHi,xiLo,yiHi,yiLo;
} boundingBoxType;

/* Functions in sw_utils.c
*/
int	mod(int a, int b);
cleanCompsType initializeCleanComps(void);
void	freeCleanComps(cleanCompsType *cleanComps);
errType	extractInfo(char *infoFileName, int *maxNumIter, double *cleanLoopGain\
  , boundingBoxType **cleanBoxList, int *nCleanBoxes, char *residFileName\
  , char *outCCFileName);
errType	convertCCFromImage(cleanInfoType *cleanInfo, double *ccImage, cleanCompsType *cleanComps);
errType	convertCCToImage(cleanInfoType *cleanInfo, cleanCompsType *cleanComps, double **ccImage);

/* Functions in sw_clean.c
*/
errType	doSWMinorClean(double *dirtyImage, double **shiftedBeams, _Bool *cleanMask\
  , cleanInfoType *cleanInfo, cleanCompsType *cleanComps);

/* Functions in sw_fits_io.c
*/
errType	readFitsArrays(char *inFileName, double **dirtyImage, double ***shiftedBeams\
  , int *yiSize, int *xiSize, unsigned short *nBeams);
errType	writeImageToFITS(double *outImage, const unsigned short numLayers, const int yiSize, const int xiSize, char *outFileName);


#endif /* SW_CLEAN_H */

