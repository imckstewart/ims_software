/*
 *  sw_utils.c
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#include "sw_clean.h"

/*....................................................................*/
int mod(int a, int b){
  /* Implements a modulo function. */

  int r = a % b;
  return r < 0 ? r + b : r;
}

/*....................................................................*/
cleanCompsType
initializeCleanComps(void){
  cleanCompsType cleanComps;

  cleanComps.comps = NULL;
  cleanComps.nComps = 0;
  cleanComps.availableNComps = 0;

return cleanComps;
}

/*....................................................................*/
void
freeCleanComps(cleanCompsType *cleanComps){
  int i;

  if(cleanComps != NULL && cleanComps->comps != NULL){
    for(i=0;i<cleanComps->nComps;i++)
      free(cleanComps->comps[i].values);

    free(cleanComps->comps);
  }
}

/*....................................................................*/
errType
extractInfo(char *infoFileName, int *maxNumIter, double *cleanLoopGain\
  , boundingBoxType **cleanBoxList, int *nCleanBoxes, char *outResidFileName\
  , char *outCCFileName){
  /*
This reads various stuff from the supplied ASCII file 'infoFileName' and returns the values in the relevant arguments. The format of the ASCII infoFile is supposed to be as follows:

maxNumIter # int
cleanLoopGain # float
outResidFileName
outCCFileName
nCleanBoxes # int
xiLo,xiHi,yiLo,yiHi # for the 1st clean box
  .
  .
xiLo,xiHi,yiLo,yiHi # for the nth clean box

As an example:

1000
0.05
../data/test_resids.fits
../data/test_CC.fits
1
461,626,383,567

  */
  errType err=init_local_err();
  FILE *fp;
  int xiLo,xiHi,yiLo,yiHi,i;


  if((fp=fopen(infoFileName, "r"))==NULL) {
    snprintf(errMessage, ERR_STR_LEN, "Error opening the info file.");
return write_local_err(1, errMessage);
  }

  if(fscanf(fp, "%d\n", maxNumIter) != 1){
    fclose(fp);
    snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: maxNumIter.");
return write_local_err(2, errMessage);
  }

  if(fscanf(fp, "%lf\n", cleanLoopGain) != 1){
    fclose(fp);
    snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: cleanLoopGain.");
return write_local_err(3, errMessage);
  }

  if(fscanf(fp, "%s\n", outResidFileName) != 1){
    fclose(fp);
    snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: outResidFileName.");
return write_local_err(4, errMessage);
  }

  if(fscanf(fp, "%s\n", outCCFileName) != 1){
    fclose(fp);
    snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: outCCFileName.");
return write_local_err(5, errMessage);
  }

  if(fscanf(fp, "%d\n", nCleanBoxes) != 1){
    fclose(fp);
    snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: nCleanBoxes.");
return write_local_err(6, errMessage);
  }

  *cleanBoxList = malloc(sizeof(**cleanBoxList)*(*nCleanBoxes));

  for(i=0;i<*nCleanBoxes;i++){
    if(fscanf(fp, "%d,%d,%d,%d\n", &xiLo, &xiHi, &yiLo, &yiHi) != 4){
      free(*cleanBoxList);
      *cleanBoxList = NULL;
      fclose(fp);
      snprintf(errMessage, ERR_STR_LEN, "Error reading from info file: bounds of %dth clean box.", i);
return write_local_err(7, errMessage);
    }
    (*cleanBoxList)[i].xiLo = xiLo;
    (*cleanBoxList)[i].xiHi = xiHi;
    (*cleanBoxList)[i].yiLo = yiLo;
    (*cleanBoxList)[i].yiHi = yiHi;
  }

  fclose(fp);

return err;
}

/*....................................................................*/
errType
convertCCFromImage(cleanInfoType *cleanInfo, double *ccImage, cleanCompsType *cleanComps){
  errType err=init_local_err();
  int xi,yi,ci,i;
  unsigned short i_us;

  /* Count the number of pixels of ccImage for which at least one beam is non-zero. This will be the number of clean components.
  */
  cleanComps->availableNComps = BUFFERSIZE;
  cleanComps->comps = malloc(sizeof(*(cleanComps->comps))*cleanComps->availableNComps);
  cleanComps->nComps = 0;
  for(yi=0;yi<cleanInfo->yiSize;yi++){
    for(xi=0;xi<cleanInfo->xiSize;xi++){
      for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
        i = (i_us*cleanInfo->yiSize + yi)*cleanInfo->xiSize + xi;
        if(ccImage[i] != 0.0){
          if(cleanComps->nComps >= cleanComps->availableNComps){
            cleanComps->availableNComps += BUFFERSIZE;
            cleanComps->comps = realloc(cleanComps->comps, sizeof(*(cleanComps->comps))*cleanComps->availableNComps);
          }

          cleanComps->comps[cleanComps->nComps].xi = xi;
          cleanComps->comps[cleanComps->nComps].yi = yi;
          cleanComps->nComps++;

      break;
        }
      }
    }
  }

  for(ci=0;ci<cleanComps->nComps;ci++){
    cleanComps->comps[ci].values = malloc(sizeof(double)*cleanInfo->nBeams);
    cleanComps->comps[ci].nValues = cleanInfo->nBeams;

    for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
      i = (i_us*cleanInfo->yiSize + cleanComps->comps[ci].yi)*cleanInfo->xiSize + cleanComps->comps[ci].xi;
      cleanComps->comps[ci].values[i_us] = ccImage[i];
    }
  }

return err;
}

/*....................................................................*/
errType
convertCCToImage(cleanInfoType *cleanInfo, cleanCompsType *cleanComps, double **ccImage){
  /*
The present function takes a list of (vector-valued) clean components and transforms these to an image cube. The jth layer of the cube corresponds to the jth elements of the clean components.

For each jth image, the image has the same size and register as the original dirty image. The output image has pixel values which are the cumulative value of the jth elements of all the clean components which have their location at that pixel.
  */

  errType err=init_local_err();
  int xi,yi,i,ci;
  unsigned short i_us;

  *ccImage = malloc(sizeof(**ccImage)*cleanInfo->nBeams*cleanInfo->yiSize*cleanInfo->xiSize);

  /* Set everything to zero as default.
  */
  i = 0;
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    for(yi=0;yi<cleanInfo->yiSize;yi++){
      for(xi=0;xi<cleanInfo->xiSize;xi++){
        (*ccImage)[i] = 0.0;
        i++;
      }
    }
  }

  for(ci=0;ci<cleanComps->nComps;ci++){
    if(cleanComps->comps[ci].nValues != cleanInfo->nBeams){
      snprintf(errMessage, ERR_STR_LEN, "CC %d has %hu values, but %hud expected.", ci, cleanComps->comps[ci].nValues, cleanInfo->nBeams);
return write_local_err(1, errMessage);
    }

    for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
      i = (i_us*cleanInfo->yiSize + cleanComps->comps[ci].yi)*cleanInfo->xiSize + cleanComps->comps[ci].xi;
      (*ccImage)[i] += cleanComps->comps[ci].values[i_us];
    }
  }

return err;
}





