/*
 *  sw_main.c
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#include "sw_clean.h"

char errMessage[];

/*....................................................................*/
int main(int argc, char *argv[]){
  /*
Arguments are expected to be fitsFileName, infoFileName, [ccImageName].
  */
  errType err=init_local_err();
  cleanInfoType cleanInfo;
  boundingBoxType *cleanBoxList=NULL;
  double *dirtyImage=NULL,**shiftedBeams=NULL,*residImage=NULL,*ccImage=NULL;
  cleanCompsType cleanComps;
  _Bool *cleanMask=NULL;
  unsigned short i_us;
  int nCleanBoxes,xiSize,yiSize,xi,yi,i,j,ci,yi1,xi1,offsetI;
  char residFileName[MAX_LEN_FNAME],outCCFileName[MAX_LEN_FNAME];

  if(argc!=3 && argc!=4)
error(1, "There should be either 2 or 3 command-line arguments.");

  err = readFitsArrays(argv[1], &dirtyImage, &shiftedBeams\
    , &cleanInfo.yiSize, &cleanInfo.xiSize, &cleanInfo.nBeams);
  if(err.status != 0)
error(err.status, err.message);

  err = extractInfo(argv[2], &cleanInfo.maxNumIter, &cleanInfo.cleanLoopGain, &cleanBoxList, &nCleanBoxes, residFileName, outCCFileName);
  if(err.status != 0)
error(err.status, err.message);

  /* Read CC if given.
  */
  if(argc==4){
    err = readFitsArrays(argv[3], &ccImage, NULL, &yiSize, &xiSize, NULL);
    if(err.status != 0)
error(err.status, err.message);

    if(xiSize != cleanInfo.xiSize || yiSize != cleanInfo.yiSize){
error(99, "The CC image is not the same size as the dirty image.");
    }

    err = convertCCFromImage(&cleanInfo, ccImage, &cleanComps);
    if(err.status != 0)
error(err.status, err.message);

    free(ccImage);

  }else
    cleanComps = initializeCleanComps();

  /* Construct the clean mask.
  */
  cleanMask = malloc(sizeof(*cleanMask)*cleanInfo.xiSize*cleanInfo.yiSize);
  if(nCleanBoxes <= 0){
    /* Set entire mask TRUE - i.e., clean over the whole image.
    */
    for(yi=0;yi<cleanInfo.yiSize;yi++){
      for(xi=0;xi<cleanInfo.xiSize;xi++){
        cleanMask[yi*cleanInfo.xiSize + xi] = TRUE;
      }
    }

  }else{
    /* Set entire mask FALSE except for those parts that lie within a box.
    */
    for(yi=0;yi<cleanInfo.yiSize;yi++){
      for(xi=0;xi<cleanInfo.xiSize;xi++){
        cleanMask[yi*cleanInfo.xiSize + xi] = FALSE;
      }
    }

    for(i=0;i<nCleanBoxes;i++){
      for(yi=cleanBoxList[i].yiLo;yi<=cleanBoxList[i].yiHi;yi++){
        for(xi=cleanBoxList[i].xiLo;xi<=cleanBoxList[i].xiHi;xi++){
          cleanMask[yi*cleanInfo.xiSize + xi] = TRUE;
        }
      }
    }
  }

  /* The main routine.
  */
  err = doSWMinorClean(dirtyImage, shiftedBeams, cleanMask, &cleanInfo, &cleanComps);
  if(err.status != 0)
error(err.status, err.message);

  /* Outputs.
  */
  err = convertCCToImage(&cleanInfo, &cleanComps, &ccImage);
  if(err.status != 0)
error(err.status, err.message);

  freeCleanComps(&cleanComps);

  err = writeImageToFITS(ccImage, cleanInfo.nBeams, cleanInfo.yiSize, cleanInfo.xiSize, outCCFileName);
  if(err.status != 0)
error(err.status, err.message);

  /* Now convert the CC back again (for quicker construction of the residuals image). In the usual case in which many CC are found at each spatial pixel, this reduces each of these to a single CC.
  */
  err = convertCCFromImage(&cleanInfo, ccImage, &cleanComps);
  if(err.status != 0)
error(err.status, err.message);

  free(ccImage);

  /* Generate the residuals image.
  */
  residImage = malloc(sizeof(*residImage)*cleanInfo.xiSize*cleanInfo.yiSize);
  i = 0;
  for(yi=0;yi<cleanInfo.yiSize;yi++){
    for(xi=0;xi<cleanInfo.xiSize;xi++){
      residImage[i] = dirtyImage[i];
      i++;
    }
  }

  for(ci=0;ci<cleanComps.nComps;ci++){
    if(cleanComps.comps[ci].nValues != cleanInfo.nBeams){
      snprintf(errMessage, ERR_STR_LEN, "CC %d has %hu values, but %hud expected.", ci, cleanComps.comps[ci].nValues, cleanInfo.nBeams);
error(1, errMessage);
    }

    i = 0;
    for(yi=0;yi<cleanInfo.yiSize;yi++){
      yi1 = mod(yi - cleanComps.comps[ci].yi, cleanInfo.yiSize);
      offsetI = yi1*cleanInfo.xiSize;
      for(xi=0;xi<cleanInfo.xiSize;xi++){
        xi1 = mod(xi - cleanComps.comps[ci].xi, cleanInfo.xiSize);
        j = offsetI + xi1;

        for(i_us=0;i_us<cleanInfo.nBeams;i_us++)
          residImage[i] -= cleanComps.comps[ci].values[i_us]*shiftedBeams[i_us][j];

        i++;
      }
    }
  }

  err = writeImageToFITS(residImage, 1, cleanInfo.yiSize, cleanInfo.xiSize, residFileName);
  if(err.status != 0)
error(err.status, err.message);

  free(residImage);
  free(cleanMask);
  freeCleanComps(&cleanComps);
  free(cleanBoxList);
  free(dirtyImage);
  for(i_us=0;i_us<cleanInfo.nBeams;i_us++)
    free(shiftedBeams[i_us]);
  free(shiftedBeams);

return 0;
}

