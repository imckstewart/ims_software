/*
 *  sw_fits_io.c
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#include "sw_clean.h"

/*....................................................................*/
void
processFitsError(int status){
  /*****************************************************/
  /* Print out cfitsio error messages and exit program */
  /*****************************************************/

  if(status){
    fits_report_error(stderr, status); /* print error report */
exit(status);    /* terminate the program, returning error status */
  }

return;
}

/*....................................................................*/
errType
readFitsArrays(char *inFileName, double **dirtyImage, double ***shiftedBeams\
  , int *yiSize, int *xiSize, unsigned short *nBeams){
  /*
For simplicity the images are all stored as separate layers of the one 3-dimensional cube. The dirty image is the 0th layer, the beams start at layer 1.

Note that this can also be used to read the CC image if shiftedBeams is set to NULL.

Note that, although the arrays are stored in FITS in column-major order rather than the row-major which is conventional in C, cfitsio magically does the transpose for us.
  */
  errType err=init_local_err();
  int status=0,bitpix,naxis,anynul=0,i,xi,yi;
  fitsfile *fptr=NULL;
/*The interface to fits_get_img_param() says long* for argument 5 but I get a seg fault unless I use an array.
  long *naxes; */
long naxes[3];
  unsigned short i_us,j_us;
  long fpixels[3],lpixels[3],inc[] = {1,1,1};
  double *bufferArray=NULL;

  fits_open_file(&fptr, inFileName, READONLY, &status);
  if(status){
    fits_report_error(stderr, status);
    snprintf(errMessage, ERR_STR_LEN, "Attempt to open FITS file for reading generated FITS error %d.", status);
return write_local_err(1, errMessage);
  }

  fits_movabs_hdu(fptr, 1, NULL, &status);
  if(status){
    fits_report_error(stderr, status);
    snprintf(errMessage, ERR_STR_LEN, "Attempt to access primary HDU generated FITS error %d.", status);
return write_local_err(2, errMessage);
  }

  fits_get_img_param(fptr, 3, &bitpix, &naxis, naxes, &status);
  if(status){
    fits_report_error(stderr, status);
    snprintf(errMessage, ERR_STR_LEN, "Attempt to get image params generated FITS error %d.", status);
return write_local_err(3, errMessage);
  }

  if(bitpix != FLOAT_IMG && bitpix != DOUBLE_IMG){
    snprintf(errMessage, ERR_STR_LEN, "FITS image had non-useable BITPIX value of %d.", bitpix);
return write_local_err(4, errMessage);
  }

  *xiSize = naxes[0];
  *yiSize = naxes[1];

  bufferArray = malloc(sizeof(*bufferArray)*(*xiSize)*(*yiSize));

  fpixels[0] = 1;
  fpixels[1] = 1;
  fpixels[2] = 1;
  lpixels[0] = naxes[0];
  lpixels[1] = naxes[1];
  lpixels[2] = 1;

  /* If bitpix==FLOAT_IMG, the values are promoted by cfitsio to double.
  */
  fits_read_subset(fptr, TDOUBLE, fpixels, lpixels, inc, 0, bufferArray, &anynul, &status);
  if(status){
    fits_report_error(stderr, status);
    free(bufferArray);

    snprintf(errMessage, ERR_STR_LEN, "Attempt to read image layer 1/%d generated FITS error %d.", (int)(*nBeams+1), status);
return write_local_err(5, errMessage);
  }

  /* Copy the buffer to dirtyImage.
  */
  *dirtyImage = malloc(sizeof(**dirtyImage)*(*xiSize)*(*yiSize));
  i = 0;
  for(yi=0;yi<*yiSize;yi++){
    for(xi=0;xi<*xiSize;xi++){
      (*dirtyImage)[i] = bufferArray[i];
      i++;
    }
  }

  if(shiftedBeams != NULL){
    *nBeams = (unsigned short)naxes[2] - 1;
    *shiftedBeams = malloc(sizeof(**shiftedBeams)*(*nBeams));

    for(i_us=0;i_us<*nBeams;i_us++){
      fpixels[2] = (long)i_us + 2;
      lpixels[2] = (long)i_us + 2;

      fits_read_subset(fptr, TDOUBLE, fpixels, lpixels, inc, 0, bufferArray, &anynul, &status);
      if(status){
        fits_report_error(stderr, status);
        free(*dirtyImage);
        *dirtyImage = NULL;

        free(bufferArray);
        for(j_us=0;j_us<i_us;j_us++)
          free((*shiftedBeams)[j_us]);
        free(*shiftedBeams);
        *shiftedBeams = NULL;

        snprintf(errMessage, ERR_STR_LEN, "Attempt to read image layer %d/%d generated FITS error %d.", (int)i_us+2, (int)(*nBeams+1), status);
return write_local_err(6, errMessage);
      }

      /* Copy the buffer to the beam.
      */
      (*shiftedBeams)[i_us] = malloc(sizeof(**dirtyImage)*(*xiSize)*(*yiSize));
      i = 0;
      for(yi=0;yi<*yiSize;yi++){
        for(xi=0;xi<*xiSize;xi++){
          (*shiftedBeams)[i_us][i] = bufferArray[i];
          i++;
        }
      }
    }
  }

  free(bufferArray);

  fits_close_file(fptr, &status);
  if(status){
    fits_report_error(stderr, status);
    free(*dirtyImage);
    *dirtyImage = NULL;

    for(i_us=0;i_us<*nBeams;i_us++)
      free((*shiftedBeams)[i_us]);
    free(*shiftedBeams);
    *shiftedBeams = NULL;

    snprintf(errMessage, ERR_STR_LEN, "Attempt to close FITS file generated FITS error %d.", status);
return write_local_err(7, errMessage);
  }

return err;
}

/*....................................................................*/
errType
writeImageToFITS(double *outImage, const unsigned short numLayers\
  , const int yiSize, const int xiSize, char *outFileName){
  /*
This is a fairly generic routine which can be used to write any 3D image of double-precision values to FITS.

Note that, although the arrays are stored in FITS in column-major order rather than the row-major which is conventional in C, cfitsio magically does the transpose for us.
  */
  errType err=init_local_err();
  const int naxis=3;
  long naxes[naxis];
  long fpixels[naxis],lpixels[naxis];
  fitsfile *fptr=NULL;
  int status=0,bitpix=DOUBLE_IMG,i,xi,yi,iOffset;
  char bangThenFile[MAX_LEN_FNAME+3]="! ";
  double *bufferArray=NULL;
  unsigned short i_us;

  naxes[0] = xiSize;
  naxes[1] = yiSize;
  naxes[2] = numLayers;

  fits_create_file(&fptr, outFileName, &status);
  if(status!=0){
    warning("Overwriting existing fits file.");
    status=0;
    strcat(bangThenFile, outFileName);
    fits_create_file(&fptr, bangThenFile, &status);
    if(status!=0){
      fits_report_error(stderr, status);
      snprintf(errMessage, ERR_STR_LEN, "Attempt to create FITS file generated FITS error %d.", status);
return write_local_err(1, errMessage);
    }
  }

  fits_create_img(fptr, bitpix, naxis, naxes, &status);
  if(status!=0){
    fits_report_error(stderr, status);
    snprintf(errMessage, ERR_STR_LEN, "Attempt to create primary image generated FITS error %d.", status);
return write_local_err(2, errMessage);
  }

  bufferArray = malloc(sizeof(*bufferArray)*xiSize*yiSize);

  fpixels[0] = 1;
  fpixels[1] = 1;
  lpixels[0] = naxes[0];
  lpixels[1] = naxes[1];

  for(i_us=0;i_us<numLayers;i_us++){
    /* Copy the outImage to the buffer.
    */
    i = 0;
    iOffset = i_us*xiSize*yiSize;
    for(yi=0;yi<yiSize;yi++){
      for(xi=0;xi<xiSize;xi++){
        bufferArray[i] = outImage[iOffset + i];
        i++;
      }
    }

    fpixels[2] = (long)i_us + 1;
    lpixels[2] = (long)i_us + 1;

    fits_write_subset(fptr, TDOUBLE, fpixels, lpixels, bufferArray, &status);
    if(status){
      fits_report_error(stderr, status);
      free(bufferArray);
      snprintf(errMessage, ERR_STR_LEN, "Attempt to write array layer %d generated FITS error %d.", (int)i_us, status);
return write_local_err(3, errMessage);
    }
  }

  free(bufferArray);

  fits_close_file(fptr, &status);
  if(status){
    fits_report_error(stderr, status);
    snprintf(errMessage, ERR_STR_LEN, "Attempt to close FITS file generated FITS error %d.", status);
return write_local_err(4, errMessage);
  }

return err;
}

