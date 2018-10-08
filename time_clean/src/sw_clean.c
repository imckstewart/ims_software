/*
 *  sw_clean.c
 *
 *  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
 *
 *  See the file ../COPYRIGHT for authorship and intellectual property rights.
 *
 */

#include "sw_clean.h"

/*....................................................................*/
boundingBoxType
_calcBoundingBox(_Bool *cleanMask, cleanInfoType *cleanInfo){
  /*
Does what it says on the box - takes a mask image and returns X and Y limits of TRUE values of the mask.
  */
  boundingBoxType boundingBox;
  int i,xi,yi;
  _Bool maskedNotFound;

  boundingBox.xiLo = -1;
  boundingBox.xiHi = -1;
  boundingBox.yiLo = -1;
  boundingBox.yiHi = -1;

  maskedNotFound = TRUE; /* default */
  i = 0;
  for(yi=0;yi<cleanInfo->yiSize;yi++){
    for(xi=0;xi<cleanInfo->xiSize;xi++){
      if(cleanMask[i]){
        if(maskedNotFound){
          maskedNotFound = FALSE;
          boundingBox.xiLo = xi;
          boundingBox.xiHi = xi;
          boundingBox.yiLo = yi;
          boundingBox.yiHi = yi;
        }else{
          if(xi<boundingBox.xiLo) boundingBox.xiLo = xi;
          if(xi>boundingBox.xiHi) boundingBox.xiHi = xi;
          if(yi<boundingBox.yiLo) boundingBox.yiLo = yi;
          if(yi>boundingBox.yiHi) boundingBox.yiHi = yi;
        }
      }

      i++;
    }
  }

return boundingBox; /* Note that -1 values indicate there were no masked pixels found. */
}

/*....................................................................*/
void
_addCleanComp(cleanCompsType *cleanComps, double *cleanValues\
  , const unsigned short nValues, const int yi, const int xi){
  /*
This appends a clean component (vector of values, image pixel location) onto 'cleanComps', which is effectively a list of clean components.
  */
  unsigned short i_us;

  if(cleanComps->nComps >= cleanComps->availableNComps){
    cleanComps->availableNComps += BUFFERSIZE; /* We are going to have a seg fault if cleanComps->bufferSize < 1 */
    cleanComps->comps = realloc(cleanComps->comps, sizeof(*(cleanComps->comps))*cleanComps->availableNComps);
  }

  cleanComps->comps[cleanComps->nComps].values = malloc(sizeof(double)*nValues);
  for(i_us=0;i_us<nValues;i_us++)
    cleanComps->comps[cleanComps->nComps].values[i_us] = cleanValues[i_us];

  cleanComps->comps[cleanComps->nComps].nValues = nValues;
  cleanComps->comps[cleanComps->nComps].xi = xi;
  cleanComps->comps[cleanComps->nComps].yi = yi;

  cleanComps->nComps++;
}

/*....................................................................*/
int
_getUpperTriangleIndex(const unsigned short iSize, const unsigned short ri\
  , const unsigned short ci){
  /*
This is for the situation that we have a square, symmetrical matrix of size iSize for which we want to avoid wasting memory by storing the upper triangle elements only, as a list. The present routine is designed to return the list index corresponding to any row or column of the matrix.
  */
  int k;

  if(ci<ri){
    k = _getUpperTriangleIndex(iSize, ci, ri);
  }else{
    k = ri*(int)iSize - (ri*((int)ri - 1))/2 + ci - ri; 
  }

return k;
}

/*....................................................................*/
errType
doSWMinorClean(double *dirtyImage, double **shiftedBeams, _Bool *cleanMask\
  , cleanInfoType *cleanInfo, cleanCompsType *cleanComps){
  /*
This implements the multi-spectral clean algorithm described in Sault R J & Wieringa M H, A&A Suppl. Ser. 108, 585-594 (1994). Note however that the beams don't have to be constructed from Taylor expansions of source spectra, but can have any form. Note also that this should be considered a minor-cycle clean; for best results, periods of minor-cycling ought to be interspersed with a major cycle in which the accumulated clean components are subtracted in the UV plane.

Some complicated things are done with the sizing and positioning of various sub-images in the present routine. Most of the work is done with N R-images and N*N A-images. The R-images are subsets of correlations between the N beams and the dirty image, whereas the A-images are subsets of correlations between the beams themselves. To explain the size and locations of their respective subsets let us consider a 1-dimensional example. We start with a dirty image of size 16 pixels (plus N dirty images of the same size), schematically indicated as follows.

DI:	|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

Since these images are the result of a discrete Fourier transform, they have periodic boundary conditions. This means that one can cyclic-shift any image by an integer number of pixels without introducing a step, since this is just the equivalent of the FT of the source image (gridded visibilities in the present case) multiplied by complex exponential. Indeed the whole idea of CLEAN is based on the (not quite accurate) assumption that the dirty image is the sum of cyclic-shifted and amplitude-scaled copies of the dirty beam images.

Usually we do not desire to clean sources from the whole image but from some subset of it. In the following example, we specify a clean bounding box from pixel 2 to pixel 10 of the dirty image:

DI:	|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
CBB:	    |_|_|_|_|_|_|_|_|_|

We select this subset of the beam-correlated dirty images for the R-images. Let's denote its width by N_b, which is equal to 9 in the present example.

The cross-correlated beam images are expected to have maximum amplitude in the 0th pixel. An example is as follows:

Full A_i,j:

	 9 7 3 1 1 0 0 0 0 0 0 0 1 1 3 7
	|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

Note that these images are real-valued and symmetrical. For the cropped or subset A-images, there are two possibilities. Either the full image width of M pixels is < 2*N_b, or it is not. If it is, as in the present example, then we must use the full A image without cropping. Schematic comparison of A_{i,j} and CBB is as follows:

A_i,j:	 9 7 3 1 1 0 0 0 0 0 0 0 1 1 3 7
	|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

CBB:	|_|_|_|_|_|_|_|_|_|

To line up A_{i,j} with a peak detected in the Lth pixel of CBB we must cyclic-shift A_{i,j} rightwards by L.

In the other case, that M >= 2*N_b, we can crop the full A images to a size of 2*N_b-1 to save memory. We need to cyclic shift them before doing so, in order to position the peak at the central pixel. Let's example this by supposing that CBB extends from pixel 2 to pixel 5, thus having a width (N_b) of 4 pixels instead of 9. The preliminary arrangement as follows:

A_i,j:	 1 3 7 9 7 3 1
	|_|_|_|_|_|_|_|

However, we can retain the previous simple shift scheme by cyclic-shifting the cropped A-images so the peak once again lies on the 0th pixel:

A_i,j:	 9 7 3 1 1 3 7
	|_|_|_|_|_|_|_|

CBB:	|_|_|_|_|

*NOTE* that since the A-images are symmetrical, we could get even greater space gains, at the cost of a more fiddly index scheme. For the time being I'll resist this temptation.

*NOTE* that the dirty beam array 'shiftedBeams' is named so because it is expected to be supplied in Fourier-natural alignment (=> beam centre at bottom LH pixel).

*NOTE* also that the input 2D arrays are assumed to be in row-major order, i.e. with values arranged in the order xi + cleanInfo->xiSize*yi.

*NOTE* that a TRUE value in the mask image indicates that cleaning is desired for that pixel.

*NOTE* that the calling routine needs to free '*residImage' after use.
  */
  errType err=init_local_err();
  boundingBoxType boundingBox;
  int bbXiSize,bbYiSize,bbISize,aXiSize,aYiSize,aISize,xi1,xi2,yi1,yi2;
  int iSize,halfXiSize,iSizeComplex,i,j,k,xi,yi;
  int status,offsetI,mainImageXiStart,ci,maxXi,maxYi,bbPeakI;
  _Bool isCropped[2],maskedNotFound,hasConverged;
  unsigned short nAImages,i_us,j_us;
  double *workingBeam=NULL,**croppedA=NULL,**croppedR=NULL;
  double *shiftedUncroppedA=NULL,*shiftedUncroppedR=NULL,*equ22Image=NULL;
  double maxPeakValue,*cleanValues=NULL,ftRescale;
  fftw_complex **beamFTs=NULL,*workingBeamFT=NULL,*dirtyImageFT=NULL;
  fftw_plan plan;
  gsl_matrix *A0Mat=NULL,*A0MatInverse=NULL;

  ftRescale = 1.0/cleanInfo->xiSize/cleanInfo->yiSize;

  boundingBox = _calcBoundingBox(cleanMask, cleanInfo);
  if(boundingBox.xiLo<0)
return write_local_err(1, "No masked pixels found.");

  bbXiSize = boundingBox.xiHi + 1 - boundingBox.xiLo;
  bbYiSize = boundingBox.yiHi + 1 - boundingBox.yiLo;
  bbISize = bbYiSize*bbXiSize;

  /* Determine the size to crop the cross-correlation images. We could just keep them uncropped, but this is wasteful of memory if the boundingBox is significantly smaller than the dirty image. What we want is to crop in both axes to (2 times the size of the bounding box) minus 1, or leave uncropped, whichever is smaller.
  */
  if(cleanInfo->xiSize < 2*bbXiSize){ /* don't crop in X */
    aXiSize = cleanInfo->xiSize;
    isCropped[1] = FALSE;
  }else{ /* crop in X */
    aXiSize = 2*bbXiSize - 1;
    isCropped[1] = TRUE;
  }

  if(cleanInfo->yiSize < 2*bbYiSize){ /* don't crop in Y */
    aYiSize = cleanInfo->yiSize;
    isCropped[0] = FALSE;
  }else{ /* crop in Y */
    aYiSize = 2*bbYiSize - 1;
    isCropped[0] = TRUE;
  }

  aISize = aYiSize*aXiSize;

  /* For N beams, technically we should have N*N A-images, but since A_{i,j} = A{j,i} we can make do with N*(N+1)/2.
  */
  nAImages = (cleanInfo->nBeams*(cleanInfo->nBeams + 1))/2;

  /* According to the FFTW notes, if the input array to a real-to-complex transform has N0 rows and N1 columns, therefore total size N = N0*N1, the output complex array should have size N' = N0*(N1/2 + 1).
  */
  iSize = cleanInfo->yiSize*cleanInfo->xiSize;
  halfXiSize = cleanInfo->xiSize/2 + 1;
  iSizeComplex = cleanInfo->yiSize*halfXiSize;

  /*
FT the beams in preparation for the next step. Note that, since the beams are real-valued (i.e. there is no imaginary component), the output images are Hermitian symmetrical, and thus need only be about half the size of the inputs (see the FFTW notes for explanation).

Two comments: (i) I copy the data for each beam to 'workingBeam' because it is recommended in the FFTW doco to create the plan before initializing the arrays; (ii) I create a new plan at each iteration rather than re-using the old because this is also recommended (at least in standard useage) in the case that one or both of the in/out arrays is not the same at each iteration (i.e. that the actual memory mapping changes between 'execute' calls).
  */
  workingBeam = fftw_alloc_real(iSize);
  beamFTs = malloc(sizeof(*beamFTs)*cleanInfo->nBeams);
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    beamFTs[i_us] = fftw_alloc_complex(iSizeComplex);

    plan = fftw_plan_dft_r2c_2d(cleanInfo->yiSize, cleanInfo->xiSize, workingBeam, beamFTs[i_us], FFTW_ESTIMATE);

    for(j=0;j<iSize;j++)
      workingBeam[j] = shiftedBeams[i_us][j];

    fftw_execute(plan);
    fftw_destroy_plan(plan);
  }

  /* Now FT the dirty image:
  */
  dirtyImageFT = fftw_alloc_complex(iSizeComplex);
  plan = fftw_plan_dft_r2c_2d(cleanInfo->yiSize, cleanInfo->xiSize, workingBeam, dirtyImageFT, FFTW_ESTIMATE);

  for(j=0;j<iSize;j++)
    workingBeam[j] = dirtyImage[j];

  fftw_execute(plan);
  fftw_destroy_plan(plan);

  fftw_free(workingBeam);

  /* Malloc a lot of necessary arrays.
  */
  croppedA = malloc(sizeof(*croppedA)*nAImages);
  i = 0;
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    for(j_us=i_us;j_us<cleanInfo->nBeams;j_us++){
      croppedA[i] = malloc(sizeof(**croppedA)*aISize);
      i++;
    }
  }

  A0Mat = gsl_matrix_alloc(cleanInfo->nBeams, cleanInfo->nBeams);

  /*
Generate correlation images, then crop and cshift them if necessary. Note that since the input images to the correlations are real, the result of each correlation must be real as well.

Note that, with the use now of a boundingBox-wide padding border, the beams (which begin with [0,0] as the MS pixel, i.e. with beam centre located at the origin) after cropping should end up with [cropYiMSP,cropXiMSP] as the MSP.

It is convenient also to load the matrix of (0,0) values at the same time.
  */
  workingBeamFT = fftw_alloc_complex(iSizeComplex);
  shiftedUncroppedA = fftw_alloc_real(iSize);
  plan = fftw_plan_dft_c2r_2d(cleanInfo->yiSize, cleanInfo->xiSize, workingBeamFT, shiftedUncroppedA, FFTW_ESTIMATE);
  i = 0; /* Index in list of A-images. */
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    for(j_us=i_us;j_us<cleanInfo->nBeams;j_us++){
      /* Cross-correlate the beams by multiplying the conjugate of one FT by the other FT, then back-transforming. 
      */
      j = 0; /* Pixel number in flat complex (i.e. half-size) image. */
      for(yi=0;yi<cleanInfo->yiSize;yi++){
        for(xi=0;xi<halfXiSize;xi++){
          /* Correlation between (xR + i*xI) and (yR + i*yI) in Fourier space is
		(xR - i*xI)*(yR + i*yI) = (xR*yR + xI*yI) + i*(xR*yI - xI*yR)
          */
          workingBeamFT[j][0] = beamFTs[i_us][j][0]*beamFTs[j_us][j][0] + beamFTs[i_us][j][1]*beamFTs[j_us][j][1];
          workingBeamFT[j][1] = beamFTs[i_us][j][0]*beamFTs[j_us][j][1] - beamFTs[i_us][j][1]*beamFTs[j_us][j][0];
          j++;
        }
      }
      fftw_execute(plan); /* The back-transform of workingBeamFT, outputs to shiftedUncroppedA. */

      gsl_matrix_set(A0Mat, (size_t)i_us, (size_t)j_us, shiftedUncroppedA[0]*ftRescale);

      j = 0; /* Now pixel number in each flat A-image. */
      for(yi=0;yi<aYiSize;yi++){
        if(isCropped[0]){
          yi1 = mod((yi + bbYiSize - 1), aYiSize);
          yi2 = mod((yi1 - (bbYiSize - 1)), cleanInfo->yiSize);
          offsetI = yi2*cleanInfo->xiSize;
        }else{
          offsetI = yi*cleanInfo->xiSize;
        }

        for(xi=0;xi<aXiSize;xi++){
          if(isCropped[1]){
            xi1 = mod((xi + bbXiSize - 1), aXiSize);
            xi2 = mod((xi1 - (bbXiSize - 1)), cleanInfo->xiSize);
            k = offsetI + xi2;
          }else
            k = offsetI + xi;

          croppedA[i][j] = shiftedUncroppedA[k]*ftRescale;

          j++;
        }
      }

    i++;
    }
  }
  fftw_destroy_plan(plan);
  fftw_free(shiftedUncroppedA);

  /* Fill in the rest of A0Mat.
  */
  for(i_us=1;i_us<cleanInfo->nBeams;i_us++){
    for(j_us=0;j_us<i_us;j_us++)
      gsl_matrix_set(A0Mat, (size_t)i_us, (size_t)j_us, gsl_matrix_get(A0Mat, (size_t)j_us, (size_t)i_us));
  }

  croppedR = malloc(sizeof(*croppedR)*cleanInfo->nBeams);
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++)
    croppedR[i_us] = malloc(sizeof(**croppedR)*bbISize);

  shiftedUncroppedR = fftw_alloc_real(iSize);
  plan = fftw_plan_dft_c2r_2d(cleanInfo->yiSize, cleanInfo->xiSize, workingBeamFT, shiftedUncroppedR, FFTW_ESTIMATE);
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    /* Correlate the beam with the dirty image. 
    */
    i = 0;
    for(yi=0;yi<cleanInfo->yiSize;yi++){
      for(xi=0;xi<halfXiSize;xi++){
        /* Correlation between (xR + i*xI) and (yR + i*yI) in Fourier space is
		(xR - i*xI)*(yR + i*yI) = (xR*yR + xI*yI) + i*(xR*yI - xI*yR)
        */
        workingBeamFT[i][0] = beamFTs[i_us][i][0]*dirtyImageFT[i][0] + beamFTs[i_us][i][1]*dirtyImageFT[i][1];
        workingBeamFT[i][1] = beamFTs[i_us][i][0]*dirtyImageFT[i][1] - beamFTs[i_us][i][1]*dirtyImageFT[i][0];
        i++;
      }
    }
    fftw_execute(plan); /* The back-transform, fills shiftedUncroppedR. */

    i = 0;
    for(yi=0;yi<bbYiSize;yi++){
      mainImageXiStart = (yi + boundingBox.yiLo)*cleanInfo->xiSize + boundingBox.xiLo;
      for(xi=0;xi<bbXiSize;xi++){
        croppedR[i_us][i] = shiftedUncroppedR[mainImageXiStart + xi]*ftRescale;
        i++;
      }
    }
  }
  fftw_destroy_plan(plan);
  fftw_free(shiftedUncroppedR);
  fftw_free(workingBeamFT);
  fftw_free(dirtyImageFT);

  /* Invert the A0 matrix. Validity of the Cholesky inversion relies on A0 being positive-definite as well as symmetric.
  */
  A0MatInverse = A0Mat; /* Just gives us another name for the same memory space, which the GSL routines will fill with the matrix inverse. */
  status = gsl_linalg_cholesky_decomp1(A0MatInverse);
  if(status){
    snprintf(errMessage, ERR_STR_LEN, "Cholesky decomposition failed with error %d.", status);
return write_local_err(2, errMessage);
  }
  status = gsl_linalg_cholesky_invert(A0MatInverse);
  if(status){
    snprintf(errMessage, ERR_STR_LEN, "Inversion of A0 matrix failed with error %d.", status);
return write_local_err(3, errMessage);
  }

  /* Malloc some arrays we need in cleaning:
  */
  equ22Image = malloc(sizeof(*equ22Image)*bbISize);
  cleanValues = malloc(sizeof(*cleanValues)*cleanInfo->nBeams);

  /* Now we start the clean iterations:
  */
  hasConverged = FALSE;
  ci = 0;
  while(!hasConverged && ci<cleanInfo->maxNumIter){
    printf("Clean iteration %d\n", ci);

    /* Evaluate 'equation 22':
    */
    i = 0;
    for(yi=0;yi<bbYiSize;yi++){
      for(xi=0;xi<bbXiSize;xi++){
        equ22Image[i] = 0.0;

        for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
          for(j_us=0;j_us<cleanInfo->nBeams;j_us++){
            equ22Image[i] += croppedR[i_us][i]*croppedR[j_us][i]*gsl_matrix_get(A0MatInverse, (size_t)i_us, (size_t)j_us);
          }
        }

        i++;
      }
    }

    /* Find the location and value of the maximum in the (masked) equ22Image:
    */
    maxPeakValue = -1.0;
    maxXi = -1;
    maxYi = -1;
    maskedNotFound = TRUE; /* default */
    i = 0;
    for(yi=0;yi<bbYiSize;yi++){
      mainImageXiStart = (yi + boundingBox.yiLo)*cleanInfo->xiSize + boundingBox.xiLo;
      for(xi=0;xi<bbXiSize;xi++){
        if(cleanMask[mainImageXiStart + xi]){
          if(maskedNotFound || equ22Image[i] > maxPeakValue){
            maxPeakValue = equ22Image[i];
            maxXi = xi;
            maxYi = yi;
            if(maskedNotFound)
              maskedNotFound = FALSE;
          }
        }
        i++;
      }
    }

//**** convergence criterion

    /* Find the location of the peak value in the croppedR image:
    */
    bbPeakI = maxYi*bbXiSize + maxXi;

    /*
Invert the equation

	A0Mat * coeffs = croppedR[:][bbPeakI]

to obtain the 'a' coefficients, then store them:
    */
    for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
      cleanValues[i_us] = 0.0;
      for(j_us=0;j_us<cleanInfo->nBeams;j_us++){
        cleanValues[i_us] += gsl_matrix_get(A0MatInverse, (size_t)i_us, (size_t)j_us)*croppedR[j_us][bbPeakI];
      }
      cleanValues[i_us] *= cleanInfo->cleanLoopGain;
    }

    _addCleanComp(cleanComps, cleanValues, cleanInfo->nBeams, maxYi + boundingBox.yiLo, maxXi + boundingBox.xiLo);

    /* Subtract cleanLoopGain * components from the vector of Rs.
    */
    i = 0;
    for(yi=0;yi<bbYiSize;yi++){
      yi1 = mod(yi - maxYi, aYiSize);
      offsetI = yi1*aXiSize;
      for(xi=0;xi<bbXiSize;xi++){
        xi1 = mod(xi - maxXi, aXiSize);
        j = offsetI + xi1;

        for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
          for(j_us=0;j_us<cleanInfo->nBeams;j_us++){
            k = _getUpperTriangleIndex(cleanInfo->nBeams, i_us, j_us);
            croppedR[i_us][i] -= cleanValues[j_us]*croppedA[k][j];
          }
        }

        i++;
      }
    }

    ci++;
  } /* End clean iteration loop. */

  free(cleanValues);
  free(equ22Image);

  for(i_us=0;i_us<cleanInfo->nBeams;i_us++)
    free(croppedR[i_us]);
  free(croppedR);

  gsl_matrix_free(A0Mat);

  i = 0;
  for(i_us=0;i_us<cleanInfo->nBeams;i_us++){
    for(j_us=i_us;j_us<cleanInfo->nBeams;j_us++){
      free(croppedA[i]);
      i++;
    }
  }
  free(croppedA);

  for(i_us=0;i_us<cleanInfo->nBeams;i_us++)
    fftw_free(beamFTs[i_us]);
  free(beamFTs);

return err;
}

