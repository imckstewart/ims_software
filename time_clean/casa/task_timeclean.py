import os

try:
  from taskinit import *
except ImportError:
  class DummyCasaLog:
    originStr=''

    def origin(self, originStr):
      self.originStr = originStr

    def post(self, errStr, priority='Message'):
      print "%s: %s" % (priority, errStr)

  casalog = DummyCasaLog()

def timeclean(vis,timeseries,imagename,niter,gain,threshold,mask,imsize,cell,weighting,robust,uvtaper,outertaper,innertaper,restoringbeam,noise,npixels,nterms):
  """
The script should perform the following operations:
	- Import the observation MS. Then:
		* Grid it.
		* FT to the sky plane (=>dirty image).

	- Import a time series of the flux, optionally also one of the spectral index.

	- Construct 1 or more 'beam visibility sets' (BVS) from the time series and other information.

	- For each BVS:
		* Grid it.
		* FT it to the sky plane (=> dirty beam).

	- Perform Sault-Wieringa clean.

	- Export clean components and/or restored image.
  """

  _doTest = False

  if _doTest: print '>>> Entering task_timeclean.timeclean()'

  #**** do stuff

  return

