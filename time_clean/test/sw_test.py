#!/usr/bin/env python

#  sw_test.py
#
#  This file is part of time_clean, a CASA task for cleaning images in which the source varies significantly in flux over the time of the observation.
#
#  See the file ../COPYRIGHT for authorship and intellectual property rights.

import os
import tempfile
import subprocess as sub
import signal


def writeInfoFile(infoFileName):
  maxNumIter    = 1000
  cleanLoopGain = 0.05
  residFileName = "../data/test_resids.fits"
  outCCFileName = "../data/test_CC.fits"
  nCleanBoxes   = 1
  xiLo = 461
  xiHi = 626
  yiLo = 383
  yiHi = 567

  fp = open(infoFileName, 'w')

  fp.write("%d\n" % maxNumIter)
  fp.write("%f\n" % cleanLoopGain)
  fp.write("%s\n" % residFileName)
  fp.write("%s\n" % outCCFileName)
  fp.write("%d\n" % nCleanBoxes)
  fp.write("%d,%d,%d,%d\n" % (xiLo,xiHi,yiLo,yiHi))

  fp.close()

def test_SW():
  tempFileObj = tempfile.NamedTemporaryFile(prefix='info_', delete=False)
  infoFileName = tempFileObj.name
  writeInfoFile(infoFileName)

  fitsFileName = "../data/image_and_beams.fits"
  execName = "../sw_clean"

  sigInterp = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

  print "Command:\n%s %s %s" % (execName, fitsFileName, infoFileName)

  try:
    p = sub.Popen([execName, fitsFileName, infoFileName], bufsize=0, stdout=sub.PIPE)
  except OSError:
    raise ValueError("Cannot execute %s (not in PATH, or wrong permissions?)" % (execName))
  except:
    raise ValueError("Execution of %s failed with an unrecognized exception." % (execName))

  while True:
    stdoutdata = p.stdout.readline()[:-1]
    if not stdoutdata is None and stdoutdata!='':
      print stdoutdata

    status = p.poll()
    if not status is None:
      break

  if status==0:
    print "Run of %s is complete." % (execName)
  elif status<0:
    signalStr = sigInterp.get(-status, "unnamed signal: %d" % -status)
    print "Run of %s terminated by signal: %s" % (execName, signalStr)
  else:
    print "Run of %s exited with nonzero status %s." % (execName, status)

  if os.path.exists(infoFileName):
    print "Unlinking %s." % (infoFileName)
    os.unlink(infoFileName)

  # If we get to here, sw_clean is finished or has crashed, so return.
  return

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__':
  test_SW()

