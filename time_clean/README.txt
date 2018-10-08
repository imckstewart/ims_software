20180912 IMS: In this directory I intend to write a CASA task which performs the CLEAN procedure described in Stewart, Fenech & Muxlow, A&A 535, A81 (2011).

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

How to test the core code:
	* CD to the directory that contains the Makefile (aka the base directory). Run 'make' with no arguments. The build should generate no errors, leaving an executable 'sw_clean' in the same directory.

	* CD to the 'test' sub-directory.

	* Run './sw_test.py'.

	* The task should clean the test image in <base dir>/data/ and write two files to the same directory: 'test_CC.fits' and 'test_resids.fits'. The first is a cube, each jth layer of which is an image of the jth element of the vector-valued clean components; the second is the residual image.

	* After you have finished, clean up by running 'make distclean' in the base directory.

