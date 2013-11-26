#! /usr/bin/env python

'''
ABOUT:
This procedure applies the pixel-area map (PAM) correction to WFC3 images.


DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
Aug. 2013: Original script (v0.1).


FUTURE IMPROVEMENTS:

USE:
python run_make_PAMcorr.py
'''

__author__='D.M. HAMMER'
__version__= 0.2


import pyraf, os, glob, argparse, pdb, pyfits, pylab, fileinput, shutil, scipy
import numpy as np

if __name__=='__main__':
        parser = argparse.ArgumentParser(description='Apply the PAM to image(s).')
        parser.add_argument('-im', '--images',default='*crclean.fits', type=str, help='Input fits file(s). \
                                Default is all CRCLEAN images in working directory.')
	parser.add_argument('-out', '--outfile', default='.pam', type=str, help='Name of the output PAM-corrected image. \
				If it starts with "." then it is assumed that we will append that extension to each input image. \
				Only single images are allowed a name that does NOT start with ".".')
        options = parser.parse_args()


	#---Assign name(s) of output PAM-corrected images
	ims = glob.glob(options.images)
	nameflag = options.outfile[0] != '.'
        if (nameflag & (len(ims) != 1)): raise Exception ('Output filenames for a stack of images must be input as an extension to input filename only, e.g. ".corr.".')


	# --Iterate over each image and apply PAM
	for im in ims:
		if nameflag: oname = options.outfile
		else: oname = im + options.outfile

	  	#---Read in fits image -- identify the SCI extensions
		shutil.copy(im,oname)
		hdulist = pyfits.open(oname,mode='update')
		instrume = hdulist[0].header['INSTRUME']
		detector  = hdulist[0].header['DETECTOR']

		sciext = []
		for ff in xrange(len(hdulist)):
			if hdulist[ff].name == 'SCI': sciext.append(ff)

		# - perform correction for each SCI extension
		for ext in sciext:
			fhdr = hdulist[ext].header
			fdata = hdulist[ext].data
			chip = fhdr['CCDCHIP']
			sizaxis1 = fhdr['SIZAXIS1']
			if sizaxis1 != fhdr['NAXIS1']: pdb.set_trace()
			sizaxis2 = fhdr['SIZAXIS2']
			if sizaxis2 != fhdr['NAXIS2']: pdb.set_trace()
			LTV1 = fhdr['LTV1']
                	LTV2 = fhdr['LTV2']

			# initialize the x/y detector coordinates of this image (will be 2048,4096 for full frame; subarray will vary).
			x0 = LTV1
			y0 = LTV2
			x1 = LTV1 + sizaxis1
			y1 = LTV2 + sizaxis2

			# apply the PAM
		        if instrume == 'WFC3':
				if detector == 'UVIS':
					if chip == 1:
						dam=pyfits.getdata('/user/hammer/FF_TEST/UVIS1wfc3_map.fits')
						fdata_cnts = fdata * dam[y0:y1,x0:x1]
					elif chip == 2:
                        			dam=pyfits.getdata('/user/hammer/FF_TEST/UVIS2wfc3_map.fits')
                        			fdata_cnts = fdata * dam[y0:y1,x0:x1]
					else: raise Exception('Case not handled.')
				elif detector == 'IR':
                                        dam=pyfits.getdata('/user/hammer/FF_TEST/ir_wfc3_map.fits')
                        		fdata_cnts = fdata * dam[y0:y1,x0:x1]
                                else: raise Exception('Case not handled.')
			elif instrume == 'ACS':
				if detector == 'WFC':
                                	if chip == 1:
                                        	dam=pyfits.getdata('/user/hammer/FF_TEST/wfc1_pam.fits')
                                        	fdata_cnts = fdata * dam[y0:y1,x0:x1]
                                	elif chip == 2:
                                                dam=pyfits.getdata('/user/hammer/FF_TEST/wfc2_pam.fits')
                                        	fdata_cnts = fdata * dam[y0:y1,x0:x1]
                                	else: raise Exception('Case not handled.')
                                else: raise Exception('Case not handled.')
                        else: raise Exception('Case not handled.')

			hdulist[ext].data = fdata_cnts

