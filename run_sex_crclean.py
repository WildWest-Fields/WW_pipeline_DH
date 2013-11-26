#! /usr/bin/env python

'''
ABOUT:
This script runs SExtractor on list of images (crclean in current directory is default).

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2013

HISTORY:
Jul 2013: Original script (v0.1).
Aug 2013: Modified to run on any instrument's image (e.g., double-chip ACS/optical or single-chip WFC/IR).


FUTURE IMPROVEMENTS:

USE:
python run_sex_crclean
'''

__author__='D.M. HAMMER'
__version__= 0.1


import os, glob, argparse, pdb, pylab, pyfits
import numpy as np
from subprocess import call


def get_wfc3_zeropoint(filter):
        # array of WFC3 zeropoints (not all here - add as needed)
        zp = {'F225W':24.0403, 'F275W':24.1305, 'F336W':24.6682,'F390W': 25.3562, 'F438W':24.8206,'F475W':25.6799,'F555W':25.7906,'F606W':26.0691,'F814W':25.0985,'F850LP':23.8338,'F105W':26.2687,'F125W':26.2303,'F140W':26.4524,'F160W':25.9463}
        if zp.has_key(filter.upper()): return zp[filter.upper()]
        else: raise Exception('Zeropoint is not specified for this filter: '+filter)


def get_acs_zeropoint(hdr):
        PHOTFLAM = hdr['PHOTFLAM']
        PHOTPLAM = hdr['PHOTPLAM']
        zeropt = -2.5*np.log10(PHOTFLAM) - 5.0*np.log10(PHOTPLAM) - 2.408
        return zeropt


if __name__=='__main__':

    # -- Parse input parameters
    #--------------------------
    parser = argparse.ArgumentParser(description='Run SExtractor on specified images.')
    parser.add_argument('-im', '--images',default='*crclean.fits', type=str, help='Input fits file(s). \
                        Default is all CR-cleaned science images in working directory.')
    parser.add_argument('-max', '--maxobj',default=100, type=int, help='Input maximum number of SExtractor detections to use for tweakreg matching (each chip). \
                        Default is 100 objects per chip.')
    parser.add_argument('-th', '--threshold',default=10, type=int, help='Input threshold used for BOTH detect_thresh and analysis_thresh. \
                        Default is 10 - tests showed better tweakreg residuals when using high thresh for even bright objects.')

    options = parser.parse_args()
    imlist = glob.glob(options.images)
    imlist.sort()
    maxobj = options.maxobj
    threshold = options.threshold


    # -- initialize variables to hold names of catalogs for each image
    tcatnamesA_global = []
    tcatnamesB_global = []


    # -- Run SExtractor on each image
    #--------------------------------
    for im in imlist:

        # -- get filter name, exposure time, zeropoint, & gain value
        fheader = pyfits.getheader(im)
        instr = fheader['INSTRUME']
        if instr == 'WFC3':
	    #  ***NOTE: WFC3 IR crclean images are in cnts, while flts are in cnts/sec
            filtname = fheader['FILTER']
	    exptime = fheader['EXPTIME']
            zeropt = get_wfc3_zeropoint(filtname) + 2.5*np.log10(exptime)
            hstgain = 1.0
        elif instr == 'ACS':
            filtname = fheader['FILTER1']
            if filtname[0] == 'C': filtname = fheader['FILTER2']
	    exptime = fheader['EXPTIME']
            zeropt = get_acs_zeropoint(fheader) + 2.5*np.log10(exptime)
            hstgain = fheader['CCDGAIN']
        else: raise Exception('Instrument '+instr+' is not supported in case list.')


        # -- assign config name
	configname = filtname.lower()+'.sex.crclean.config'


	# -- determine the science extension for each image & store chip information; assign output catalog name based on chip info
	sciext = []
	chipid = []
	catname = []
	tcatname = []
        hdulist = pyfits.open(im)
        for ff in xrange(len(hdulist)):
            if hdulist[ff].name == 'SCI':
                sciext.append(ff)
                fheader = pyfits.getheader(im,ext=sciext[-1])
		chipid.append(fheader.get('CCDCHIP',default=1))
		catname.append(im.split('_')[0] +'_sci'+str(len(chipid))+'.sex.all')
                tcatname.append(catname[-1].split('.all')[0])
		if len(sciext) == 1: tcatnamesA_global.append(tcatname[-1])
		elif len(sciext) == 2: tcatnamesB_global.append(tcatname[-1])
		else: raise Exception('Unexpected number of SCI extensions (>2).')
	if len(sciext) == 1: tcatnamesB_global.append('')


        # -- run SExtractor
	for ff in xrange(len(sciext)):

		# -- create SE catalogs & make corresponding ds9 region file with Kron apertures.
		#    ***NOTE: SExtractor does not consider the zeroth extension to be the primary extension if it doesn't contain science data.
		call(['sex','-c',configname,im+'['+str(sciext[ff]-1)+']','-CATALOG_NAME',catname[ff],'-GAIN',str(hstgain),'-MAG_ZEROPOINT',str(zeropt),'-DETECT_THRESH',str(threshold),'-ANALYSIS_THRESH',str(threshold)])
		call(['cat2reg',catname[ff]])

		# -- trim catalog of faint sources and those with bad flags (29 keeps 0 & 2). Found that APER flux inside 5pix diameter gives best tweakreg residuals.
		id,x,y,a,b,theta,kron,xwin,ywin,ra,dec,flux,fluxerr,mag,magerr,faper1,faper2,faper1err,faper2err,fwhm,flag,sclass = np.loadtxt(catname[ff],unpack=True)
		fluxmin = sorted(faper2,reverse=True)[np.min([len(faper2),maxobj])-1]
        	gd = np.where(((np.int32(flag) & 29) == 0) & (faper2 > fluxmin))[0]
		np.savetxt(tcatname[ff],zip(id[gd],x[gd],y[gd],a[gd],b[gd],theta[gd],kron[gd],xwin[gd],ywin[gd],ra[gd],dec[gd],flux[gd],fluxerr[gd],mag[gd],magerr[gd],faper1[gd], faper2[gd],faper1err[gd], faper2err[gd],flag[gd],sclass[gd]))

		# -- create ds9 region files for trimmed catalog
		call(['cat2reg',tcatname[ff]])


    # -- Create "catfile" for input to tweakreg
    #----------------------------------------
    fltlist = [imlist[ff].split('crclean')[0] + 'flt.fits' for ff in xrange(len(imlist))]
    np.savetxt('catfile.sex',zip(fltlist,tcatnamesA_global,tcatnamesB_global),fmt='%s')
