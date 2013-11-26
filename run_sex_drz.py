#! /usr/bin/env python

'''
ABOUT:
This script runs SExtractor on a list of drizzled images (no default - must input filenames).

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2013

HISTORY:
Jul 2013: Original script (v0.1).
Aug 2013: Modified to run on any instrument's image (e.g., double-chip ACS/optical or single-chip WFC/IR).


FUTURE IMPROVEMENTS:
Currently, we select "maxobj" brightest objects, then keep those with good flags/stellarity/fwhm. 
**Edit code to keep max obj brightest after selecting on flags/stellarity/fwhm.


USE:
python run_sex_drz.py
'''

__author__='D.M. HAMMER'
__version__= 0.2


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
    parser.add_argument('-im', '--images',default='NONE', type=str, help='Input drizzled fits image(s). \
                        Deafult is NONE - must input image(s).')
    parser.add_argument('-ivm', '--ivmname',default='NONE', type=str, help='Input IVM weight map image(s). \
                        Default is NONE. Number of weight maps must match number of drizzled images.')
    parser.add_argument('-exp', '--expname',default='NONE', type=str, help='Input effective exposure map. \
                        Default is NONE. Use either ivm or exp (not both). Number of exposure maps must match number of drizzled images.')
    parser.add_argument('-th', '--threshold',default=10, type=int, help='Input threshold used for BOTH detect_thresh and analysis_thresh. \
                        Default is 10 - tests showed better tweakreg residuals when using high thresh for even bright objects.')
    parser.add_argument('-max', '--maxobj',default=100, type=int, help='Input maximum number of SExtractor detections to use for tweakreg matching (each chip). \
                        Default is 100 objects per chip.')
    parser.add_argument('-cs', '--stellarity',default=0.0, type=float, help='Input lower limit for selecting objects via SExtractor "CLASS_STAR." \
                        Default is 0.0, i.e., does not select on stellarity.')
    parser.add_argument('-f', '--fwhm',default=100.0, type=float, help='Input upper limit for selecting objects via SExtractor "FWHM_IMAGE." \
                        Default is 100, i.e., does not select on FWHM.')
    options = parser.parse_args()
    

    # -- initialize image list
    if options.images == 'NONE': raise Exception('Must input drizzled images.')
    imlist = glob.glob(options.images)
    imlist.sort()

   # -- initialize weight image
    if options.ivmname != 'NONE':
        whtlist = glob.glob(options.ivmname)
        whtlist.sort()
	if len(whtlist) != len(imlist): raise Exception('Number of weight images must equal number of drizzled images.')
	if options.expname != 'NONE': raise Exception('Must input either ivm maps or exposure maps, not both.')
	weighttype ='MAP_RMS'
    elif options.expname != 'NONE':
        whtlist = glob.glob(options.expname)
        whtlist.sort()
        if len(whtlist) != len(imlist): raise Exception('Number of weight images must equal number of drizzled images.')
	weighttype ='MAP_WEIGHT'
    else:
        weighttype ='NONE'
	whtlist = []
	for ff in xrange(len(imlist)): whtlist[ff] = 'NONE'

    threshold = options.threshold
    maxobj = options.maxobj
    minstellar = options.stellarity
    maxfwhm = options.fwhm


    tcatname_all = []
  
    # -- Run SExtractor on each image
    #--------------------------------
    for im,wht in zip(imlist,whtlist):

	# -- assign various image properties (gain,exposure,zeropt,filter name, pixel scale)
        fheader = pyfits.getheader(im)
        instr = fheader['INSTRUME']
        if instr == 'WFC3':
            filtname = fheader['FILTER']
	    exptime = fheader['EXPTIME']
	    zeropt = get_wfc3_zeropoint(filtname)
	    hstgain = 1.0
	    gain = hstgain * exptime
	    saturation = 1000000.	# assign large value - IR doesn't saturate
	    pscale = fheader['D001SCAL']
	    configname = 'wfc3.sex.drz.config'
        elif instr == 'ACS':
            filtname = fheader['FILTER1']
            if filtname[0] == 'C': filtname = fheader['FILTER2']
            exptime = fheader['EXPTIME']
	    zeropt = get_acs_zeropoint(fheader)
	    hstgain = fheader['CCDGAIN']
	    gain = hstgain * exptime
	    # use median exposure time to estimate saturation level
	    ftable = pyfits.getdata(im,ext=1)
	    medexp = np.median(ftable['EXPTIME'])
	    saturation = 80000./medexp
	    pscale = fheader['D001SCAL']
	    configname = 'acs.sex.drz.config'
	else: raise Exception('Instrument '+instr+' is not supported in case list.')

	# -- assign output catalog name
	catname = im.split('fits')[0]+'sex.all'
	tcatname = im.split('fits')[0]+'sex'
	tcatname_all.append(tcatname)

	# -- construct rms map if "IVM" weight map was requested
	if weighttype == 'WEIGHT_RMS':
	    shutil.copy(wht,'./rms.fits')
	    rms = pyfits.open('rms.fits',mode='update')
	    data = rms[0].data
	    rms[0].data = 1.0/np.sqrt(data)
	    rms.close()
 	    wht = 'rms.fits'

        # -- run SExtractor
        call(['sex','-c',configname,im,'-WEIGHT_IMAGE',wht+','+wht,'-WEIGHT_TYPE',weighttype+','+weighttype,'-CATALOG_NAME',catname, \
	      '-PIXEL_SCALE',str(pscale),'-GAIN',str(gain),'-MAG_ZEROPOINT',str(zeropt),'-SATUR_LEVEL',str(saturation), \
	      '-DETECT_THRESH',str(threshold),'-ANALYSIS_THRESH',str(threshold)])
	call(['cat2reg',catname])

	# -- trim catalog of faint extended sources, objects with bad flags (29 keeps 0 & 2), AND OPTIONALLY select on size/profile via fwhm and class_star.
	id,x,y,a,b,theta,kron,xwin,ywin,ra,dec,flux,fluxerr,mag,magerr,faper1,faper2,faper1err,faper2err,fwhm,flag,sclass = np.loadtxt(catname,unpack=True)
	fluxmin = sorted(faper2,reverse=True)[np.min([len(faper2),maxobj])-1]
       	gd = np.where(((np.int32(flag) & 29) == 0) & (faper2 > fluxmin) & (sclass > minstellar) & (fwhm < maxfwhm))[0]	# USING APER FLUX
	np.savetxt(tcatname,zip(id[gd],x[gd],y[gd],a[gd],b[gd],theta[gd],kron[gd],xwin[gd],ywin[gd],ra[gd],dec[gd],flux[gd],fluxerr[gd],mag[gd],magerr[gd],faper1[gd], faper2[gd],faper1err[gd], faper2err[gd],flag[gd],sclass[gd]))

	# -- create ds9 region files for trimmed catalog
       	call(['cat2reg',tcatname])


    # -- Create "catfile" for input to run_tweakreg_drz
    #--------------------------------------------------
    np.savetxt('catfile.sex',zip(imlist,tcatname_all),fmt='%s')
