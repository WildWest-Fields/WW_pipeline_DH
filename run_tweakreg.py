#! /usr/bin/env python

'''
ABOUT:
This program runs tweakreg on a list of _fl? images or drizzled-to-drizzled images using external catalogs to compare object centroids.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
July 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_tweakreg.py
'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
from pyraf import iraf
import numpy as np
from drizzlepac import tweakreg


if __name__=='__main__':
    # Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-im', '--images',default='*fl?.fits', type=str, help='Input image file(s). \
				 Default is all _fl? images in current directory.')
    parser.add_argument('-rim', '--refim',default='', type=str, help='Input refererence image. \
                                 There is no default - must be entered.')
    parser.add_argument('-rcat', '--refcat',default='', type=str, help='Input image file(s). \
                                 There is no default - must be entered.')
    options = parser.parse_args()


    # -- parse input to hold image and reference info
    im = glob.glob(options.images)
    im.sort()
    irefim = glob.glob(options.refim)
    irefim.sort()
    irefcat = glob.glob(options.refcat)
    irefcat.sort()


    if len(irefim) > 0:
        irefim = irefim[0]
        irefcat = irefcat[0]
	USE_REF = True
    else:
        irefim = ''
        irefcat = ''
        USE_REF = False


    # print image list to file
    f = open('imlist.dat', 'w')
    for ff in im: f.write(ff+'\n')
    f.close()
    catfilename = glob.glob('catfile.???')[0]


    # -- set conv width based on filter name (IR = 2.5; optical=3.5)
    fheader = pyfits.getheader(im[0])
    instr = fheader['INSTRUME']
    if instr == 'WFC3':
        filtname = fheader['FILTER']
    else: #assuming ACS
        filtname = fheader['FILTER1']
        if filtname[0] == 'C': filtname = fheader['FILTER2']  

    if filtname[0:2] == 'F1': conv_wid = 2.5
    else: conv_wid = 3.5


    # -- run tweakreg (SExtractor catalog: XWIN cols 8/9; XORG cols 2/3)
    iraf.unlearn('tweakreg')
    if USE_REF:
        tweakreg.TweakReg('@imlist.dat',catfile=catfilename,xcol=8,ycol=9,refimage=irefim,refcat=irefcat,refxyunits='pixels',refxcol=8,refycol=9,\
			conv_width=conv_wid,searchrad=1.0,updatehdr=True,nclip=7,shiftfile=True,outshifts='shift.dat', see2dplot=False,residplot='No plot')
    else:
        tweakreg.TweakReg('@imlist.dat',catfile=catfilename,xcol=8,ycol=9,conv_width=conv_wid,searchrad=1.0,\
			updatehdr=True,nclip=7,shiftfile=True,outshifts='shift.dat',see2dplot=False,residplot='No plot')


    # -- remove unwanted tweakreg files
    tmp=np.concatenate((glob.glob('*.coo'),glob.glob('*catalog.match'),glob.glob('shift*.fits')))
    for ff in tmp: os.remove(ff)

