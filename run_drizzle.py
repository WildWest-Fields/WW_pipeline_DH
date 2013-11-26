#! /usr/bin/env python

'''
ABOUT:
This program drizzles a stack of flt images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
July 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_drizzle.py
'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
from pyraf import iraf
import numpy as np
from drizzlepac import astrodrizzle
from stsci.tools import teal


if __name__=='__main__':
    # Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-im', '--images',default='*fl?.fits', type=str, help='Input image file(s). \
				 Default is all _fl? images in current directory.')
    options = parser.parse_args()


    # -- initialize variables to hold filenames & print image list to file
    imlist = glob.glob(options.images)
    imlist.sort()
    f = open('imlist.dat', 'w')
    for ff in imlist: f.write(ff+'\n')
    f.close()


    # -- get instrument and filter name
    instrum = pyfits.getheader(imlist[0])['INSTRUME']
    if instrum == 'WFC3':
            filtname = pyfits.getheader(imlist[0])['FILTER']
    elif instrum == 'ACS':
            filtname = pyfits.getheader(imlist[0])['FILTER1']
            if filtname[0] == 'C': filtname = pyfits.getheader(imlist[0])['FILTER2']


    # -- run AstroDrizzle
    teal.unlearn('astrodrizzle')
    iraf.unlearn('astrodrizzle')
    #astrodrizzle.AstroDrizzle('@imlist.dat',configobj='drizzle_'+instrum+'.cfg')

    if instrum == 'WFC3':
        astrodrizzle.AstroDrizzle('@imlist.dat',output=filtname.lower(),num_cores=4,final_bits='64', in_memory=True,clean=True, combine_type='imedian',preseve=False, \
    				   final_wcs=True,final_rot=0.0,final_scale=0.05,final_pixfrac=0.4)
    elif instrum == 'ACS':
    	astrodrizzle.AstroDrizzle('@imlist.dat',output=filtname.lower(),num_cores=4,final_bits='64,32',in_memory=True,clean=True, combine_type=imedian,preserve=False,\
                                   final_wcs=True,final_rot=0.0,final_scale=0.05,final_pixfrac=0.7, final_wht_type='IVM') # weigh type dependent on avail. of rms maps
    else: raise Exception('Instrument '+instrum+' not covered in our case list.')


    # -- remove unwanted astrodrizzle files
    tmp=np.concatenate((glob.glob('*single*fits'),glob.glob('*mask.fits'),glob.glob('tmp*.fits'), glob.glob(*blt.fits')))
    for ff in tmp: os.remove(ff)

