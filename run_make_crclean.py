#! /usr/bin/env python

'''
ABOUT:
This programs constructs crclean images from a list of ASN files or flat-fielded HST images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER, 2013

HISTORY:
August 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_make_crclean.py		--> uses *asn.fits files to make crcleans for each group of flt (DEFAULT).
python run_make_crclean.py -na      	--> makes crcleans using all flts in current directory in one call.

'''

__author__='D.M. HAMMER'
__version__= 0.1

import os, glob, argparse, pyfits, pdb
from pyraf import iraf
import numpy as np
from drizzlepac import astrodrizzle


if __name__=='__main__':
    # -- Parse input parameters
    parser = argparse.ArgumentParser(description='Run tweakreg on input images using custom source catalogs.')
    parser.add_argument('-a', '--asn',default='*asn.fits', type=str, help='Input association fits file(s). \
                                 Default is all _asn files in current directory.')
    parser.add_argument('-na', '--no_asn',default=False, action='store_true', help='Should we construct crcleans using all flts in current directory in one call? \
                                 Default is False, i.e., by default, we construct crleans spearately for flts in each _asn.fits file.')
    parser.add_argument('-c', '--cores',default=4, type=int, help='Input number of processing cores used by AD. \
                                 Default is 4.')
    parser.add_argument('-i', '--images',default='*fl?.fits', type=str, help='Input fl? fits images. \
                                 Default is all _fl? images in current directory.')


    options = parser.parse_args()
    USE_ASN = not(options.no_asn)
    NCORES = options.cores
    asnlist = glob.glob(options.asn)
    asnlist.sort()
    imlist = glob.glob(options.images)
    imlist.sort()


    # -- construct crclean images
    if USE_ASN:
	# -- make list of assoc files
	if len(asnlist) == 0: raise Exception('No asn files located. Use switch "-na" if not using asn files.') 
	for asn in asnlist: astrodrizzle.AstroDrizzle(asn,driz_combine=False,num_cores=NCORES,driz_sep_bits='256,64,32',driz_cr_corr=True,preserve=False)
    else:
	astrodrizzle.AstroDrizzle('*fl?.fits',driz_combine=False,num_cores=NCORES,driz_sep_bits='256,64,32',driz_cr_corr=True,preserve=False)

    # -- remove unwanted files
    bad = np.concatenate((glob.glob('*single*fits'),glob.glob('*crmask.fits'),glob.glob('tmp*.fits'),glob.glob('*blt.fits'), glob.glob('*med.fits')))
    for tmp in bad: os.remove(tmp)
