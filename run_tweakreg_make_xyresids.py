#! /usr/bin/env python

'''
ABOUT:
This script reads in the matched catalogs from tweakreg and recreates the four
diagnostic delta-x/y vs. x/y residual diagrams, and the single vectorgram, constructed by tweakreg in interactive mode.
This is necessary as tweakreg does not allow us to save the diagnostic diagrams in batch mode.


DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2013

HISTORY:
July 2013: Original script (v0.1).

FUTURE IMPROVEMENTS:

USE:
python run_tweakreg_make_xyresids.py
'''

__author__='D.M. HAMMER'
__version__= 0.1


import glob, argparse, pdb, pylab, pyfits
import numpy as np
from matplotlib.ticker import MultipleLocator


if __name__=='__main__':
    # -- parse input parameters
    parser = argparse.ArgumentParser(description='Read in xy residuals from tweakreg matched catalog files and plot.')
    parser.add_argument('-cat', '--catalogs',default='*fit.match', type=str, help='Input tweakreg catalogs with residuals. \
                        Default is all *fit.match files in working directory.')
    options = parser.parse_args()


    # -- save images and tweakreg catalogs in arrays
    fitcats = glob.glob(options.catalogs)
    fitcats.sort()
    images = [fitcats[i].split('_catalog')[0]+'.fits' for i in xrange(len(fitcats))]
    impref = [fitcats[i].split('_')[0] for i in xrange(len(fitcats))]


    # Iterate over each catalog and make diagrams
    # --------------------------------------------
    for cat,im,pref,ctr in zip(fitcats,images,impref,xrange(len(fitcats))):

        # -- read in matched catalog (note that residuals are fit-input)
        xref,yref,x,y,xs,ys,xreforg,yreforg,xorg,yorg,chip = np.loadtxt(cat,unpack=True,usecols=(0,1,2,3,6,7,8,9,10,11,14))


        # -- generate fit statistics
        xrms = np.std(xs)
        yrms = np.std(ys)


        # -- record filter name
        checkim = glob.glob(im)
        if len(checkim) > 0:
            fheader = pyfits.getheader(im)
            instr = fheader['INSTRUME']
            if instr == 'WFC3':
                filtname = fheader['FILTER']
            else: #assuming ACS
                filtname = fheader['FILTER1']
                if filtname[0] == 'C': filtname = fheader['FILTER2']
        else: filtname=''

	
	# -- establish x-range & xtick intervals for diagrams
	if filtname[0:2] == 'F1':
            x0 = -520
	    x1 = 520
	    xtmin = 50
	    xtmax = 200
            xtmin_vec = 50
            xtmax_vec = 200
	else:
            x0 = -2100
            x1 = 2100
	    xtmin = 100
	    xtmax = 1000
            xtmin_vec = 200
            xtmax_vec = 1000


        # -- establish y-range for diagrams
        ymax = np.max(np.abs([np.min(np.concatenate((xs,ys))),np.max(np.concatenate((xs,ys)))]))
        y0 = -1.0 * np.max([(ymax + 0.1),0.5])
        y1 = np.max([(ymax + 0.1),0.5])


	# MAKE X/Y RESIDUALS VS X/Y POSITION DIAGRAM
        # ------------------------------------------
	if ctr == 0:
        	fig1=pylab.figure(figsize=(11,8))
        	fig1.subplots_adjust(wspace=0.3, hspace=0.2)
        	fig1.suptitle(filtname+' XY Residuals in UDF ('+pref+').  RMS(X)={:1.2f}'.format(xrms)+'  RMS(Y)={:1.2f}'.format(yrms)+'  # objects='+str(len(xs)), \
                		ha='center', color='black', weight='normal')
	else:
		pylab.figure(fig1.number)
		pylab.clf()
		fig1.suptitle(filtname+' XY Residuals in UDF ('+pref+').  RMS(X)={:1.2f}'.format(xrms)+'  RMS(Y)={:1.2f}'.format(yrms)+'  # objects='+str(len(xs)), \
				ha='center', color='black', weight='normal')

	# -- delta-X vs. X
        ax1=pylab.subplot(2,2,1)
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

	pylab.scatter(xref,xs,s=4)
	pylab.axhline(0.0,color='red')
	pylab.xlim(x0,x1)
	pylab.ylim(y0,y1)
	pylab.xlabel('X (pixels)')
	pylab.ylabel('$\Delta$X (pixels)')


        # -- delta-X vs. Y
        ax1=pylab.subplot(2,2,2)
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

        pylab.scatter(yref,xs,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('Y (pixels)')
        pylab.ylabel('$\Delta$X (pixels)')


        # -- delta-Y vs. X
        ax1=pylab.subplot(2,2,3)
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))
        
        pylab.scatter(xref,ys,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('X (pixels)')
        pylab.ylabel('$\Delta$Y (pixels)')


        # -- delta-Y vs. Y
        ax1=pylab.subplot(2,2,4)
        ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_major_locator(MultipleLocator(0.2))
        ax1.xaxis.set_minor_locator(MultipleLocator(xtmin))
        ax1.xaxis.set_major_locator(MultipleLocator(xtmax))

        pylab.scatter(yref,ys,s=4)
        pylab.axhline(0.0,color='red')
        pylab.xlim(x0,x1)
        pylab.ylim(y0,y1)
        pylab.xlabel('Y (pixels)')
        pylab.ylabel('$\Delta$Y (pixels)')

	# -- save figure
	pylab.savefig(pref+'_tweakreg_xyresids.pdf')


        # MAKE VECTORGRAM
	# ----------------
        if ctr == 0: fig2=pylab.figure()
	else:
            pylab.figure(fig2.number)
            pylab.clf()

        # -- plot vector diagram
        ax2=pylab.subplot(1,1,1,aspect='equal')
        Q=ax2.quiver(xref,yref,xs,ys,color='r', angles='xy',units='xy',scale=1.0/2000.)
        pylab.xlabel('X')
        pylab.ylabel('Y')
        pylab.xlim(x0-100,x1+100)
        pylab.ylim(x0-100,x1+100)
        ax2.yaxis.set_minor_locator(MultipleLocator(xtmin_vec))
        ax2.yaxis.set_major_locator(MultipleLocator(xtmax_vec))
        ax2.xaxis.set_minor_locator(MultipleLocator(xtmin_vec))
        ax2.xaxis.set_major_locator(MultipleLocator(xtmax_vec))
        pylab.title(filtname+' XY Residuals in UDF ('+pref+').  RMS(X)={:1.2f}'.format(xrms)+'  RMS(Y)={:1.2f}'.format(yrms)+'  # objects='+str(len(xs)),size='small')
        qk = pylab.quiverkey(Q, 0.915, 0.937,0.05,'0.05 pix', labelpos='N', coordinates='axes', fontproperties={'weight': 'bold'}, color='black')
        pylab.show()
        pylab.savefig(pref+'_tweakreg_xyvector.pdf')

