astrodrizzle.AstroDrizzle('*flt.fits',driz_combine=False,num_cores=4,driz_sep_bits='256,64,32',clean=False,driz_cr_corr=True)
bad = np.concatenate((glob.glob(*mask*fits),glob.glob('*single*.fits',glob.glob('tmp*.fits'),glob.glob('*blt*fits')))
for tmp in bad: os.remove(tmp)

tweakreg.TweakReg('*crclean.fits',conv_width=2.5,searchrad=1.0, updatehdr=False,nclip=7,shiftfile=True,outshifts='shift.dat',threshold=28,see2dplot=False,residplot='No plot')
bad = np.concatenate((glob.glob('*coo'),glob.glob('*list'), glob.glob('shift*fits'),glob.glob('*.log'),glob.glob('*catalog.match')))
for tmp in bad: os.remove(tmp)

astrodrizzle.AstroDrizzle('*crclean.fits',driz_combine=False,num_cores=4,driz_sep_bits='256,64,32',clean=False,driz_cr_corr=True)

