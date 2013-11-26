#!/bin/bash

#Panayiotis Tzanavaris 1/2008 -- original script
#Derek Hammer 07/2013 -- modified to handle MEF extensions to images

#Run sextractor and plot ellipses

#EG CALL: sexp wmtot.fits config.sex ext

fitsfile=$1
configfile=$2
#ext=$3

sex $fitsfile -c $configfile

catfile=`gawk '{if($1~/CATALOG_NAME/) {print $2}}' $configfile`
regionfile=$catfile.reg

cat2reg $catfile

ds9 $fitsfile -region $regionfile &
