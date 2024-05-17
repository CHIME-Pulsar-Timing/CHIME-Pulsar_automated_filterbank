#!/bin/sh

#this script will take all the arguments and convert the fits file to filterbank file

AFP="$(dirname $(readlink -f "$0"))"
#loop over all arguments
for fitsfile in "$@"
do
    #get the base name of the fits file
    base=$(basename $fitsfile .fits)
    #convert the fits file to filterbank file
    sbatch $AFP/fits2fil.sh $fitsfile $base
done
