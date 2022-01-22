#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
FILFILES=*.fil
BATCH=false
FETCH=false
for FIL in $FILFILES;
do
    #strip the extension
    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    if [ -d $PULSAR ]; then
        #lists all the splits in the pulsar's directory
        SPLITS=$PULSAR/*/
        for SPLIT in $SPLITS;
        do
            SP="$SPLIT${PULSAR}"*"singlepulse.ps"
            if [ -f $SP ]; then
                #check that the SPEG file has more than 1 line
                SPEGL=$(cat 0_SPEG_all.csv)
                if (( SPEGL > 1 )); then
                    #untar all the filterbank files
                    tar -xf $SPLIT/filfiles.tar.gz
                    #now need to do matching between the files
                fi
            fi
        done
    fi
done
