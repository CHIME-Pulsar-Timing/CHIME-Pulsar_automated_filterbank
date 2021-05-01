#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
FILFILES=*.fil
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
                #now finally check if results has been run
                FP="${SPLIT}data/results_a.csv"

                if [ ! -f $FP ]; then
                    echo $FP
                    echo "$FIL never ran FETCH"
                fi
            else
                echo "$FIL never finished running single_pulse_search.py"
            fi
        done
    else
        echo "$FIL has no directory"
    fi

done
