#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
while getopts "bfd:" flag
do
    case "${flag}" in
        b) RBATCH=true;;
        f) RFETCH=true;;
        d) DM=$OPTARG;;
    esac
done

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
                #now finally check if results has been run
                FP="${SPLIT}nsub_128/results_a.csv"
                if [ ! -f $FP ]; then
                    echo $FP
                    echo "$FIL never ran FETCH"
                    FETCH=true
                fi
            else
                echo "$FIL never finished running single_pulse_search.py"
                BATCH=true
            fi
        done
    else
        echo "$FIL has no directory"
        BATCH=true
    fi
    #run the a batch for this pulsar, should probably use the base job script... but it's easier to do it this way
    if [ "$RBATCH" = true ]; then
	    if [ "$BATCH" = true ]; then
            echo "submitting batch job for $PULSAR"
            #find the directory that the script belongs to
            SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
            $SCRIPT_DIR/process_all_fil.sh 1 $DM $FIL
	    fi
    fi
    if [ "$RFETCH" = true ]; then
	    if [ "$FETCH" = true ]; then
            echo "submitting FETCH job for $PULSAR"
            #find the directory that the script belongs to
            SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
            $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -a -i $PULSAR
	    fi
    fi

    BATCH=false
    FETCH=false
done
