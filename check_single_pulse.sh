#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
#set default value for prep_ts
prep_ts=0
while getopts "bfd:t:" flag
do
    case "${flag}" in
        b) RBATCH=true;;
        f) RFETCH=true;;
        d) DM=$OPTARG;;
        t) prep_ts=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))
FILFILES=$@
echo $FILFILES
BATCH=false
FETCH=false
for FIL in $FILFILES;
do
    #strip the extension
    PULSAR=$(echo "$FIL" | rev | cut -f2- -d '.' | rev)
    if [ -d $PULSAR ]; then
        SP="${PULSAR}/"*"singlepulse.ps"
        if [ -f $SP ]; then
            #now finally check if results has been run
            #Check FETCH 1 has been run
            FP="${PULSAR}/nsub_128_0/results_a.csv"
            echo $PULSAR
            if [ ! -f $FP ]; then
                #echo $FP
                #echo "$FIL never ran FETCH missing 0"
                #ls -lHd $FIL
                FETCH=true
            fi

            #check FETCH 2 has been run
            FP="${PULSAR}/nsub_128_1/results_a.csv"
            if [ ! -f $FP ]; then
                #echo $FP
                #echo "$FIL never ran FETCH missing 1"
                #ls -lHd $FIL
                FETCH=true
            fi

            #check FETCH 3 has been run
            FP="${PULSAR}/nsub_128_0_short/results_a.csv"
            if [ ! -f $FP ]; then
                #echo $FP
                #ls -lHd $FIL
                #echo "$FIL never ran FETCH missing short"
                FETCH=true
            fi

            if [ "$FETCH" = false ]; then
                echo "$FIL finished everything nothing to see here..." >> completed.csv
        else
        #check if cands is empty
        if [ -s ${PULSAR}/cands128_0.csv ]
        then
            FETCH=true
            echo "**** printing cands *****"
            cat "${PULSAR}"/*cands*.csv
            echo "****end cands*****"
            echo "$FIL never ran FETCH"
            ls -lHd $FIL
        else
            FETCH=false
            echo "${PULSAR} - cands file empty"
        fi
            fi

        else
            echo "$FIL never finished running single_pulse_search.py"
            ls -hlHd $FIL
            BATCH=true
        fi
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
            #this will send the batch job and after it's done sent the fetch job
            $SCRIPT_DIR/process_all_fil.sh -d $DM -f $FIL -t $prep_ts
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
