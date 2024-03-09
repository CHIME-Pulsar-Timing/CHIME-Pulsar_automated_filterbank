#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
ALL=false
while getopts "bflad:" flag
do
    case "${flag}" in
        b) RBATCH=true;;
        f) RFETCH=true;;
        l) LOCAL=true;;
        a) ALL=true;;
        d) DM=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))
FILFILES=$@
BATCH=false
FETCH=false
FETCH_0_5=false
FETCH_S_0_5=false
FETCH_0_1=false
FETCH_S_0_1=false
FETCH_1=false
for FIL in $FILFILES;
do
    #strip the extension
    PULSAR=$(echo "$FIL" | rev | cut -f2- -d '.' | rev)
    if [ -d $PULSAR ]; then
        SP="${PULSAR}/"*"cands.csv"
        if [ -f $SP ]; then
            #now finally check if results has been run
            if [ "$ALL" = true ]; then
                FP="${PULSAR}/nsub_0_5/results_a.csv"
                if [ ! -f $FP ]; then
                    #echo $FP
                    #echo "$FIL never ran FETCH missing 0"
                    #ls -lHd $FIL
                    FETCH_0_5=true
                    FETCH=true
                fi
                FP="${PULSAR}/nsub_short_0_5/results_a.csv"
                if [ ! -f $FP ]; then
                    #echo $FP
                    echo "$FIL never ran FETCH missing short 0 5"
                    #ls -lHd $FIL
                    FETCH_S_0_5=true
                    FETCH=true
                fi

                FP="${PULSAR}/nsub_0_1/results_a.csv"
                if [ ! -f $FP ]; then
                    #echo $FP
                    #echo "$FIL never ran FETCH missing 1"
                    #ls -lHd $FIL
                    FETCH_0_1=true
                    FETCH=true
                fi

                FP="${PULSAR}/nsub_0_1_short/results_a.csv"
                if [ ! -f $FP ]; then
                    #echo $FP
                    #echo "$FIL never ran FETCH missing 1"
                    #ls -lHd $FIL
                    FETCH_S_0_1=true
                    FETCH=true
                fi
            fi
            FP="${PULSAR}/nsub_1/results_a.csv"
            if [ ! -f $FP ]; then
                #echo $FP
                #echo "$FIL never ran FETCH missing 1"
                #ls -lHd $FIL
                FETCH_1=true
                FETCH=true
            fi

            if [ "$FETCH" = false ]; then
                echo "$FIL finished everything nothing to see here..."
            else
                #check if cands is empty
                LINES=$(cat "$PULSAR"/cands.csv | wc -l)
                if [ "$LINES" -gt 1 ]
                then
                    FETCH=true
                    #echo "**** printing cands *****"
                    #cat "${PULSAR}"/*cands*.csv
                    #echo "****end cands*****"
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
            SCRIPT_DIR="$(dirname $(readlink -f $0))"
            #this will send the batch job and after it's done sent the fetch job
            if [ "$LOCAL" = true ]; then
                $SCRIPT_DIR/process_all_fil.sh -l -d $DM -f $FIL
            else
                if [ "$ALL" = true ]; then
                    $SCRIPT_DIR/process_all_fil.sh -d $DM -f $FIL -a
                else
                    $SCRIPT_DIR/process_all_fil.sh -d $DM -f $FIL
                fi
            fi
        fi
    fi
    if [ "$RFETCH" = true ]; then
        if [ "$FETCH" = true ]; then
            echo "submitting FETCH job for $PULSAR"
            #find the directory that the script belongs to
            SCRIPT_DIR="$(dirname $(readlink -f $0))"
            AP=$(readlink -f $PULSAR)
            #lets find all directories where we've run prep_fetch
            PROCESSED=$(find $AP -name 'cands.csv' -printf '%h\n' | sort -u)
            cd $PROCESSED

            if [ "$LOCAL" = true ]; then
                $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -l -i $PROCESSED -p $SCRIPT_DIR -t 1 -g 0 -n 20
                # $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -l -i $PROCESSED -p $SCRIPT_DIR -t 0.5 -g 0 -n 5
                # $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -l -i $PROCESSED -p $SCRIPT_DIR -t 0.5 -s -g 0 -n 5
                # $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -l -i $PROCESSED -p $SCRIPT_DIR -t 0.1 -g 0 -n 5
                # $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -l -i $PROCESSED -p $SCRIPT_DIR -t 0.1 -g 0 -s -n 5
            else
                if [ "$FETCH_0_5" = true ]; then
            echo $FETCH_0_5
            echo submitting job for FETCH 0.5
                    sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR -t 0.5
        fi
                if [ "$FETCH_S_0_5" = true ]; then
            echo $FETCH_S_0_5
            echo submitting job for FETCH S 0.5
                    sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR -t 0.5 -s
        fi
                if [ "$FETCH_0_1" = true ]; then
            echo $FETCH_0_1
            echo submitting job for FETCH 0.1
                    sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR -t 0.1
        fi
                if [ "$FETCH_S_0_1" = true ]; then
            echo $FETCH_S_0_1
            echo submitting job for FETCH S 0.1
                    sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR -t 0.1 -s
        fi
                if [ "$FETCH_1" = true ]; then
            echo $FETCH_1
            echo submitting job for FETCH 1
                    sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR -t 1
                fi
            fi
            cd ..
        fi
    fi
    BATCH=false
    FETCH=false
    FETCH_0_5=false
    FETCH_S_0_5=false
    FETCH_0_1=false
    FETCH_S_0_1=false
    FETCH_1=false
done
