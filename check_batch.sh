#!/bin/bash
#this file will run to check which filterbank files have been run and which have not
FILFILES=*.fil
BATCH=false
FETCH=false
SCRIPT_DIR="$(dirname $(readlink -f $0))"
for FIL in $FILFILES;
do
    #strip the extension
    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    if [ -d $PULSAR ]; then
        #lists all the splits in the pulsar's directory
        SPLITS=$PULSAR/*/
        for SPLIT in $SPLITS;
        do
            #just make the file for now
            # touch "$SPLIT${PULSAR}_singlepulse.ps"
            SP="$SPLIT${PULSAR}"*"singlepulse.ps"
            if [ -f $SP ]; then
                #check that the SPEG file has more than 1 line
                SPEGL=$(cat $SPLIT/0_SPEG_all.csv | wc -l)
                if (( SPEGL > 1 )); then
                    #untar all the filterbank files
                    tar -xf $SPLIT/filfiles.tar.gz --directory $SPLIT
                    python $SCRIPT_DIR/check_batch.py $SPLIT/0_SPEG.csv $1 $SPLIT/*.fil
                    prev=$?
                    if (( $prev > 0 )); then
                        echo $SPLIT
                        #rm $SP
                    fi
                    rm $SPLIT/*.fil
                    #now need to do matching between the files
                fi
            fi
        done
    fi
done
