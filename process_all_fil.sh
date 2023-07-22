#!/bin/bash
#this function makes the folder structures, furthermore it will submit the job
while getopts "lad:f:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        a) ALL=true;;
        d) DM=${OPTARG};;
        f) FIL=${OPTARG};;
    esac
done

AFP="$(dirname $(readlink -f $0))"
FN=$(echo "$FIL" | rev | cut -f2- -d '.' | rev)
#make a new folder and set up to run main pipeline
if [ ! -d $FN ]; then
    mkdir $FN
fi
#copy the symbolic link into the folder we made
cp -d $FIL $FN
cd $FN
if [ "$LOCAL" = true ]; then
    $AFP/automated_filterbank_batch.sh -l -d $DM -p $FIL -a $AFP &
else
    echo "$AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP"
    jbid_batch=$(sbatch $AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP)
    #batch job submit string
    jbid_batch=${jbid_batch#*job }
    sleep 1
    #sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR
    echo "sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP"
    if [ "$ALL" = true ]; then
        jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP -t 0.1)
        jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP -t 0.5)
        jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP -t 0.1 -s)
        jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP -t 0.5 -s)
    fi
    jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP -t 1)
    jbid_fetch=${jbid_fetch#*job }
fi
