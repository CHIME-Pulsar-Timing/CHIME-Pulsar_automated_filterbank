#!/bin/bash
#this function makes the folder structures, furthermore it will submit the job
while getopts "ld:f:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
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
    jbid_batch=$(sbatch $AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP)
    #batch job submit string
    jbid_batch=${jbid_batch#*job }
    sleep 1
    #sbatch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -i $PROCESSED -p $SCRIPT_DIR

    echo "sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP"
    jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -i $FN -p $AFP)
    jbid_fetch=${jbid_fetch#*job }
fi
