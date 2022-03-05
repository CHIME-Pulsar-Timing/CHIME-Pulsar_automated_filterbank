#!/bin/bash
#this function makes the folder structures, furthermore it will submit the job
while getopts d:f: flag
do
    case "${flag}" in
        d) DM=${OPTARG};;
        f) FIL=${OPTARG};;
    esac
done

AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
FN=$(echo "$FIL" | cut -f 1 -d '.')
#make a new folder and set up to run main pipeline
if [ ! -d $FN ]; then
    mkdir $FN
fi
#copy the symbolic link into the folder we made
cp -d $FIL $FN
cd $FN
jbid_batch=$(sbatch $AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP)
# $AFP/automated_filterbank_batch.sh -d $DM -a $AFP -p $FIL
# echo $?
cd ..
jbid_fetch=$(sbatch --dependency=afterok$jbid_batch $SCRIPT_DIR/automated_filterbank_FETCH_single.sh -a -i $FN)
