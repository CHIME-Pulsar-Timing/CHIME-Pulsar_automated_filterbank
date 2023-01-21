#!/bin/bash
#this function makes the folder structures, furthermore it will submit the job
while getopts d:f:t: flag
do
    case "${flag}" in
        d) DM=${OPTARG};;
        f) FIL=${OPTARG};;
    esac
done

AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
FN=$(echo "$FIL" | rev | cut -f2- -d '.' | rev)
#make a new folder and set up to run main pipeline
if [ ! -d $FN ]; then
    mkdir $FN
fi
#copy the symbolic link into the folder we made
cp -d $FIL $FN
cd $FN
#jbid_batch=$(sbatch $AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP)
$AFP/automated_filterbank_batch.sh -d $DM -p $FIL -a $AFP

#batch job submit string
# jbid_batch=${jbid_batch#*job }
# cd ..
# sleep 1
# echo "sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -a -i $FN"
# jbid_fetch=$(sbatch --dependency=afterok:$jbid_batch $AFP/automated_filterbank_FETCH_single.sh -a -i $FN)
# jbid_fetch=${jbid_fetch#*job }
