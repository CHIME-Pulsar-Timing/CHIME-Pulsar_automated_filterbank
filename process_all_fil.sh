#!/bin/bash
i=0
#first argument is split size
#second argument is DM
SPLIT_SIZE=$1
AFP="/home/adamdong/scratch/CHIME-Pulsar_automated_filterbank/"
#AFP="/home/adam/Documents/automated_filterbank/"
DM=$2
if [ "$#" -le 2 ]; then
    echo "not enough arguments, need 3 at least first is split size, second is DM, rest are filterbank files"
    exit 1
fi
for FIL in $@;
do
    if [ $i -gt 1 ]; then
        FN=$(echo "$FIL" | cut -f 1 -d '.')
        if [ ! -d $FN ]; then
            mkdir $FN
        fi
        cp -d $FIL $FN
        cd $FN
        sbatch $AFP/automated_filterbank_batch.sh $SPLIT_SIZE $DM $FIL $AFP
        #$AFP/automated_filterbank_batch.sh $SPLIT_SIZE $DM $FIL $AFP
        cd ..	
    fi
    i=$((i+1))
done
