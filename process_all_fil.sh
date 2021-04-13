#!/bin/bash
SPLIT_SIZE=3
AFP="/home/adamdong/scratch/CHIME-Pulsar_automated_filterbank/"
DM=36
for FIL in $@;
do
    FN=$(echo "$FIL" | cut -f 1 -d '.')
    if [ ! -d $FN ]; then
    	mkdir $FN
    fi
    cp -d $FIL $FN
    cd $FN
    sbatch $AFP/automated_filterbank_batch.sh $SPLIT_SIZE $FIL $DM $AFP
    cd ..
done
