#!/bin/bash
split_size=3
afp="/home/adamdong/scratch/CHIME-Pulsar_automated_filterbank/"
dm=130
for filfile in $@;
do
    foldername=$(echo "$filfile" | cut -f 1 -d '.')
    mkdir $foldername
    cp -d $filfile $foldername
    cd $foldername
    sbatch $afp/automated_filterbank_batch.sh $split_size $filfile $dm $afp
    cd ..
done
