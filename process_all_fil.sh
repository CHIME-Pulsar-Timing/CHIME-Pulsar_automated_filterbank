#!/bin/bash
split_size=8
for filfile in $@;
do
    foldername=$(echo "$filfile" | cut -f 1 -d '.')
    mkdir $foldername
    mv $filfile $foldername
    cd $foldername
    $(dirname "$BASH_SOURCE")/split_files.sh $split_size $filfile
    cd ..
done
