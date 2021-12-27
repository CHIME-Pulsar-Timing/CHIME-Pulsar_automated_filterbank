#!/usr/bin/env bash
#!/bin/bash
while getopts "i:" flag
do
    case "${flag}" in
        i) MY_PATH=$OPTARG;;
    esac
done
AP=$(readlink -f $MY_PATH)
#lets find all directories where we've run prep_fetch
PROCESSED=$(find $AP -name 'results_a.csv' -printf '%h\n' | sort -u)
PROCESSED=$(readlink -f $PROCESSED)
mkdir -p $MY_PATH/"all_bursts"
for RESULT_PATH in $PROCESSED;
do
    cp $RESULT_PATH/*.png all_bursts
done
