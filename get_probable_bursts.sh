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
rm positive_bursts.csv
mkdir -p positive_bursts
for RESULT_PATH in $PROCESSED;
do
    RESULT_FILE="$RESULT_PATH/results_a.csv"
    while IFS=, read -r cand filepath probability score
    do
        ROOT=$(dirname $RESULT_PATH)
        echo $ROOT
        if [[ "$probability" != *"e"* ]]; then
            CAND_PATH=$ROOT/$filepath
            CAND_PATH=$(echo "$CAND_PATH" | sed 's/.h5//')
            PATH_PNG=$CAND_PATH.png
            cp $PATH_PNG positive_bursts
            echo "$CAND_PATH,$probability,$score" >> positive_bursts.csv
        fi
    done < $RESULT_FILE
done
