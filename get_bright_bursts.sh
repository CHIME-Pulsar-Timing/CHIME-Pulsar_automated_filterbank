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
mkdir -p positive_bursts_short

for RESULT_PATH in $PROCESSED;
do
    RESULT_FILE="$RESULT_PATH/results_a.csv"
    echo $RESULT_PATH
    while IFS=, read -r cand filepath probability score
    do
        ROOT=$(dirname $RESULT_PATH)
        if [ $score = "1.0" ]; then
            CAND_PATH=$ROOT/$filepath
            CAND_PATH=$(echo "$CAND_PATH" | sed 's/.h5//')
            PATH_PNG=$CAND_PATH.png

            if [[ $PATH_PNG == *"short"* ]];
            then
                # code if found
                cp $PATH_PNG positive_bursts_short
            else
                # code if not found
                cp $PATH_PNG positive_bursts
            fi

            echo "$CAND_PATH,$probability,$score" >> positive_bursts.csv
        fi
    done < $RESULT_FILE
    #once it has finished everything, tar all the files up
    # tar -zcvf $RESULT_PATH/../filfiles.tar.gz $RESULT_PATH/../*.fil
    # rm $RESULT_PATH/../*.fil
done
