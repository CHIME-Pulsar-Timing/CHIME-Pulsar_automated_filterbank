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

for RESULT_PATH in $PROCESSED;
do
    #only really care about the pngs
    rm $RESULT_PATH/*.h5
    RESULT_FILE="$RESULT_PATH/results_a.csv"
    echo $RESULT_PATH
    while IFS=, read -r cand filepath probability score
    do
        if [ $filepath != "candidate" ]; then
            ROOT=$(dirname $RESULT_PATH)
            if [ $score != "1.0" ]; then
                CAND_PATH=$ROOT/$filepath
                CAND_PATH=$(echo "$CAND_PATH" | sed 's/.h5//')
                PATH_PNG=${CAND_PATH}.png
                if [ -f $PATH_PNG ]; then
                    rm $PATH_PNG
                fi

            fi
        fi
    done < $RESULT_FILE
    #once it has finished everything, tar all the files up
    # tar -zcvf $RESULT_PATH/../filfiles.tar.gz $RESULT_PATH/../*.fil
    # rm $RESULT_PATH/../*.fil
done
