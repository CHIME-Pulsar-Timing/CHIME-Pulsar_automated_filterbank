#!/bin/bash
while getopts "i:m:" flag
do
    case "${flag}" in
        i) MY_PATH=$OPTARG;;
	m) MIN=$OPTARG;;
    esac
done
AP=$(readlink -f $MY_PATH)
#lets find all directories where we've run prep_fetch
PROCESSED=$(find $AP -name 'results_a.csv' -printf '%h\n' | sort -u)
PROCESSED=$(readlink -f $PROCESSED)
rm probable_bursts.csv
mkdir -p probable_bursts
for RESULT_PATH in $PROCESSED;
do
    RESULT_FILE="$RESULT_PATH/results_a.csv"
    while IFS=, read -r cand filepath probability score
    do
        if [[ "$probability" != "probability" ]]; then
            ROOT=$(dirname $RESULT_PATH)
            if [[ "$probability" != *"e"* ]]; then
                if [ 1 -eq "$(echo "${probability} > ${MIN}" | bc)" ]; then
                    CAND_PATH=$ROOT/$filepath
                    CAND_PATH=$(echo "$CAND_PATH" | sed 's/.h5//')
                    PATH_PNG=$CAND_PATH.png
                    cp $PATH_PNG probable_bursts
		    echo "copying path:"
	   	    echo $PATH_PNG
	            echo $probability
                    echo "$CAND_PATH,$probability,$score" >> probable_bursts.csv
                fi
            fi
        fi
    done < $RESULT_FILE
done
