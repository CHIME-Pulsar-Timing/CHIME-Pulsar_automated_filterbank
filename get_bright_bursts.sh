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
rm positive_bursts_short.csv
rm positive_bursts_1.csv

mkdir -p positive_bursts
mkdir -p positive_bursts_short
mkdir -p positive_bursts_1
mkdir -p positive_bursts_0_1
mkdir -p positive_bursts_0_1_short


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
            PATH_PNG=${CAND_PATH}.png
            PATH_H5=${CAND_PATH}.h5
            if [[ $PATH_PNG == *"nsub_short_0_5/"* ]];
            then
                # code if found
                ln -s $PATH_PNG positive_bursts_short
                echo "$CAND_PATH,$probability,$score" >> positive_bursts_short.csv

            elif [[ $PATH_PNG == *"nsub_1/"* ]];
            then
                ln -s $PATH_PNG positive_bursts_1
                echo "$CAND_PATH,$probability,$score" >> positive_bursts_1.csv

            elif [[ $PATH_PNG == *"nsub_0_5/"* ]];
            then
                # code if not found
                ln -s $PATH_PNG positive_bursts
                echo "$CAND_PATH,$probability,$score" >> positive_bursts.csv

            elif [[ $PATH_PNG == *"nsub_0_1/"* ]];
            then
                # code if not found
                ln -s $PATH_PNG positive_bursts_0_1
                echo "$CAND_PATH,$probability,$score" >> positive_bursts_0_1.csv
            elif [[ $PATH_PNG == *"nsub_0_1_short/"* ]];
            then
                # code if not found
                ln -s $PATH_PNG positive_bursts_0_1_short
                echo "$CAND_PATH,$probability,$score" >> positive_bursts_0_1_short.csv
            fi

        fi
    done < $RESULT_FILE
    #once it has finished everything, tar all the files up
    # tar -zcvf $RESULT_PATH/../filfiles.tar.gz $RESULT_PATH/../*.fil
    # rm $RESULT_PATH/../*.fil
done
