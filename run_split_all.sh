#!/bin/bash
local=false
while getopts "l" flag
do
    case "${flag}" in
        l) LOCAL=true;;
    esac
done
shift $(($OPTIND - 1))
FILFILES=$@
echo $FILFILES
SCRIPT_DIR="$(dirname $(readlink -f $0))"

for FIL in $FILFILES;
do
    PULSAR=$(echo "$FIL" | rev | cut -f2- -d '.' | rev)
    SP="${PULSAR}"_split_complete
    if [ -f $SP ]; then
        #split has run successfully
        echo $FIL has been split successfully
    else
        if [ "$LOCAL" = true ]; then
            $SCRIPT_DIR/split_files.sh -l -a $SCRIPT_DIR -p $FIL
        else
            sbatch $SCRIPT_DIR/split_files.sh -a $SCRIPT_DIR -p $FIL
        fi
    fi
done
