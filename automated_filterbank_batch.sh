#!/bin/bash
while getopts "d:p:" flag
do
    case "${flag}" in
        d) DM=${OPTARG};;
        p) p=${OPTARG};;
    esac
done

AFP="$(dirname $(readlink -f $0))"
PULSAR=$(echo "$p" | rev | cut -f2- -d '.' | rev)
EXT="${p##*.}"
#source the pulsar environtment
source /home/pulsar_rhel8/pulsar.bash
#make sure that $p is a file
if test -f "$p"; then
    #rename it FIL
    FIL=$p
    n=0
    #basically try catch
    until [ "$n" -ge 1 ]
    do
        echo $EXT
        if [ $EXT == "fits" ]; then
            python $AFP/gwg_cand_search_pipeline.py --dm $DM --rfifind --sk_mask --kc_iqrm --dedisp --fil $FIL && break
        else

        #ANYTHING BEYOND HERE IS DEBUGGING
        #FETCH and SPEGID are slow so lets like ignore that for now
        n=$((n+1))
        echo "ERRORS IN PROCESSING CHECK LOGS OR ENABLE DEBUGGING"
        exit 1
        #remove the extra fil files
    done
    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    echo "ALL DONE"
    exit 0
fi
#didn't find file, throw error
echo "filterbank file not found"
exit 1
