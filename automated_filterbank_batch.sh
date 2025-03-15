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
FIL=$p
echo $EXT
if [ $EXT == "fits" ]; then
    python $AFP/gwg_cand_search_pipeline.py --dm $DM --rfifind --sk_mask --kc_iqrm --dedisp --fil $FIL
fi

PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
echo "ALL DONE"
exit 0
