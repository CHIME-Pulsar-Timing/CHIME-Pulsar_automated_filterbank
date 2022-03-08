#!/bin/bash
module use /project/6004902/modulefiles
module load presto
AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
#comment out on cc
# SLURM_TMPDIR="/media/adam/1c126a4b-fb16-4471-909f-4b0fda74a5d2/tmpdir"
SCRATCH_DIR=$(pwd)
for FIL in $@
do
    #strip the fil at the end
    PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
    FDP_fn="${PULSAR}"_fdp.fil
    if [ ! -f $FDP_fn ]; then
        echo "submitting job for "${PULSAR}
        sbatch $AFP/fdp.sh $FIL
    fi
done
