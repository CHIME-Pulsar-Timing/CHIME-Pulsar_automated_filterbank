#!/bin/bash
module use /project/6004902/modulefiles
module load presto
AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
#comment out on cc
# SLURM_TMPDIR="/media/adam/1c126a4b-fb16-4471-909f-4b0fda74a5d2/tmpdir"
SCRATCH_DIR=$(pwd)
for FIL in $@
do
    sbatch $AFP/fdp.sh $FIL
done
