#!/bin/bash
#SBATCH --account=rrg-istairs-ad
#SBATCH --export=NONE
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fix_dropped
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#This script aims to fix the dropped packets
#load presto
#uncomment on CC
module use /project/6004902/modulefiles
module load presto
AFP="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
#comment out on cc
# SLURM_TMPDIR="/media/adam/1c126a4b-fb16-4471-909f-4b0fda74a5d2/tmpdir"
SCRATCH_DIR=$(pwd)
FIL=$1
cp $FIL ${SLURM_TMPDIR}
#go to compute node
cd ${SLURM_TMPDIR}
python $AFP/fdp.py $FIL
#come back from compute node
cd $SCRATCH_DIR
PULSAR=$(echo "$FIL" | cut -f 1 -d '.')
cp ${SLURM_TMPDIR}/"$PULSAR"_fdp.fil $SCRATCH_DIR
#clear tmpdir
rm ${SLURM_TMPDIR}/*
