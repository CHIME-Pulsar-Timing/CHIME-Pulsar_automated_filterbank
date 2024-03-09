#!/bin/bash
#SBATCH --account=rrg-istairs-ad
#SBATCH --export=NONE
#SBATCH --time=3:00:00
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --job-name=split_observation
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
set -euo pipefail
LOCAL=false
while getopts "la:p:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        a) AFP=${OPTARG};;
        p) p=${OPTARG};;
    esac
done

PULSAR=$(echo "$p" | rev | cut -f2- -d '.' | rev)
EXT="${p##*.}"
if [ "$LOCAL" != true ]; then
    source ~/util/load_presto.sh
    module load presto
    module load chime-psr
    source ~/rrat_or_not_width/bin/activate 
else
    SLURM_JOB_ID=1
fi
python $AFP/split_fits.py -fil $p
#touch a file to indicate that it has run successfully
touch "$PULSAR"_split_complete
