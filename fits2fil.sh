#!/bin/bash
#SBATCH --account=rrg-istairs-ad
#SBATCH --export=NONE
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=2
#SBATCH --job-name=fix_dropped
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

source ~/util/load_presto.sh
psrfits2fil.py -o "$2".fil -n 16 $1
#downsample the data as well
downsample_filterbank.py 16 $2.fil
