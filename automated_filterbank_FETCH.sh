#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=4:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fetch
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --gres=gpu:p100:1
#the first argument is the tree to search down
#run FETCH

#the following code is only valid for Adam's personal computer comment out if on CC
#source ~/anaconda3/etc/profile.d/conda.sh    
#conda activate fetch
#the following is valid for CC
source ~/afp2/bin/activate
#work in absolute paths, CC is weird when launching batch script
AP=$(readlink -f $1)
#lets find all directories where we've run prep_fetch
PROCESSED=$(find $AP -name 'cands.csv' -printf '%h\n' | sort -u)
echo $PROCESSED
for CAND_PATH in $PROCESSED;
do
    FP=cands.csv
    DATA=data/
    cd $CAND_PATH
    if [ ! -d $DATA ]; then
        mkdir $DATA
    fi
    candmaker.py --frequency_size 256 --time_size 256 --cand_param_file $FP --plot --fout $DATA
    #don't do predict as we don't have GPU allocation... this can be done in seperate script
    predict.py --data_dir $DATA --model a
    cd $AP
    echo $CAND_PATH >> combined_results.csv
    cat "$CAND_PATH/${DATA}"results_*.csv >> combined_results.csv
done
