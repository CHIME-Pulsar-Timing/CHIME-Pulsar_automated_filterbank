#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=40:00:00
#SBATCH --mem=4GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fetch
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --gres=gpu:v100l:1
#the first argument is the tree to search down
#run FETCH

#the following code is only valid for Adam's personal computer comment out if on CC
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate fetch
#the following is valid for CC
# source ~/projects/rrg-istairs-ad/GWG2/environments/AFP/bin/activate
#work in absolute paths, CC is weird when launching batch script
while getopts "ai:" flag
do
    case "${flag}" in
        a) ADDITIONAL=true;;
        i) MY_PATH=$OPTARG;;
    esac
done
AP=$(readlink -f $MY_PATH)
#lets find all directories where we've run prep_fetch
PROCESSED=$(find $AP -name 'cands128_0.csv' -printf '%h\n' | sort -u)

for CAND_PATH in $PROCESSED;
do
    cd $CAND_PATH
    if test -f filfiles.tar.gz; then
        tar -xzf filfiles.tar.gz
    fi

    #don't do predict as we don't have GPU allocation... this can be done in seperate script
    #if we have the second argument then
    if [ "$ADDITIONAL" = true ] ; then
        #make plots and do a predict for general pulses
        PLOT=nsub_128_0/
        FP128=cands128_0.csv
        if [ ! -d $PLOT ]; then
            mkdir $PLOT
        fi
        candmaker.py --frequency_size 256 --time_size 256 --cand_param_file $FP128 --plot --fout $PLOT
        predict.py --data_dir $PLOT --model a --probability 0.1

        #do the small dm_range one for very short timescales pulses
        PLOT=nsub_128_0_short/
        FP128=cands128_0.csv
        if [ ! -d $PLOT ]; then
            mkdir $PLOT
        fi
        candmaker.py --frequency_size 256 --time_size 256 --dm_range 5 --cand_param_file $FP128 --plot --fout $PLOT
        predict.py --data_dir $PLOT --model a --probability 0.1

        #do the 1 second one for long timescales pulses
        PLOT=nsub_128_1/
        FP128=cands128_1.csv
        if [ ! -d $PLOT ]; then
            mkdir $PLOT
        fi
        candmaker.py --frequency_size 256 --time_size 256 --cand_param_file $FP128 --plot --fout $PLOT
        predict.py --data_dir $PLOT --model a --probability 0.1

    fi
    #once it has finished everything, tar all the files up
    tar -zcvf filfiles.tar.gz *.fil
    rm *.fil
    #combine the results if we have split things
    cd $AP
    echo $CAND_PATH >> combined_results.csv
    cat "$CAND_PATH/${PLOT}"results_*.csv >> combined_results.csv
done
