#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=1:00:00
#SBATCH --mem=4096M
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fetch
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --gres=gpu:v100l:1
#the first argument is the tree to search down
#run FETCH

LOCAL=false
GPU=0
n=5
TIME=1
SHORT=false
while getopts "li:p:g:n:t:s" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        i) CAND_PATH=$OPTARG;;
        p) AFP=$OPTARG;;
        g) GPU=$OPTARG;;
        n) n=$OPTARG;;
        t) TIME=$OPTARG;;
        s) SHORT=true;;
    esac
done
echo $GPU

if [ "$LOCAL" != true ]; then
    #module use /project/6004902/chimepsr-software/v1/environment-modules
    module use /project/6004902/chimepsr-software/v2/environment-modules
    #module use /project/6004902/chimepsr-software/v1/environment-modules
    #module load presto
    module load cuda
    module load chime-psr
    source ~/projects/rrg-istairs-ad/CHIPSPIPE_CANDMAKER/bin/activate
    #module use /project/6004902/chimepsr-software/v1/environment-modules
    #source ~/projects/rrg-istairs-ad/adamdong/Your_161123/bin/activate
fi

#work in absolute paths, CC is weird when launching batch script
CAND_PATH=$(pwd)
echo $PWD $CAND_PATH
if [ "$LOCAL" != true ]; then
    ls
    cp -r ./* $SLURM_TMPDIR
    cd $SLURM_TMPDIR
fi
# Set the maximum number of attempts
max_attempts=5

# Set a counter for the number of attempts
attempt_num=1

# Set a flag to indicate whether the command was successful
success=false

# Loop until the command is successful or the maximum number of attempts is reached
while [ $success = false ] && [ $attempt_num -le $max_attempts ]; do
    if [ "$TIME" == 0.1 ]; then
        #see if $SHORT is true
        if [ "$SHORT" == true ]; then
            mkdir -p nsub_short_0_1
            echo "short 0.1"
            python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_short_0_1 -r -n $n -ws 100 --gpu_id $GPU --range_dm 5 --no_log_file
            # Check the exit code of the command
            if [ $? -eq 0 ]; then
                # The command was successful
                success=true
            else
                # The command was not successful
                echo "Attempt $attempt_num failed. Trying again..."
                # Increment the attempt counter
                attempt_num=$(( attempt_num + 1 ))
            fi
        else
            mkdir -p nsub_0_1
            echo "not short 0.1"
            python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_1 -r -n $n -ws 100 --gpu_id $GPU --no_log_file
            # Check the exit code of the command
            if [ $? -eq 0 ]; then
                # The command was successful
                success=true
            else
                # The command was not successful
                echo "Attempt $attempt_num failed. Trying again..."
                # Increment the attempt counter
                attempt_num=$(( attempt_num + 1 ))
            fi
        fi

    elif [ "$TIME" == 0.5 ]; then
        if [ "$SHORT" == true ]; then
            mkdir -p nsub_short_0_5
            echo "short 0.5"
            echo here 0_5
            python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_short_0_5 -r -n $n -ws 500 --gpu_id $GPU --range_dm 5 --no_log_file
            # Check the exit code of the command
            if [ $? -eq 0 ]; then
                # The command was successful
                success=true
            else
                # The command was not successful
                echo "Attempt $attempt_num failed. Trying again..."
                # Increment the attempt counter
                attempt_num=$(( attempt_num + 1 ))
            fi
        else
            mkdir -p nsub_0_5
            echo "not short 0.5"
            python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_5 -r -n $n -ws 500 --gpu_id $GPU --no_log_file
            # Check the exit code of the command
            if [ $? -eq 0 ]; then
                # The command was successful
                success=true
            else
                # The command was not successful
                echo "Attempt $attempt_num failed. Trying again..."
                # Increment the attempt counter
                attempt_num=$(( attempt_num + 1 ))
            fi
        fi

    elif [ "$TIME" == 1 ]; then
        mkdir -p nsub_1
        python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_1 -r -n $n -ws 1000 --gpu_id $GPU --no_log_file
        # Check the exit code of the command
        if [ $? -eq 0 ]; then
            # The command was successful
            success=true
        else
            # The command was not successful
            echo "Attempt $attempt_num failed. Trying again..."
            # Increment the attempt counter
            attempt_num=$(( attempt_num + 1 ))
        fi
    fi
done
# Check if the command was successful
if [ $success = true ]; then
  # The command was successful
  echo "The command was successful after $attempt_num attempts."
else
  # The command was not successful
  echo "The command failed after $max_attempts attempts."
  ls
  echo "nsub*/*"
  ls nsub*/*
  exit 1
fi
echo "candidates made"


#make plots and do a predict for general pulses
echo "grading candidates"
if [ "$LOCAL" != true ]; then
    #go backwards, deactivate everything
    #
    module unload chime-psr
    module unload cuda
    deactivate
    module unuse /project/6004902/chimepsr-software/v2/environment-modules

    #module unload chime-psr
    #module unload presto
    #module use /project/6004902/chimepsr-software/v1/environment-modules
    source ~/projects/rrg-istairs-ad/CHIPSPIPE_FETCH/bin/activate
else
    echo "running locally, make sure FETCH, Your and sigpyproc are installed!!"
    source ~/anaconda3/etc/profile.d/conda.sh
    echo "activating fetch"
    conda activate fetch
    echo "Fetch activated"
fi
echo "predicting candidates"
if [ "$TIME" == 0.1 ]; then
    #see if $SHORT is true
    if [ "$SHORT" == true ]; then
        predict.py --data_dir nsub_short_0_1 --model a --probability 0.1 -g $GPU
    else
        predict.py --data_dir nsub_0_1 --model a --probability 0.1 -g $GPU
    fi
elif [ "$TIME" == 0.5 ]; then
    if [ "$SHORT" == true ]; then
        predict.py --data_dir nsub_short_0_5 --model a --probability 0.1 -g $GPU
    else
        predict.py --data_dir nsub_0_5 --model a --probability 0.1 -g $GPU
    fi
elif [ "$TIME" == 1 ]; then
    predict.py --data_dir nsub_1 --model a --probability 0.1 -g $GPU
fi

#remove all the files that are irrelevant
#$AFP/clear_false_candidates.sh -i .
#chown -R adamdong:rrg-istairs-ad $SLURM_TMPDIR/nsub*
if [ "$LOCAL" != true ]; then
    cp -r $SLURM_TMPDIR/nsub* $CAND_PATH/
fi
