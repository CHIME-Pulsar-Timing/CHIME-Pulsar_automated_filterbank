#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=16:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=5
#SBATCH --job-name=fetch
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --gres=gpu:v100l:1
#the first argument is the tree to search down
#run FETCH

LOCAL=false
GPU=0
n=5

while getopts "li:p:g:n:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        i) CAND_PATH=$OPTARG;;
        p) AFP=$OPTARG;;
        g) GPU=$OPTARG;;
        n) n=$OPTARG;;
    esac
done
echo $GPU

if [ "$LOCAL" != true ]; then
    module use /project/6004902/modulefiles
    module load presto
    module load chime-psr
    source ~/projects/rrg-istairs-ad/Your/bin/activate
    module load cuda
fi

#work in absolute paths, CC is weird when launching batch script
CAND_PATH=$(pwd)
echo $PWD $CAND_PATH
if [ "$LOCAL" != true ]; then
    ls
    cp -r ./* $SLURM_TMPDIR
    cd $SLURM_TMPDIR
fi
mkdir -p nsub_0_5
mkdir -p nsub_1
mkdir -p nsub_short_0_5
mkdir -p nsub_0_1
mkdir -p nsub_0_1_short
echo "making candidates"
# python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_5 -r -n $n -ws 500 --gpu_id $GPU
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_1 -r -n $n -ws 1000 --gpu_id $GPU
# python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_short_0_5 -r -n $n -ws 500 --gpu_id $GPU --range_dm 5
# python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_1 -r -n $n -ws 100 --gpu_id $GPU
# python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_1_short -r -n $n -ws 100 --gpu_id $GPU --range_dm 5


#make plots and do a predict for general pulses
echo "grading candidates"
if [ "$LOCAL" != true ]; then
    #go backwards, deactivate everything
    module unload cuda
    deactivate
    module unload chime-psr
    module unload presto
    module unuse /project/6004902/modulefiles
    source ~/projects/rrg-istairs-ad/GWG2/environments/AFP/bin/activate
else
    echo "running locally, make sure FETCH, Your and sigpyproc are installed!!"
    source ~/anaconda3/etc/profile.d/conda.sh
    echo "activating fetch"
    conda activate fetch
    echo "Fetch activated"
fi
echo "predicting candidates"
# predict.py --data_dir nsub_0_5 --model a --probability 0.1 -g $GPU
predict.py --data_dir nsub_1 --model a --probability 0.1 -g $GPU
# predict.py --data_dir nsub_short_0_5 --model a --probability 0.1 -g $GPU
# predict.py --data_dir nsub_0_1 --model a --probability 0.1 -g $GPU
# predict.py --data_dir nsub_0_1_short --model a --probability 0.1 -g $GPU

if [ "$LOCAL" != true ]; then
    cp -r $SLURM_TMPDIR/nsub* $CAND_PATH/
fi
