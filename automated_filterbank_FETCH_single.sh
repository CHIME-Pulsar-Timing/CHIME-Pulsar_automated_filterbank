#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=2:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=5
#SBATCH --job-name=fetch
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --gres=gpu:v100l:1
#the first argument is the tree to search down
#run FETCH
while getopts "li:p:" flag
do
    case "${flag}" in
        l) LOCAL=true;;
        i) CAND_PATH=$OPTARG;;
        p) AFP=$OPTARG;;
    esac
done


#the following code is only valid for Adam's personal computer comment out if on CC
# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate fetch
#the following is valid for CC
if [ "$LOCAL" != true ]; then
    module use /project/6004902/modulefiles
    module load presto
    module load chime-psr
    source ~/projects/rrg-istairs-ad/Your/bin/activate
    module load cuda
else
    PULSAR=$(basename $CAND_PATH)
fi

#work in absolute paths, CC is weird when launching batch script
echo $AFP
echo "copying files to tmpdir"
echo $PWD $CAND_PATH
ls

mkdir -p nsub_0_5
mkdir -p nsub_1
echo "making candidates"
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_5 -r -n 5 -ws 500 --gpu_id 0
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_1 -r -n 5 -ws 1000 --gpu_id 0
#make plots and do a predict for general pulses
echo "grading candidates"
if [ "$LOCAL" != true ]; then
    source ~/projects/rrg-istairs-ad/GWG2/environments/AFP/bin/activate
else
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate fetch
fi

predict.py --data_dir nsub_0_5 --model a --probability 0.1
predict.py --data_dir nsub_1 --model a --probability 0.1
#remove all the fil files
