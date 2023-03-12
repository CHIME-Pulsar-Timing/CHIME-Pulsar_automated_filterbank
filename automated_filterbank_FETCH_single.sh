#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=10:00:00
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
    SLURM_TMPDIR='/media/adam/d0fdb915-c69f-4fba-9759-ed1844c4685b/tmpdir/'$PULSAR
    mkdir -p $SLURM_TMPDIR
fi

#work in absolute paths, CC is weird when launching batch script
echo $AFP
CAND_PATH=$(pwd)
echo $PWD $CAND_PATH
ls
cp -r ./* $SLURM_TMPDIR
cd $SLURM_TMPDIR
mkdir -p nsub_0_5
mkdir -p nsub_1
mkdir -p nsub_short_0_5
mkdir -p nsub_0_1
mkdir -p nsub_0_1_short
echo "making candidates"
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_5 -r -n 5 -ws 500 --gpu_id 0
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_1 -r -n 5 -ws 1000 --gpu_id 0
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_short_0_5 -r -n 5 -ws 500 --gpu_id 0 --range_dm 5
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_1 -r -n 5 -ws 100 --gpu_id 0
python $AFP/your_candmaker.py -fs 256 -ts 256 -c cands.csv -o nsub_0_1_short -r -n 5 -ws 100 --gpu_id 0 --range_dm 5


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
    source ~/anaconda3/etc/profile.d/conda.sh
    conda activate fetch
fi

predict.py --data_dir nsub_0_5 --model a --probability 0.1
predict.py --data_dir nsub_1 --model a --probability 0.1
predict.py --data_dir nsub_short_0_5 --model a --probability 0.1
predict.py --data_dir nsub_0_1 --model a --probability 0.1
predict.py --data_dir nsub_0_1_short --model a --probability 0.1
cp -r $SLURM_TMPDIR/* $CAND_PATH/
