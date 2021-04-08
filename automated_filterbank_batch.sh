#!/bin/bash
#SBATCH --account=def-istairs
#SBATCH --export=NONE
#SBATCH --time=10:00:00
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=1
#SBATCH --job-name=automated_filterbank
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#1 is number of splits
#2 is filterbank file
#3 is DM
#if we are not splitting then just set out as the fil file

#source ~/afp/bin/activate
#path to automated filterbank file script locations
afp=$4

if [ $1 -gt 1 ]
then
    OUT=`python $afp/split_filterbank.py $1 $2 ${SLURM_TMPDIR}`
else
    OUT=$2
fi
i=0
echo $out
for FIL in $OUT;
do
    echo $FIL
    SPFILES="${SLURM_TMPDIR}/$i"
    if [ !-d $SPFILES ]; then
        mkdir $SPFILES
    fi
    mv $FIL $SPFILES
    #run pipeline and prep_fetch prep spegID
    FIL="$SPFILES/$FIL"
    python $afp/gwg_cand_search_pipeline.py --dm $3 --speg --fetch --no_fft --rfifind --sk_mad --dedisp --sp --fil $FIL --slurm "${SLURM_TMPDIR}/$i"
    $DATA="$SPFILES/data"
    if [ !-d $DATA ]; then
        mkdir "$DATA"
    fi
    #remove all the large files that we no longer need
    rm "$SPFILES/*.dat"
    if [ $1 -gt 1 ]
    then
        #if we made the SK file then delete it too
        rm "$SPFILES/$FIL"
    fi

    #run FETCH
    #the following code is only valid for Adam's personal computer
    #source ~/anaconda3/etc/profile.d/conda.sh    
    #conda activate fetch

    candmaker.py --frequency_size 256 --time_size 256 --cand_param_file "$SPFILES/cands.csv" --plot --fout $DATA
    
    #don't do predict as we don't have GPU allocation... this can be done in seperate script
    #predict.py --data_dir data/ --model a
    #conda deactivate
    ((i=i+1))
done
#now copy all the files back
cp -r ${SLURM_TMPDIR}/* .
